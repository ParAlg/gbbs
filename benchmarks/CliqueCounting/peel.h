#pragma once

#include "ligra/bucket.h"
#include "ligra/edge_map_reduce.h"
#include "ligra/ligra.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "ligra/pbbslib/sparse_table.h"
#include "pbbslib/list_allocator.h"
#include "pbbslib/integer_sort.h"

#include "induced_hybrid.h"
#include "intersect.h"

struct hashtup {
inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
};

template <class bucket_t, class Graph, class Graph2>
bucket_t _Peel_serial(Graph& G, Graph2& DG, size_t k, size_t* cliques, bool label, sequence<uintE> &rank, char* still_active) {
  bheapLLU* heap=mkheapLLU(cliques,still_active,G.n);
  auto heap_init_size = heap->n-1;
  auto D_update = sequence<bucket_t>(G.n);
  size_t num_updates = 0;
  size_t finished=0;
  size_t max_deg = induced_hybrid::get_max_deg(G);
  /*parallel_for(0, DG.n, [&] (size_t i) {
    if (still_active[i] != 2) {
      size_t deg = DG.get_vertex(i).getOutDegree();
      pbbs::write_min(&max_deg, deg, std::greater<size_t>());
    }
  });*/
  HybridSpace_lw* induced = new HybridSpace_lw();
  induced->alloc(max_deg, k, G.n, label, true);
  auto update_idxs = sequence<uintE>(max_deg);
  bucket_t c = 0;
  while (finished != heap_init_size) {
    num_updates = 0;
    auto kv=popminLLU(heap);
    auto v = kv.key;
    if (kv.value > static_cast<long>(c)){
			c = kv.value;
		}
    if (G.get_vertex(v).getOutDegree() != 0) {
      auto ignore_f = [&](const uintE& a, const uintE& b) {
        if (still_active[a] == 2 || still_active[b] == 2) return false;
        return true;
      };
      auto update_d = [&](uintE vtx, size_t count) {
        if (D_update[vtx] == 0 && count > 0) {
          update_idxs[num_updates] = vtx;
          num_updates++;
        }
        D_update[vtx] += count;
      };
      induced->setup(G, DG, k, v, ignore_f, still_active);
      induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
    }

    for (size_t j=0; j < num_updates; j++) {
      auto u = update_idxs[j];
      updateLLU(heap,u,D_update[u]);
      D_update[u] = 0;
    }

    still_active[v] = 2;
    finished++;
  }

  freeheapLLU(heap);
  if (induced != nullptr) { delete induced; }

  return c;

}



template <typename bucket_t, class Graph, class Graph2>
sequence<bucket_t> Peel(Graph& G, Graph2& DG, size_t k, size_t* cliques, bool label, sequence<uintE> &rank, bool par_serial, size_t num_buckets=16) {
auto stats = sequence<size_t>(G.n);
timer t2; t2.start();
  size_t n = G.n;
  //const size_t eltsPerCacheLine = 64/sizeof(long);
  auto D = sequence<bucket_t>(G.n, [&](size_t i) { return cliques[i]; });
  //auto D_update = sequence<long>(eltsPerCacheLine*G.n);
  //parallel_for(0, G.n, [&](size_t j){D_update[eltsPerCacheLine*j] = 0;});
  auto D_filter = sequence<std::tuple<uintE, bucket_t>>(G.n);

  auto b = make_vertex_custom_buckets<bucket_t>(G.n, D, increasing, num_buckets);

  auto per_processor_counts = sequence<size_t>(n*num_workers(), static_cast<size_t>(0));

  char* still_active = (char*) calloc(G.n, sizeof(char));
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active?
  auto update_idxs = sequence<uintE>(max_deg);

  size_t rounds = 0;
  size_t finished = 0;
  bucket_t cur_bkt = 0;
  bucket_t max_bkt = 0;
  timer updct_t, bkt_t, filter_t;
  timer next_b;
  bool log = false;
  // Peel each bucket
  while (finished != G.n) {
    timer round_t; round_t.start();
    // Retrieve next bucket
    next_b.start();
    auto bkt = b.next_bucket();
    next_b.stop();
    auto active = vertexSubset(G.n, bkt.identifiers);
    stats[rounds] = active.size();
    cur_bkt = bkt.id;
    if (log) {
      std::cout << "bkt = " << cur_bkt << " size = " << active.size() << std::endl;
    }
    if (par_serial && active.size() <= 200) { /* switch */
      break;
    }
    finished += active.size();
    max_bkt = std::max(cur_bkt, max_bkt);

    size_t active_deg = 0;
    auto degree_map = pbbslib::make_sequence<size_t>(active.size(), [&] (size_t i) { return G.get_vertex(active.vtx(i)).getOutDegree(); });
    active_deg += pbbslib::reduce_add(degree_map);

    parallel_for (0, active.size(), [&] (size_t j) {still_active[active.vtx(j)] = 1;}, 2048);
    // here, update D[i] if necessary
    // for each vert in active, just do the same kickoff, but we drop neighbors if they're earlier in the active set
    // also drop if already peeled -- check using D

    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, G.n, label, true); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };

    size_t granularity = (cur_bkt * active.size() < 10000) ? 1024 : 1;
    size_t filter_size = 0;
if (active.size() > 1) {
    size_t edge_table_size = std::min((size_t) cur_bkt*k*active.size(), (size_t) (active_deg < G.n ? active_deg : G.n));
    auto edge_table = sparse_table<uintE, bool, hashtup>(edge_table_size, std::make_tuple(UINT_E_MAX, false), hashtup());
    updct_t.start();
    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active.size(), [&](size_t i, HybridSpace_lw* induced) {
      if (G.get_vertex(active.vtx(i)).getOutDegree() != 0) {
        auto ignore_f = [&](const uintE& u, const uintE& v) {
          auto status_u = still_active[u]; auto status_v = still_active[v];
          if (status_u == 2 || status_v == 2) return false; /* deleted edge */
          if (status_v == 0) return true; /* higher edge */
          return rank[u] < rank[v]; /* orient edge within bucket */
        };
        induced->setup(G, DG, k, active.vtx(i), ignore_f, still_active);
        auto update_d = [&](uintE vtx, size_t count) {
          size_t worker = worker_id();
          size_t ct = per_processor_counts[worker*n + vtx];
          per_processor_counts[worker*n + vtx] += count;
          if (ct == 0) { /* only need bother trying hash table if our count is 0 */
            edge_table.insert(std::make_tuple(vtx, true));
          }
        };
        induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
      }
    }, granularity, false);
    updct_t.stop();

    /* extract the vertices that had their count changed */
    auto changed_vtxs = edge_table.entries();
    edge_table.del();

    /* Aggregate the updated counts across all worker's local arrays, into the
     * first worker's array. Also zero out the other worker's updated counts. */
    parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
      size_t nthreads = num_workers();
      uintE v = std::get<0>(changed_vtxs[i]);
      for (size_t j=1; j<nthreads; j++) {
        per_processor_counts[v] += per_processor_counts[j*n + v];
        per_processor_counts[j*n + v] = 0;
      }
    }, 128);

    filter_t.start();
    parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
      const uintE v = std::get<0>(changed_vtxs[i]);
      /* Update the clique count for v, and zero out first worker's count */
      cliques[v] -= per_processor_counts[v];
      per_processor_counts[v] = 0;
      bucket_t deg = D[v];
      if (deg > cur_bkt) {
        bucket_t new_deg = std::max((bucket_t) cliques[v], (bucket_t) cur_bkt);
        D[v] = new_deg;
        // store (v, bkt) in an array now, pass it to apply_f below instead of what's there right now -- maybe just store it in D_filter?
        D_filter[i] = std::make_tuple(v, b.get_bucket(deg, new_deg));
      } else D_filter[i] = std::make_tuple(UINT_E_MAX, 0);
    }, 2048);
    filter_t.stop();
    filter_size = changed_vtxs.size();
} else {
    updct_t.start();
    size_t num_updates = 0;
    HybridSpace_lw* induced = new HybridSpace_lw();
    induced->alloc(max_deg, k, G.n, label, true);
    for (size_t i=0; i < active.size(); i++) {
      if (G.get_vertex(active.vtx(i)).getOutDegree() != 0) {
        auto ignore_f = [&](const uintE& u, const uintE& v) {
          auto status_u = still_active[u]; auto status_v = still_active[v];
          if (status_u == 2 || status_v == 2) return false; /* deleted edge */
          if (status_v == 0) return true; /* higher edge */
          return rank[u] < rank[v]; /* orient edge within bucket */
        };
        induced->setup(G, DG, k, active.vtx(i), ignore_f, still_active);
        auto update_d = [&](uintE vtx, size_t count) {
          size_t ct = per_processor_counts[vtx];
          per_processor_counts[vtx] += count;
          if (ct == 0) { /* only need bother trying hash table if our count is 0 */
            update_idxs[num_updates] = vtx;
            num_updates++;
          }
        };
        induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
      }
    }
    if (induced != nullptr) { delete induced; }
    updct_t.stop();

    filter_t.start();
    filter_size = 0;
    for (size_t i=0; i < num_updates; i++) {
      const uintE v = update_idxs[i];
      /* Update the clique count for v, and zero out first worker's count */
      cliques[v] -= per_processor_counts[v];
      per_processor_counts[v] = 0;
      bucket_t deg = D[v];
      if (deg > cur_bkt) {
        bucket_t new_deg = std::max((bucket_t) cliques[v], (bucket_t) cur_bkt);
        D[v] = new_deg;
        D_filter[filter_size] = std::make_tuple(v, b.get_bucket(deg, new_deg));
        filter_size++;
      }
    }
    filter_t.stop();
}

    /* mark all as deleted */
    parallel_for (0, active.size(), [&] (size_t j) {still_active[active.vtx(j)] = 2;}, 2048);

    auto apply_f = [&](size_t i) -> Maybe<std::tuple<uintE, bucket_t>> {
      uintE v = std::get<0>(D_filter[i]);
      bucket_t bucket = std::get<1>(D_filter[i]);
      if (v != UINT_E_MAX && still_active[v] != 2) return wrap(v, bucket);
      return Maybe<std::tuple<uintE, bucket_t> >();
    };
    bkt_t.start();
    b.update_buckets(apply_f, filter_size);
    bkt_t.stop();

    active.del();

    rounds++;
    round_t.stop();
    if (log) {
      round_t.reportTotal("round time");
    }
  }

timer ser_t; ser_t.start();
  if (finished != G.n) {
    auto bkt = _Peel_serial<bucket_t>(G, DG, k, cliques, label, rank, still_active);
    max_bkt = std::max(max_bkt, bkt);
  }
ser_t.stop();
double tt2 = t2.stop();
std::cout << "### Peel Running Time: " << tt2 << std::endl;
ser_t.reportTotal("serial time");
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "max_bkt: " << max_bkt << std::endl;

  bkt_t.reportTotal("bkt time");
  next_b.reportTotal("next bucket time");
  filter_t.reportTotal("filter time");
  updct_t.reportTotal("update count time");

auto const Q1 = rounds / 4;
auto const Q2 = rounds / 2;
auto const Q3 = Q1 + Q2;
pbbs::integer_sort_inplace(stats.slice(0, rounds), [&] (size_t x) {
    return stats[x];
  });
std::cout << "Q1: " << stats[Q1] << std::endl;
std::cout << "Q2: " << stats[Q2] << std::endl;
std::cout << "Q3: " << stats[Q3] << std::endl;

  b.del();
  free(still_active);

  return D;
}


/*template <typename bucket_t>
sequence<bucket_t> TriPeel(symmetric_graph<csv_bytepd_amortized, pbbs::empty>& G, symmetric_graph<csv_byte, pbbs::empty>& DG, size_t* cliques, sequence<uintE>& rank){
  assert(false);
  return sequence<bucket_t>(G.n);
}*/



template <class Graph, class Graph2, class F, class H>
void vtx_intersect(Graph& G, Graph2& DG, F f, H ignore_f, uintE vg, uintE vdg) {
    size_t v_deg = G.get_vertex(vg).getOutDegree();
    size_t i_deg = DG.get_vertex(vdg).getOutDegree();
    auto i_iter = DG.get_vertex(vdg).getOutIter(vdg);
    auto v_iter = G.get_vertex(vg).getOutIter(vg);
    size_t i_iter_idx = 0;
    size_t v_iter_idx = 0;
    size_t j = 0;
    while (i_iter_idx < i_deg && v_iter_idx < v_deg) {
      if (std::get<0>(i_iter.cur()) == std::get<0>(v_iter.cur())) {
        if (ignore_f(vg, std::get<0>(i_iter.cur())) && ignore_f(vdg, std::get<0>(i_iter.cur()))) {
          f(std::get<0>(i_iter.cur()), 1);
          j++;
        }
        i_iter_idx++; v_iter_idx++;
        if (i_iter.has_next()) i_iter.next();
        if (v_iter.has_next()) v_iter.next();
      } else if (std::get<0>(i_iter.cur()) < std::get<0>(v_iter.cur())) {
        i_iter_idx++;
        if (i_iter.has_next()) i_iter.next();
      }
      else {
        v_iter_idx++;
        if (v_iter.has_next()) v_iter.next();
      }
    }
    f(vdg, j);
}


template <typename bucket_t, class Graph, class Graph2>
sequence<bucket_t> TriPeel(Graph& G, Graph2& DG, size_t* cliques, sequence<uintE> &rank, size_t num_buckets=16) {
  using W = typename Graph::weight_type;
  auto n = G.n;
//size_t k = 2;
auto stats = sequence<size_t>(G.n);
timer t2; t2.start();
  auto D = sequence<bucket_t>(G.n, [&](size_t i) { return cliques[i]; });
  auto D_filter = sequence<std::tuple<uintE, bucket_t>>(G.n);
  auto b = make_vertex_custom_buckets<bucket_t>(G.n, D, increasing, num_buckets);
  auto per_processor_counts = sequence<size_t>(n*num_workers(), static_cast<size_t>(0));
  char* still_active = (char*) calloc(G.n, sizeof(char));

  size_t rounds = 0;
  size_t finished = 0;
  bucket_t cur_bkt = 0;
  bucket_t max_bkt = 0;
  bool log = false;
  timer updct_t, bkt_t, filter_t;
  timer next_b; timer round_t;
  // Peel each bucket
  while (finished != G.n) {
    round_t.start();
    // Retrieve next bucket
    auto bkt = b.next_bucket();

    auto active = vertexSubset(G.n, bkt.identifiers);
    stats[rounds] = active.size();
    cur_bkt = bkt.id;
    finished += active.size();
    max_bkt = std::max(cur_bkt, max_bkt);

    size_t active_deg = 0;
    auto degree_map = pbbslib::make_sequence<size_t>(active.size(), [&] (size_t i) { return G.get_vertex(active.vtx(i)).getOutDegree(); });
    active_deg += pbbslib::reduce_add(degree_map);

    parallel_for (0, active.size(), [&] (size_t j) {still_active[active.vtx(j)] = 1;}, 2048);

    size_t granularity = (cur_bkt * active.size() < 10000) ? 1024 : 1;
    size_t filter_size = 0;

    size_t edge_table_size = G.n; //std::min((size_t) cur_bkt*k*active.size(), (size_t) (active_deg < G.n ? active_deg : G.n))
    auto edge_table = sparse_table<uintE, bool, hashtup>(edge_table_size, std::make_tuple(UINT_E_MAX, false), hashtup());

    auto ignore_f = [&](const uintE& u, const uintE& v) {
      auto status_u = still_active[u]; auto status_v = still_active[v];
      if (status_u == 2 || status_v == 2) return false; /* deleted edge */
      if (status_v == 0) return true; /* higher edge */
      return rank[u] < rank[v]; /* orient edge within bucket */
    };
    auto update_d = [&](uintE vtx, size_t count) {
      size_t worker = worker_id();
      size_t ct = per_processor_counts[worker*n + vtx];
      per_processor_counts[worker*n + vtx] += count;
      if (ct == 0) edge_table.insert(std::make_tuple(vtx, true));
    };

    parallel_for (0, active.size(), [&](size_t i) {
      auto j = active.vtx(i);
      auto map_label_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        if (ignore_f(u, v)) vtx_intersect(G, DG, update_d, ignore_f, u, v);
      };
      G.get_vertex(j).mapOutNgh(j, map_label_f, false);
      // we want to intersect vtx's neighbors minus !ignore_f with v's out neighbors // TODO TODO TODO
    }, granularity, false);

    /* extract the vertices that had their count changed */
    auto changed_vtxs = edge_table.entries();
    edge_table.del();

    /* Aggregate the updated counts across all worker's local arrays, into the
     * first worker's array. Also zero out the other worker's updated counts. */
    parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
      size_t nthreads = num_workers();
      uintE v = std::get<0>(changed_vtxs[i]);
      for (size_t j=1; j<nthreads; j++) {
        per_processor_counts[v] += per_processor_counts[j*n + v];
        per_processor_counts[j*n + v] = 0;
      }
    }, 128);

    parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
      const uintE v = std::get<0>(changed_vtxs[i]);
      /* Update the clique count for v, and zero out first worker's count */
      cliques[v] -= per_processor_counts[v];
      per_processor_counts[v] = 0;
      bucket_t deg = D[v];
      if (deg > cur_bkt) {
        bucket_t new_deg = std::max((bucket_t) cliques[v], (bucket_t) cur_bkt);
        D[v] = new_deg;
        // store (v, bkt) in an array now, pass it to apply_f below instead of what's there right now -- maybe just store it in D_filter?
        D_filter[i] = std::make_tuple(v, b.get_bucket(deg, new_deg));
      } else D_filter[i] = std::make_tuple(UINT_E_MAX, 0);
    }, 2048);
    filter_size = changed_vtxs.size();

    /* mark all as deleted */
    parallel_for (0, active.size(), [&] (size_t j) {still_active[active.vtx(j)] = 2;}, 2048);

    auto apply_f = [&](size_t i) -> Maybe<std::tuple<uintE, bucket_t>> {
      uintE v = std::get<0>(D_filter[i]);
      bucket_t bucket = std::get<1>(D_filter[i]);
      if (v != UINT_E_MAX && still_active[v] != 2) return wrap(v, bucket);
      return Maybe<std::tuple<uintE, bucket_t> >();
    };
    b.update_buckets(apply_f, filter_size);

    active.del();

    rounds++;
    round_t.stop();
    if (log) {
      round_t.reportTotal("round time");
    }
  }

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "max_bkt: " << max_bkt << std::endl;

  bkt_t.reportTotal("bkt time");
  next_b.reportTotal("next bucket time");
  filter_t.reportTotal("filter time");
  updct_t.reportTotal("update count time");

//auto const Q1 = rounds / 4;
//auto const Q2 = rounds / 2;
//auto const Q3 = Q1 + Q2;
//pbbs::integer_sort_inplace(stats.slice(0, rounds), [&] (size_t x) {
//    return stats[x];
//  });
//std::cout << "Q1: " << stats[Q1] << std::endl;
//std::cout << "Q2: " << stats[Q2] << std::endl;
//std::cout << "Q3: " << stats[Q3] << std::endl;

  b.del();
  free(still_active);

  return D;
}





template <class Graph, class Graph2>
double ApproxPeel(Graph& G, Graph2& DG, size_t k, size_t* cliques, size_t num_cliques,
  bool label, sequence<uintE> &rank, double eps) {
    std::cout << "eps: " << eps << "\n";
    timer t2; t2.start();
  const size_t n = G.n;
  auto D = sequence<size_t>(n, [&](size_t i) { return cliques[i]; });
  auto vertices_remaining = pbbs::delayed_seq<uintE>(n, [&] (size_t i) { return i; });

  size_t round = 1;
  uintE* last_arr = nullptr;
  size_t remaining_offset = 0;
  size_t num_vertices_remaining = n;
  double density_multiplier = (k+1)*(1.+eps);

  double max_density = 0.0;
  char* still_active = (char*) calloc(G.n, sizeof(char));
  size_t max_deg = induced_hybrid::get_max_deg(G);
  auto per_processor_counts = sequence<size_t>(n*num_workers(), static_cast<size_t>(0));

  // First round
  {
    size_t edges_remaining = num_cliques;
    // Update density
    double current_density = ((double)edges_remaining) / ((double) (n));
    double target_density = (density_multiplier*((double)edges_remaining)) / ((double)vertices_remaining.size());
    auto rho = target_density;
    if (current_density > max_density) max_density = current_density;

    auto keep_seq = pbbs::delayed_seq<bool>(n, [&] (size_t i) {
      return !(D[i] <= target_density);
    });

    auto split_vtxs_m = pbbs::split_two(vertices_remaining, keep_seq);
    uintE* this_arr = split_vtxs_m.first.to_array();
    size_t num_removed = split_vtxs_m.second;
    size_t active_size = num_removed;

// remove this_arr vertices ************************************************
    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, G.n, label, true); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };

    size_t granularity = (rho * active_size < 10000) ? 1024 : 1;

    size_t active_deg = 0;
    auto degree_map = pbbslib::make_sequence<size_t>(active_size, [&] (size_t i) { return G.get_vertex(this_arr[i]).getOutDegree(); });
    active_deg += pbbslib::reduce_add(degree_map);

    parallel_for (0, active_size, [&] (size_t j) {still_active[this_arr[j]] = 1;}, 2048);
  
    size_t edge_table_size = (size_t) (active_deg < G.n ? active_deg : G.n); //std::min((size_t) rho*k*active_size, 
    auto edge_table = sparse_table<uintE, bool, hashtup>(edge_table_size, std::make_tuple(UINT_E_MAX, false), hashtup());
    //sequence<size_t> tots = sequence<size_t>(n, [&](size_t i) { return 0; });

    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active_size, [&](size_t i, HybridSpace_lw* induced) {
      if (G.get_vertex(this_arr[i]).getOutDegree() != 0) {
        auto ignore_f = [&](const uintE& u, const uintE& v) {
          auto status_u = still_active[u]; auto status_v = still_active[v];
          if (status_u == 2 || status_v == 2) return false; /* deleted edge */
          if (status_v == 0) return true; /* higher edge */
          return rank[u] < rank[v]; /* orient edge within bucket */
        };
        induced->setup(G, DG, k, this_arr[i], ignore_f, still_active);
        auto update_d = [&](uintE vtx, size_t count) {
          size_t worker = worker_id();
          size_t ct = per_processor_counts[worker*n + vtx];
          per_processor_counts[worker*n + vtx] += count;
          if (ct == 0) edge_table.insert(std::make_tuple(vtx, true));
        };
        induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
      }
    }, granularity, false);
    //edges_remaining -= pbbslib::reduce_add(tots);

    /* extract the vertices that had their count changed */
    auto changed_vtxs = edge_table.entries();
    edge_table.del();

    /* Aggregate the updated counts across all worker's local arrays, into the
     * first worker's array. Also zero out the other worker's updated counts. */
    parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
      size_t nthreads = num_workers();
      uintE v = std::get<0>(changed_vtxs[i]);
      for (size_t j=0; j<nthreads; j++) {
        D[v] -= per_processor_counts[j*n + v];
        per_processor_counts[j*n + v] = 0;
      }
    }, 128);

    /* mark all as deleted */
    parallel_for (0, active_size, [&] (size_t j) {D[this_arr[j]] = 0; still_active[this_arr[j]] = 2;}, 2048);
  //***************

    round++;
    last_arr = this_arr;
    remaining_offset = num_removed;
    num_vertices_remaining -= num_removed;
  }

  while (num_vertices_remaining > 0) {
    uintE* start = last_arr + remaining_offset;
    uintE* end = start + num_vertices_remaining;
    auto vtxs_remaining = pbbs::make_range(start, end);

    auto degree_f = [&] (size_t i) {
      uintE v = vtxs_remaining[i];
      return static_cast<size_t>(D[v]);
    };
    auto degree_seq = pbbslib::make_sequence<size_t>(vtxs_remaining.size(), degree_f);
    long edges_remaining = pbbslib::reduce_add(degree_seq);

    // Update density
    double current_density = ((double)edges_remaining) / ((double)vtxs_remaining.size());
    double target_density = (density_multiplier*((double)edges_remaining)) / ((double)vtxs_remaining.size());
    auto rho = target_density;
    //debug(std::cout << "Target density on round " << round << " is " << target_density << " erm = " << edges_remaining << " vrm = " << vtxs_remaining.size() << std::endl;
    //std::cout << "Current density on round " << round << " is " << current_density << std::endl;);
    if (current_density > max_density) max_density = current_density;

    auto keep_seq = pbbs::delayed_seq<bool>(vtxs_remaining.size(), [&] (size_t i) {
      return !(D[vtxs_remaining[i]] <= target_density);
    });

    auto split_vtxs_m = pbbs::split_two(vtxs_remaining, keep_seq);
    uintE* this_arr = split_vtxs_m.first.to_array();
    size_t num_removed = split_vtxs_m.second;
    //auto vs = vertexSubset(n, num_removed, this_arr);
    //debug(std::cout << "removing " << num_removed << " vertices" << std::endl;);

    num_vertices_remaining -= num_removed;
    if (num_vertices_remaining > 0) {
    size_t active_size = num_removed;

// remove this_arr vertices ************************************************
    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, G.n, label, true); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };

    size_t granularity = (rho * active_size < 10000) ? 1024 : 1;

    size_t active_deg = 0;
    auto degree_map = pbbslib::make_sequence<size_t>(active_size, [&] (size_t i) { return G.get_vertex(this_arr[i]).getOutDegree(); });
    active_deg += pbbslib::reduce_add(degree_map);

    parallel_for (0, active_size, [&] (size_t j) {still_active[this_arr[j]] = 1;}, 2048);
  
    size_t edge_table_size = (size_t) (active_deg < G.n ? active_deg : G.n); //std::min((size_t) rho*k*active_size, 
    auto edge_table = sparse_table<uintE, bool, hashtup>(edge_table_size, std::make_tuple(UINT_E_MAX, false), hashtup());
    //sequence<size_t> tots = sequence<size_t>(n, [&](size_t i) { return 0; });

    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active_size, [&](size_t i, HybridSpace_lw* induced) {
      if (G.get_vertex(this_arr[i]).getOutDegree() != 0) {
        auto ignore_f = [&](const uintE& u, const uintE& v) {
          auto status_u = still_active[u]; auto status_v = still_active[v];
          if (status_u == 2 || status_v == 2) return false; /* deleted edge */
          if (status_v == 0) return true; /* higher edge */
          return rank[u] < rank[v]; /* orient edge within bucket */
        };
        induced->setup(G, DG, k, this_arr[i], ignore_f, still_active);
        auto update_d = [&](uintE vtx, size_t count) {
          size_t worker = worker_id();
          size_t ct = per_processor_counts[worker*n + vtx];
          per_processor_counts[worker*n + vtx] += count;
          if (ct == 0) edge_table.insert(std::make_tuple(vtx, true));
        };
        induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
      }
    }, granularity, false);
    //edges_remaining -= pbbslib::reduce_add(tots);

    /* extract the vertices that had their count changed */
    auto changed_vtxs = edge_table.entries();
    edge_table.del();

    /* Aggregate the updated counts across all worker's local arrays, into the
     * first worker's array. Also zero out the other worker's updated counts. */
    parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
      size_t nthreads = num_workers();
      uintE v = std::get<0>(changed_vtxs[i]);
      for (size_t j=0; j<nthreads; j++) {
        D[v] -= per_processor_counts[j*n + v];
        per_processor_counts[j*n + v] = 0;
      }
    }, 128);

    /* mark all as deleted */
    parallel_for (0, active_size, [&] (size_t j) {D[this_arr[j]] = 0; still_active[this_arr[j]] = 2;}, 2048);
  //***************
    }

    round++;
    pbbs::free_array(last_arr);
    last_arr = this_arr;
    remaining_offset = num_removed;
    //if (vs.dense()) {
    //  pbbs::free_array(vs.d);
   // }
  }

  if (last_arr) {
    pbbs::free_array(last_arr);
  }
 
double tt2 = t2.stop();
std::cout << "### Peel Running Time: " << tt2 << std::endl;
std::cout << "rho: " << round << std::endl;
 std::cout << "### Density of (2(1+\eps))-Densest Subgraph is: " << max_density << endl;

  free(still_active);

  return max_density;
}
