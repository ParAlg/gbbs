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

// Wrapper for a hash function
struct hashtup {
  inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
};

// For triangle peeling
template <class Graph, class Graph2, class F, class H>
void vtx_intersect(Graph& G, Graph2& DG, F f, H ignore_f, uintE vg, uintE vdg) {
  // Intersect the neighbors of vg (vertex in G) with the neighbros of vdg (vertex in DG)
  size_t v_deg = G.get_vertex(vg).getOutDegree();
  size_t i_deg = DG.get_vertex(vdg).getOutDegree();
  auto i_iter = DG.get_vertex(vdg).getOutIter(vdg);
  auto v_iter = G.get_vertex(vg).getOutIter(vg);
  size_t i_iter_idx = 0;
  size_t v_iter_idx = 0;
  // Count of number of neighbors in the intersection
  size_t j = 0;
  while (i_iter_idx < i_deg && v_iter_idx < v_deg) {
    if (std::get<0>(i_iter.cur()) == std::get<0>(v_iter.cur())) {
      if (ignore_f(vg, std::get<0>(i_iter.cur())) && ignore_f(vdg, std::get<0>(i_iter.cur()))) {
        // If the neighbor is valid (ignore_f), one triangle is added to the neighbor
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
  // j triangles are added to the neighbor vdg
  f(vdg, j);
}


template <class Graph, class Graph2, class F, class H, class I>
inline size_t cliqueUpdate(Graph& G, Graph2& DG, size_t k, size_t max_deg, bool label, F get_active, size_t active_size,
  size_t granularity, char* still_active, sequence<uintE> &rank, sequence<size_t>& per_processor_counts, H update,
  bool do_update_changed, I update_changed) {
  const size_t n = G.n;
  auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, n, label, true); };
  auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };

  size_t active_deg = 0;
  auto degree_map = pbbslib::make_sequence<size_t>(active_size, [&] (size_t i) { return G.get_vertex(get_active(i)).getOutDegree(); });
  active_deg += pbbslib::reduce_add(degree_map);

  parallel_for (0, active_size, [&] (size_t j) {still_active[get_active(j)] = 1;}, 2048);
  
  size_t edge_table_size = (size_t) (active_deg < n ? active_deg : n); 
  auto edge_table = sparse_table<uintE, bool, hashtup>(edge_table_size, std::make_tuple(UINT_E_MAX, false), hashtup());

  parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active_size, [&](size_t i, HybridSpace_lw* induced) {
    auto vert = get_active(i);
    if (G.get_vertex(vert).getOutDegree() != 0) {
      auto ignore_f = [&](const uintE& u, const uintE& v) {
        auto status_u = still_active[u]; auto status_v = still_active[v];
        if (status_u == 2 || status_v == 2) return false; /* deleted edge */
        if (status_v == 0) return true; /* higher edge */
        return rank[u] < rank[v]; /* orient edge within bucket */
      };
      induced->setup(G, DG, k, vert, ignore_f, still_active);
      auto update_d = [&](uintE vtx, size_t count) {
        size_t worker = worker_id();
        size_t ct = per_processor_counts[worker*n + vtx];
        per_processor_counts[worker*n + vtx] += count;
        if (ct == 0) edge_table.insert(std::make_tuple(vtx, true));
      };
      induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
    }
  }, granularity, false);

  /* extract the vertices that had their count changed */
  auto changed_vtxs = edge_table.entries();
  edge_table.del();

  /* Aggregate the updated counts across all worker's local arrays, into the
   * first worker's array. Also zero out the other worker's updated counts. */
  parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
    size_t nthreads = num_workers();
    uintE v = std::get<0>(changed_vtxs[i]);
    for (size_t j=0; j<nthreads; j++) {
      update(per_processor_counts, j, v);
    }
  }, 128);

  if (do_update_changed) {
    parallel_for(0, changed_vtxs.size(), [&] (size_t i) { update_changed(per_processor_counts, i, std::get<0>(changed_vtxs[i])); });
  }

  /* mark all as deleted */
  parallel_for (0, active_size, [&] (size_t j) {still_active[get_active(j)] = 2;}, 2048);

  return changed_vtxs.size();
}

template <class Graph, class Graph2, class F, class I>
inline size_t cliqueUpdate_serial(Graph& G, Graph2& DG, size_t k, size_t max_deg, bool label, F get_active, size_t active_size,
  char* still_active, sequence<uintE> &rank, sequence<size_t>& per_processor_counts, 
  bool do_update_changed, I update_changed, sequence<uintE>& update_idxs) {
  const size_t n = G.n;

  parallel_for (0, active_size, [&] (size_t j) {still_active[get_active(j)] = 1;}, 2048);

  size_t num_updates = 0;
  HybridSpace_lw* induced = new HybridSpace_lw();
  induced->alloc(max_deg, k, n, label, true);
  for (size_t i=0; i < active_size; i++) {
    auto vtx = get_active(i);
    if (G.get_vertex(vtx).getOutDegree() != 0) {
      auto ignore_f = [&](const uintE& u, const uintE& v) {
        auto status_u = still_active[u]; auto status_v = still_active[v];
        if (status_u == 2 || status_v == 2) return false; /* deleted edge */
        if (status_v == 0) return true; /* higher edge */
        return rank[u] < rank[v]; /* orient edge within bucket */
      };
      induced->setup(G, DG, k, vtx, ignore_f, still_active);
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

  if (do_update_changed) {
    for (size_t i=0; i < num_updates; i++) {
      update_changed(per_processor_counts, i, update_idxs[i]);
    }
  }

  /* mark all as deleted */
  parallel_for (0, active_size, [&] (size_t j) {still_active[get_active(j)] = 2;}, 2048);

  return num_updates;
}

template <typename bucket_t, class Graph, class Graph2>
sequence<bucket_t> Peel(Graph& G, Graph2& DG, size_t k, size_t* cliques, bool label, sequence<uintE> &rank,
  size_t num_buckets=16) {
  timer t2; t2.start();
  const size_t n = G.n;
  auto D = sequence<bucket_t>(n, [&](size_t i) { return cliques[i]; });
  auto D_filter = sequence<std::tuple<uintE, bucket_t>>(n);

  auto b = make_vertex_custom_buckets<bucket_t>(n, D, increasing, num_buckets);

  auto per_processor_counts = sequence<size_t>(n*num_workers(), static_cast<size_t>(0));

  char* still_active = (char*) calloc(n, sizeof(char));
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active?
  auto update_idxs = sequence<uintE>(max_deg);

  size_t rounds = 0;
  size_t finished = 0;
  bucket_t cur_bkt = 0;
  bucket_t max_bkt = 0;
  // Peel each bucket
  auto update_clique = [&](sequence<size_t>& ppc, size_t j, uintE v) {
    if (j == 0) return;
    ppc[v] += ppc[j*n + v];
    ppc[j*n + v] = 0;
  };
  while (finished != n) {
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    auto active = vertexSubset(n, bkt.identifiers);
    cur_bkt = bkt.id;
    
    finished += active.size();
    max_bkt = std::max(cur_bkt, max_bkt);

    auto get_active = [&](size_t j){ return active.vtx(j); };

    size_t granularity = (cur_bkt * active.size() < 10000) ? 1024 : 1;

    size_t filter_size = 0;

    if (active.size() > 1) {
      auto update_changed = [&](sequence<size_t>& ppc, size_t i, uintE v){
        /* Update the clique count for v, and zero out first worker's count */
        cliques[v] -= ppc[v];
        ppc[v] = 0;
        bucket_t deg = D[v];
        if (deg > cur_bkt) {
          bucket_t new_deg = std::max((bucket_t) cliques[v], (bucket_t) cur_bkt);
          D[v] = new_deg;
          // store (v, bkt) in an array now, pass it to apply_f below instead of what's there right now -- maybe just store it in D_filter?
          D_filter[i] = std::make_tuple(v, b.get_bucket(deg, new_deg));
        } else D_filter[i] = std::make_tuple(UINT_E_MAX, 0);
      };
      filter_size = cliqueUpdate(G, DG, k, max_deg, label, get_active, active.size(), granularity, still_active, rank, per_processor_counts, update_clique, true, update_changed);
    }
    else {
      size_t filter_size_serial = 0;
      auto update_changed_serial = [&](sequence<size_t>& ppc, size_t i, uintE v) {
        /* Update the clique count for v, and zero out first worker's count */
        cliques[v] -= ppc[v];
        ppc[v] = 0;
        bucket_t deg = D[v];
        if (deg > cur_bkt) {
          bucket_t new_deg = std::max((bucket_t) cliques[v], (bucket_t) cur_bkt);
          D[v] = new_deg;
          D_filter[filter_size_serial] = std::make_tuple(v, b.get_bucket(deg, new_deg));
          filter_size_serial++;
        }
      };
      cliqueUpdate_serial(G, DG, k, max_deg, label, get_active, active.size(), still_active, rank, per_processor_counts, true, update_changed_serial, update_idxs);
      filter_size = filter_size_serial;
    }

    auto apply_f = [&](size_t i) -> Maybe<std::tuple<uintE, bucket_t>> {
      uintE v = std::get<0>(D_filter[i]);
      bucket_t bucket = std::get<1>(D_filter[i]);
      if (v != UINT_E_MAX && still_active[v] != 2) return wrap(v, bucket);
      return Maybe<std::tuple<uintE, bucket_t> >();
    };
  
    b.update_buckets(apply_f, filter_size);

    active.del();

    rounds++;
  }

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  std::cout << "rho: " << rounds << std::endl;
  std::cout << "max_bkt: " << max_bkt << std::endl;

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
  char* still_active = (char*) calloc(n, sizeof(char));
  size_t max_deg = induced_hybrid::get_max_deg(G);
  auto per_processor_counts = sequence<size_t>(n*num_workers(), static_cast<size_t>(0));

  auto update_clique = [&](sequence<size_t>& ppc, size_t j, uintE v) {
    D[v] -= ppc[j*n + v];
    ppc[j*n + v] = 0;
  };
  auto nop = [&](sequence<size_t>& ppc, size_t i, uintE v) { return; };

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
    size_t granularity = (rho * active_size < 10000) ? 1024 : 1;
    auto get_active = [&](size_t j) { return this_arr[j]; };
    cliqueUpdate(G, DG, k, max_deg, label, get_active, active_size, granularity, still_active, rank, per_processor_counts, update_clique, false, nop);
    parallel_for (0, active_size, [&] (size_t j) {D[this_arr[j]] = 0;}, 2048);
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
    if (current_density > max_density) max_density = current_density;

    auto keep_seq = pbbs::delayed_seq<bool>(vtxs_remaining.size(), [&] (size_t i) {
      return !(D[vtxs_remaining[i]] <= target_density);
    });

    auto split_vtxs_m = pbbs::split_two(vtxs_remaining, keep_seq);
    uintE* this_arr = split_vtxs_m.first.to_array();
    size_t num_removed = split_vtxs_m.second;

    num_vertices_remaining -= num_removed;
    if (num_vertices_remaining > 0) {
    size_t active_size = num_removed;

// remove this_arr vertices ************************************************
    size_t granularity = (rho * active_size < 10000) ? 1024 : 1;
    auto get_active = [&](size_t j) { return this_arr[j]; };
    cliqueUpdate(G, DG, k, max_deg, label, get_active, active_size, granularity, still_active, rank, per_processor_counts, update_clique, false, nop);
    parallel_for (0, active_size, [&] (size_t j) {D[this_arr[j]] = 0;}, 2048);
  //***************
    }

    round++;
    pbbs::free_array(last_arr);
    last_arr = this_arr;
    remaining_offset = num_removed;
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



template <class Graph, class Graph2, class F, class H, class I>
inline size_t triUpdate(Graph& G, Graph2& DG, F get_active, size_t active_size, size_t granularity, char* still_active, 
  sequence<uintE> &rank, sequence<size_t>& per_processor_counts, H update, bool do_update_changed, I update_changed) {
  using W = typename Graph::weight_type;
  auto n = G.n;
  parallel_for (0, active_size, [&] (size_t j) {still_active[get_active(j)] = 1;}, 2048);

  auto edge_table = sparse_table<uintE, bool, hashtup>(G.n, std::make_tuple(UINT_E_MAX, false), hashtup());

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

  parallel_for (0, active_size, [&](size_t i) {
    auto j = get_active(i);
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
    for (size_t j=0; j<nthreads; j++) {
      update(per_processor_counts, j, v);
    }
  }, 128);

  if (do_update_changed) {
    parallel_for(0, changed_vtxs.size(), [&] (size_t i) { update_changed(per_processor_counts, i, std::get<0>(changed_vtxs[i])); });
  }

  /* mark all as deleted */
  parallel_for (0, active_size, [&] (size_t j) {still_active[get_active(j)] = 2;}, 2048);

  return changed_vtxs.size();
}



template <typename bucket_t, class Graph, class Graph2>
sequence<bucket_t> TriPeel(Graph& G, Graph2& DG, size_t* cliques, sequence<uintE> &rank, size_t num_buckets=16) {
  auto n = G.n;
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
  timer updct_t, bkt_t, filter_t;
  timer next_b; timer round_t;

  auto update_tri = [&](sequence<size_t>& ppc, size_t j, uintE v) {
    if (j == 0) return;
    ppc[v] += ppc[j*n + v];
    ppc[j*n + v] = 0;
  };
  // Peel each bucket
  while (finished != G.n) {
    round_t.start();
    // Retrieve next bucket
    auto bkt = b.next_bucket();

    auto active = vertexSubset(G.n, bkt.identifiers);
    cur_bkt = bkt.id;
    finished += active.size();
    max_bkt = std::max(cur_bkt, max_bkt);

    auto update_changed = [&](sequence<size_t>& ppc, size_t i, uintE v){
      /* Update the clique count for v, and zero out first worker's count */
      cliques[v] -= ppc[v];
      ppc[v] = 0;
      bucket_t deg = D[v];
      if (deg > cur_bkt) {
        bucket_t new_deg = std::max((bucket_t) cliques[v], (bucket_t) cur_bkt);
        D[v] = new_deg;
        // store (v, bkt) in an array now, pass it to apply_f below instead of what's there right now -- maybe just store it in D_filter?
        D_filter[i] = std::make_tuple(v, b.get_bucket(deg, new_deg));
      } else D_filter[i] = std::make_tuple(UINT_E_MAX, 0);
    };

    auto get_active = [&](size_t j){ return active.vtx(j); };
    size_t granularity = (cur_bkt * active.size() < 10000) ? 1024 : 1;

    size_t filter_size = triUpdate(G, DG, get_active, active.size(), granularity, still_active, rank, per_processor_counts, update_tri, true, update_changed);

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
  }

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "max_bkt: " << max_bkt << std::endl;

  b.del();
  free(still_active);

  return D;
}


template <class Graph, class Graph2>
double ApproxTriPeel(Graph& G, Graph2& DG, size_t* cliques, size_t num_cliques,
  sequence<uintE> &rank, double eps) {
  //using W = typename Graph::weight_type;
  timer t2; t2.start();
  const size_t n = G.n;
  auto D = sequence<size_t>(n, [&](size_t i) { return cliques[i]; });
  auto vertices_remaining = pbbs::delayed_seq<uintE>(n, [&] (size_t i) { return i; });

  size_t round = 1;
  uintE* last_arr = nullptr;
  size_t remaining_offset = 0;
  size_t num_vertices_remaining = n;
  double density_multiplier = 3*(1.+eps);

  double max_density = 0.0;
  char* still_active = (char*) calloc(G.n, sizeof(char));
  //size_t max_deg = induced_hybrid::get_max_deg(G);
  auto per_processor_counts = sequence<size_t>(n*num_workers(), static_cast<size_t>(0));
  auto update_tri = [&](sequence<size_t>& ppc, size_t j, uintE v) {
    D[v] -= ppc[j*n + v];
    ppc[j*n + v] = 0;
  };
  auto nop = [&](sequence<size_t>& ppc, size_t i, uintE v) { return; };

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
    size_t active_size = split_vtxs_m.second;

// remove this_arr vertices ************************************************
    size_t granularity = (rho * active_size < 10000) ? 1024 : 1;
    auto get_active = [&](size_t j) { return this_arr[j]; };
    triUpdate(G, DG, get_active, active_size, granularity, still_active, rank, per_processor_counts, update_tri, false, nop);

    /* mark all as deleted */
    parallel_for (0, active_size, [&] (size_t j) {D[this_arr[j]] = 0;}, 2048);
  //***************

    round++;
    last_arr = this_arr;
    remaining_offset = active_size;
    num_vertices_remaining -= active_size;
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
    size_t active_size = split_vtxs_m.second;
    num_vertices_remaining -= active_size;
    if (num_vertices_remaining > 0) {

// remove this_arr vertices ************************************************
    size_t granularity = (rho * active_size < 10000) ? 1024 : 1;
    auto get_active = [&](size_t j) { return this_arr[j]; };
    triUpdate(G, DG, get_active, active_size, granularity, still_active, rank, per_processor_counts, update_tri, false, nop);

    /* mark all as deleted */
    parallel_for (0, active_size, [&] (size_t j) {D[this_arr[j]] = 0;}, 2048);
  //***************
    }

    round++;
    pbbs::free_array(last_arr);
    last_arr = this_arr;
    remaining_offset = active_size;
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
