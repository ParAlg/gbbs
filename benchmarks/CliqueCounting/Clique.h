// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <math.h>

// Library dependencies
#include "ligra/bucket.h"
#include "ligra/edge_map_reduce.h"
#include "ligra/ligra.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "ligra/pbbslib/sparse_table.h"
#include "pbbslib/list_allocator.h"
#include "pbbslib/integer_sort.h"

// Ordering files
#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

// Clique files
#include "intersect.h"
#include "induced_intersection.h"
#include "induced_neighborhood.h"
#include "induced_hybrid.h"
#include "induced_split.h"
#include "relabel.h"

#define SIMD_STATE 4

template <class Graph>
inline uintE* rankNodes(Graph& G, size_t n) {
  uintE* r = pbbslib::new_array_no_init<uintE>(n);
  sequence<uintE> o(n);

  timer t;
  t.start();
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) { o[i] = i; });

  pbbs::integer_sort_inplace(o.slice(), [&] (size_t p) {
    return G.get_vertex(p).getOutDegree();
  });
//  pbbslib::sample_sort_inplace(o.slice(), [&](const uintE u, const uintE v) {
//    return G.get_vertex(u).getOutDegree() < G.get_vertex(v).getOutDegree();
//  });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { r[o[i]] = i; });
  t.stop();
  t.reportTotal("Rank time");
  return r;
}

template <class Graph>
pbbs::sequence<uintE> get_ordering(Graph& GA, long order_type, double epsilon = 0.1) {
  if (order_type == 0) return goodrichpszona_degen::DegeneracyOrder_intsort(GA, epsilon);
  else if (order_type == 1) return barenboimelkin_degen::DegeneracyOrder(GA, epsilon);
  else if (order_type == 2) {
    auto rank = sequence<uintE>(GA.n, [&](size_t i) { return i; });
    auto kcore = KCore(GA);
    auto get_core = [&](uintE& p) -> uintE { return kcore[p]; };
    pbbs::integer_sort_inplace(rank.slice(), get_core);
    return rank;
  }
  else if (order_type == 3) return pbbslib::make_sequence(rankNodes(GA, GA.n), GA.n);
  else if (order_type == 4) {
    auto rank = sequence<uintE>(GA.n, [&](size_t i) { return i; });
    return rank;
  }
}


struct hashtup {
inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
};
template <typename bucket_t, class Graph, class Graph2>
sequence<long> Peel(Graph& G, Graph2& DG, size_t k, long* cliques, bool label, sequence<uintE> &rank, bool par_serial, size_t num_buckets=16) {
auto stats = sequence<size_t>(G.n);
timer t2; t2.start();
  size_t n = G.n;
  const size_t eltsPerCacheLine = 64/sizeof(long);
  auto D = sequence<long>(G.n, [&](size_t i) { return cliques[i]; });
  //auto D_update = sequence<long>(eltsPerCacheLine*G.n);
  //parallel_for(0, G.n, [&](size_t j){D_update[eltsPerCacheLine*j] = 0;});
  auto D_filter = sequence<std::tuple<uintE, long>>(G.n);

  auto b = make_vertex_custom_buckets<bucket_t>(G.n, D, increasing, num_buckets);

  auto per_processor_counts = sequence<size_t>(n*num_workers(), static_cast<size_t>(0));

  char* still_active = (char*) calloc(G.n, sizeof(char));
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active?
  auto update_idxs = sequence<long>(max_deg);

  size_t rounds = 0;
  size_t finished = 0;
  long cur_bkt = 0;
  long max_bkt = 0;
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
      for (size_t i=1; i<nthreads; i++) {
        per_processor_counts[v] += per_processor_counts[i*n + v];
        per_processor_counts[i*n + v] = 0;
      }
    }, 128);

    filter_t.start();
    parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
      const uintE v = std::get<0>(changed_vtxs[i]);
      /* Update the clique count for v, and zero out first worker's count */
      cliques[v] -= per_processor_counts[v];
      per_processor_counts[v] = 0;
      uintE deg = D[v];
      if (deg > cur_bkt) {
        long new_deg = std::max(cliques[v], (long) cur_bkt);
        D[v] = new_deg;
        long bkt = b.get_bucket(deg, new_deg);
        // store (v, bkt) in an array now, pass it to apply_f below instead of what's there right now -- maybe just store it in D_filter?
        D_filter[i] = std::make_tuple(v, bkt);
      } else D_filter[i] = std::make_tuple(UINT_E_MAX, LONG_MAX);
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
      uintE deg = D[v];
      if (deg > cur_bkt) {
        long new_deg = std::max(cliques[v], (long) cur_bkt);
        D[v] = new_deg;
        long bkt = b.get_bucket(deg, new_deg);
        D_filter[filter_size] = std::make_tuple(v, bkt);
        filter_size++;
      }
    }
    filter_t.stop();
}

    /* mark all as deleted */
    parallel_for (0, active.size(), [&] (size_t j) {still_active[active.vtx(j)] = 2;}, 2048);

    auto apply_f = [&](size_t i) -> Maybe<std::tuple<uintE, uintE>> {
      uintE v = std::get<0>(D_filter[i]);
      uintE bkt = std::get<1>(D_filter[i]);
      if (v != UINT_E_MAX && still_active[v] != 2) return wrap(v, bkt);
      return Maybe<std::tuple<uintE, uintE> >();
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
    auto bkt = _Peel_serial(G, DG, k, cliques, label, rank, still_active);
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


// induced
// generated
// -i 0 (simple gbbs intersect), -i 2 (simd intersect), -i 1 (graph set inter)
// -o 0 (goodrich), 1 (barnboimelkin approx), 2 (barenboimelkin exact)

// todo approx work and do some kind of break in gen if too much
// TODO get rid of duplicates in edge lists????
template <class Graph>
inline size_t Clique(Graph& GA, size_t k, long order_type, double epsilon, long space_type, bool label, bool filter, bool use_base,
  long recursive_level, bool par_serial) {
  std::cout << "### Starting clique counting" << std::endl;
  const size_t eltsPerCacheLine = 64/sizeof(long);
  long* per_vert = use_base ? (long*) calloc(GA.n*num_workers(), sizeof(long)) : nullptr;

  using W = typename Graph::weight_type;
  assert (k >= 3);
  // TODO put in triangle counting here

  timer t_rank; t_rank.start();
  sequence<uintE> rank = get_ordering(GA, order_type, epsilon);
  double tt_rank = t_rank.stop();
  std::cout << "### Rank Running Time: " << tt_rank << std::endl;

  timer t_filter; t_filter.start();
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return (rank[u] < rank[v]) && GA.get_vertex(u).getOutDegree() >= k-1 && GA.get_vertex(v).getOutDegree() >= k-1;
  };
  auto DG = filter ? filter_graph(GA, pack_predicate) : relabel_graph(GA, rank.begin(), pack_predicate);
  double tt_filter = t_filter.stop();
  std::cout << "### Filter Graph Running Time: " << tt_filter << std::endl;

  // Done preprocessing
  timer t; t.start();
  size_t count = 0;
  auto use_f = [&](const uintE& src, const uintE& u) { return true; };

  if (!use_base) {
    auto base_f = [&](uintE vtx, size_t count) {};
  if (space_type == 2) {
    count = induced_intersection::CountCliques(DG, k-1, use_f, base_f, use_base);
  }
  else if (space_type == 3) {
    count = induced_neighborhood::CountCliques(DG, k-1);
  }
  else if (space_type == 5) {
    count = induced_hybrid::CountCliques(DG, k-1, base_f, use_base, label, recursive_level);
  }
  else if (space_type == 6) {
    count = induced_split::CountCliques(DG, k-1, base_f, use_base, label, recursive_level);
  }
  } else {
    auto base_f = [&](uintE vtx, size_t count) {
      //pbbslib::xadd(&(per_vert[eltsPerCacheLine*(vtx+worker_id()*GA.n)]), (long) count);
      per_vert[(vtx+worker_id()*GA.n)] += count;
    }; // TODO problem with relabel not being consistent; but if using filter should be ok
  if (space_type == 2) {
    count = induced_intersection::CountCliques(DG, k-1, use_f, base_f, use_base);
  }
  else if (space_type == 3) {
    count = induced_neighborhood::CountCliques(DG, k-1);
  }
  else if (space_type == 5) {
    count = induced_hybrid::CountCliques(DG, k-1, base_f, use_base, label, recursive_level);
  }
  else if (space_type == 6) {
    count = induced_split::CountCliques(DG, k-1, base_f, use_base, label, recursive_level);
  }
  }

  double tt = t.stop();
  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << k << " cliques = " << count << "\n";

  if (!use_base) return count;

  for (size_t j=1; j < num_workers(); j++) {
    parallel_for(0,GA.n,[&](size_t l) {
      per_vert[l] += per_vert[(l + j*GA.n)];
    });
  }


  long* inverse_per_vert = use_base && !filter ? (long*) malloc(GA.n*sizeof(long)) : nullptr;
  if (!filter) {
    parallel_for(0, GA.n, [&] (size_t i) { inverse_per_vert[i] = per_vert[rank[i]]; });
    free(per_vert);
    per_vert = inverse_per_vert;
  }

  auto per_vert_seq = pbbslib::make_sequence<long>(GA.n, [&] (size_t i) { return per_vert[i]; });
  auto max_per_vert = pbbslib::reduce_max(per_vert_seq);

  sequence<long> cores;
  std::cout << "Max per-vertex count is: " << max_per_vert << std::endl;
  if (max_per_vert >= std::numeric_limits<uintE>::max()) {
    std::cout << "Calling peeling with bucket_t = size_t (8-byte bucket types)" << std::endl;
    Peel<size_t>(GA, DG, k-1, per_vert, label, rank, par_serial);
  } else {
    std::cout << "Calling peeling with bucket_t = uintE (4-byte bucket types)" << std::endl;
    Peel<uintE>(GA, DG, k-1, per_vert, label, rank, par_serial);
  }

  free(per_vert);

  return count;
}


template <class Graph, class Graph2>
long _Peel_serial(Graph& G, Graph2& DG, size_t k, long* cliques, bool label, sequence<uintE> &rank, char* still_active) {
  bheapLLU* heap=mkheapLLU(cliques,still_active,G.n);
  auto heap_init_size = heap->n-1;
  auto D_update = sequence<long>(G.n);
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
  auto update_idxs = sequence<long>(max_deg);
  long c = 0;
  while (finished != heap_init_size) {
    num_updates = 0;
    auto kv=popminLLU(heap);
    auto v = kv.key;
    if (kv.value > c){
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