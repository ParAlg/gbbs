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
#include "ligra/pbbslib/sparse_additive_map.h"
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


// induced
// generated
// -i 0 (simple gbbs intersect), -i 2 (simd intersect), -i 1 (graph set inter)
// -o 0 (goodrich), 1 (barnboimelkin approx), 2 (barenboimelkin exact)

// todo approx work and do some kind of break in gen if too much
// TODO get rid of duplicates in edge lists????
template <class Graph>
inline size_t Clique(Graph& GA, size_t k, long order_type, double epsilon, long space_type, bool label, bool filter, bool use_base) {
  std::cout << "### Starting clique counting" << std::endl;
  const size_t eltsPerCacheLine = 64/sizeof(long);
  long* per_vert = use_base ? (long*) calloc(eltsPerCacheLine*GA.n*num_workers(), sizeof(long)) : nullptr;

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
    count = induced_hybrid::CountCliques(DG, k-1, base_f, use_base, label);
  }
  else if (space_type == 6) {
    count = induced_split::CountCliques(DG, k-1, base_f, use_base, label);
  }
  } else {
    auto base_f = [&](uintE vtx, size_t count) {
      //pbbslib::xadd(&(per_vert[eltsPerCacheLine*(vtx+worker_id()*GA.n)]), (long) count);
      per_vert[eltsPerCacheLine*(vtx+worker_id()*GA.n)] += count;
    }; // TODO problem with relabel not being consistent; but if using filter should be ok
  if (space_type == 2) {
    count = induced_intersection::CountCliques(DG, k-1, use_f, base_f, use_base);
  }
  else if (space_type == 3) {
    count = induced_neighborhood::CountCliques(DG, k-1);
  }
  else if (space_type == 5) {
    count = induced_hybrid::CountCliques(DG, k-1, base_f, use_base, label);
  }
  else if (space_type == 6) {
    count = induced_split::CountCliques(DG, k-1, base_f, use_base, label);
  }
  }

  double tt = t.stop();
  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << k << " cliques = " << count << "\n";

  if (!use_base) return count;

  for (size_t j=1; j < num_workers(); j++) {
    parallel_for(0,GA.n,[&](size_t l) {
      per_vert[eltsPerCacheLine*l] += per_vert[eltsPerCacheLine*(l + j*GA.n)];
    });
  }

  timer t2; t2.start();
  long* inverse_per_vert = use_base && !filter ? (long*) malloc(eltsPerCacheLine*GA.n*sizeof(long)) : nullptr;
  if (!filter) {
    parallel_for(0, GA.n, [&] (size_t i) { inverse_per_vert[eltsPerCacheLine*i] = per_vert[eltsPerCacheLine*rank[i]]; });
    free(per_vert);
    per_vert = inverse_per_vert;
  }
  sequence<long> cores = Peel(GA, k-1, per_vert, label, rank);
  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;
  free(per_vert);

  return count;
}

struct hashtup {
inline size_t operator () (const uintE & a) {return pbbs::hash64_2(a);}
};
template <class Graph>
sequence<long> Peel(Graph& G, size_t k, long* cliques, bool label, sequence<uintE> &rank, size_t num_buckets=16) {
  using W = typename Graph::weight_type;
  const size_t eltsPerCacheLine = 64/sizeof(long);
  auto D = sequence<long>(G.n, [&](size_t i) { return cliques[eltsPerCacheLine*i]; });
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active?

  auto D_update = sequence<long>(num_workers()*G.n);
  parallel_for(0, num_workers()*G.n, [&](size_t j){D_update[j] = 0;});
  auto used_vert = sequence<uintE>(num_workers()*max_deg);
  auto used_vert_size = sequence<size_t>(num_workers()+1);
  parallel_for(0, num_workers()+1, [&](size_t j){used_vert_size[j] = 0;});

  auto D_filter = sequence<std::tuple<uintE, long>>(std::min(G.n, (size_t) 500));
  auto b = make_vertex_buckets(G.n, D, increasing, num_buckets);

  char* still_active = (char*) calloc(G.n, sizeof(char));
  

  size_t rounds = 0;
  size_t finished = 0;
  long cur_bkt = 0;
  long max_bkt = 0;
  timer updct_t, bkt_t, filter_t;
  // Peel each bucket
  while (finished != G.n) {
    // Retrieve next bucket
    bkt_t.start();
    auto bkt = b.next_bucket();
    bkt_t.stop();
    auto active = vertexSubset(G.n, bkt.identifiers);
    finished += active.size();
    cur_bkt = bkt.id;
    max_bkt = std::max(cur_bkt, (long) bkt.id);

    //size_t active_deg = 0;
    //for (size_t i=0; i < active.size(); i++) { active_deg += G.get_vertex(active.vtx(i)).getOutDegree(); }
    //size_t edge_table_size = std::min((size_t) cur_bkt*k*active.size(), (size_t) (active_deg < G.n ? active_deg : G.n));
    //auto edge_table = sparse_table<uintE, bool, hashtup>(edge_table_size, std::make_tuple(UINT_E_MAX, false), hashtup());

    parallel_for (0, active.size(), [&] (size_t j) {still_active[active.vtx(j)] = 1;});
    // here, update D[i] if necessary
    // for each vert in active, just do the same kickoff, but we drop neighbors if they're earlier in the active set
    // also drop if already peeled -- check using D

    updct_t.start();
    auto update_d = [&](uintE vtx, size_t count) {
      auto worker_id2 = worker_id();
      if (D_update[vtx+worker_id2*G.n] == 0) {
        used_vert[used_vert_size[worker_id2] + worker_id2*max_deg] = vtx;
        used_vert_size[worker_id2]++;
      }
      D_update[vtx+worker_id2*G.n] += count;
    };
    
    /*auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, G.n, label, true); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };
    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active.size(), [&](size_t i, HybridSpace_lw* induced) {
      if (G.get_vertex(active.vtx(i)).getOutDegree() != 0) {
        auto ignore_f = [&](const uintE& u, const uintE& v) {
          if (still_active[u] == 2 || still_active[v] == 2) return false;
          if (still_active[u] == 1 && still_active[v] == 0) return true;
          if (still_active[u] == 0 && still_active[v] == 1) return false;
          return rank[u] < rank[v];
        };
        induced->setup(G, k, active.vtx(i), ignore_f);
        induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
      }
    }, 1, false);*/
    auto init_induced = [&](SimpleSpace* induced) { induced->alloc(max_deg, k, G.n); };
    auto finish_induced = [&](SimpleSpace* induced) { if (induced != nullptr) { delete induced; } };
    parallel_for_alloc<SimpleSpace>(init_induced, finish_induced, 0, active.size(), [&](size_t l, SimpleSpace* induced) {
      auto i = active.vtx(l);
      if (G.get_vertex(i).getOutDegree() != 0) {
        auto use_f = [&](const uintE& u, const uintE& v) {
          if (still_active[u] == 2 || still_active[v] == 2) return false;
          if (still_active[u] == 1 && still_active[v] == 0) return true;
          if (still_active[u] == 0 && still_active[v] == 1) return false;
          return rank[u] < rank[v];
        };
        size_t j = 0;
        auto map_intersect_f = [&] (const uintE& src, const uintE& nbhr, const W& wgh) {
          if (use_f(i, nbhr)) {
            induced->induced[j] = nbhr;
            j++;
          }
        };
        G.get_vertex(i).mapOutNgh(i, map_intersect_f, false);
        induced->num_induced[0] = (uintE) j;
        if (j > 0) induced_intersection::KCliqueDir_simple(G, 1, k, induced, update_d, true);
      }
    } );
    updct_t.stop();

    parallel_for (0, active.size(), [&] (size_t j) {still_active[active.vtx(j)] = 2;});

    size_t filter_size = pbbslib::scan_add_inplace(used_vert_size.slice());

if (filter_size < 500) {

    filter_t.start();
    for (size_t i=0; i < num_workers(); i++) {
      for (size_t j=0; j < used_vert_size[i+1] - used_vert_size[i]; j++) {
        const uintE v = used_vert[j + i*max_deg];
        auto update_val = D_update[v + i*G.n];
        D_update[v + i*G.n] = 0;
      
      cliques[eltsPerCacheLine*v] -= update_val;
      uintE deg = D[v];
      if (deg > cur_bkt) {
        long new_deg = std::max(cliques[eltsPerCacheLine*v], (long) cur_bkt);
        D[v] = new_deg;
        long bkt = b.get_bucket(deg, new_deg);
        D_filter[used_vert_size[i] + j] = std::make_tuple(v, bkt);
      } else D_filter[used_vert_size[i] + j] = std::make_tuple(UINT_E_MAX, LONG_MAX);
      }
    }
    parallel_for (0, num_workers()+1, [&](size_t j) {used_vert_size[j] = 0;});

    filter_t.stop();

    auto apply_f = [&](size_t i) -> Maybe<std::tuple<uintE, uintE>> {
      uintE v = std::get<0>(D_filter[i]);
      uintE bkt = std::get<1>(D_filter[i]);
      if (v != UINT_E_MAX && still_active[v] != 2) return wrap(v, bkt);
      return Maybe<std::tuple<uintE, uintE> >();
    };
    bkt_t.start();
    b.update_buckets(apply_f, filter_size);
    bkt_t.stop();
} else {
      auto edge_table = sparse_additive_map<uintE, long>(filter_size, std::make_tuple(UINT_E_MAX, LONG_MAX));

    filter_t.start();
    parallel_for(0, num_workers(), [&] (size_t i) {
      for (size_t j=0; j < used_vert_size[i+1] - used_vert_size[i]; j++) {
        const uintE v = used_vert[j + i*max_deg];
        auto update_val = D_update[v + i*G.n];
        D_update[v + i*G.n] = 0;
        edge_table.insert(std::make_tuple(v, update_val));
      }
    });
    parallel_for (0, num_workers()+1, [&](size_t j) {used_vert_size[j] = 0;});
    auto edge_table_entries = edge_table.entries();
    edge_table.del();
    filter_size = edge_table_entries.size();

    parallel_for (0, filter_size,[&](size_t i) {
      const uintE v = std::get<0>(edge_table_entries[i]);
      cliques[eltsPerCacheLine*v] -= std::get<1>(edge_table_entries[i]);
      uintE deg = D[v];
      if (deg > cur_bkt) {
        long new_deg = std::max(cliques[eltsPerCacheLine*v], (long) cur_bkt);
        D[v] = new_deg;
        long bkt = b.get_bucket(deg, new_deg);
        D_filter[i] = std::make_tuple(v, bkt);
      } else D_filter[i] = std::make_tuple(UINT_E_MAX, LONG_MAX);
    });
    filter_t.stop();

    auto apply_f = [&](size_t i) -> Maybe<std::tuple<uintE, uintE>> {
      uintE v = std::get<0>(D_filter[i]);
      uintE bkt = std::get<1>(D_filter[i]);
      if (v != UINT_E_MAX && still_active[v] != 2) return wrap(v, bkt);
      return Maybe<std::tuple<uintE, uintE> >();
    };
    bkt_t.start();
    b.update_buckets(apply_f, filter_size);
    bkt_t.stop();
}
    active.del();

    rounds++;
  }
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "max_bkt: " << max_bkt << std::endl;

  bkt_t.reportTotal("bkt time");
  filter_t.reportTotal("filter time");
  updct_t.reportTotal("update count time");

  b.del();
  free(still_active);

  return D;
}