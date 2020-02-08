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

  if (!use_base) {
    auto base_f = [&](uintE vtx, size_t count) {};
  if (space_type == 2) {
    count = induced_intersection::CountCliques(DG, k-1);
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
    count = induced_intersection::CountCliques(DG, k-1);
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



template <class Graph>
sequence<long> Peel(Graph& G, size_t k, long* cliques, bool label, sequence<uintE> &rank, size_t num_buckets=16) {
  const size_t eltsPerCacheLine = 64/sizeof(long);
  auto D = sequence<long>(G.n, [&](size_t i) { return cliques[eltsPerCacheLine*i]; });
  //auto ER = sequence<uintE>(G.n, [&](size_t i) { return 0; });
  auto D_update = sequence<long>(eltsPerCacheLine*G.n);
  parallel_for(0, G.n, [&](size_t j){D_update[eltsPerCacheLine*j] = 0;});
  auto D_filter = sequence<std::tuple<uintE, long>>(G.n);
  auto b = make_vertex_buckets(G.n, D, increasing, num_buckets);

  char* still_active = (char*) calloc(G.n, sizeof(char));

  size_t rounds = 0;
  size_t finished = 0;
  long cur_bkt = 0;
  long max_bkt = 0;
  // Peel each bucket
  while (finished != G.n) {
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    auto active = vertexSubset(G.n, bkt.identifiers);
    finished += active.size();
    cur_bkt = bkt.id;
    max_bkt = std::max(cur_bkt, (long) bkt.id);
    //std::cout << "Fetching bucket: " << cur_bkt << std::endl;
    //active.toSparse();

  for (size_t j=0; j < active.size(); j++) { still_active[active.vtx(j)] = 1; }

// here, update D[i] if necessary
// for each vert in active, just do the same kickoff, but we drop neighbors if they're earlier in the active set
// also drop if already peeled -- check using D

  //sequence<size_t> tots = sequence<size_t>::no_init(active.size());
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active
  auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, G.n, label, true); };
  auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del();

if (active.size() < 30) {
  HybridSpace_lw* induced = new HybridSpace_lw();
  induced->alloc(max_deg, k, G.n, label, true); 
    for (size_t i=0; i < active.size(); i++) {
    if (G.get_vertex(active.vtx(i)).getOutDegree() != 0) {
      auto ignore_f = [&](const uintE& u, const uintE& v) {
        if (still_active[u] == 2 || still_active[v] == 2) return false;
        if (still_active[u] == 1 && still_active[v] == 0) return true;
        if (still_active[u] == 0 && still_active[v] == 1) return false;
        return rank[u] < rank[v];
        //return still_active[u] != 2 && (still_active[u] != 1 || u > active.vtx(i));
      }; // false if u is dead, false if u is in active and u < active.vtx(i), true otherwise
      induced->setup(G, k, active.vtx(i), ignore_f);
      auto update_d = [&](uintE vtx, size_t count) { pbbslib::xadd(&(D_update[eltsPerCacheLine*vtx]), (long) count); };
      induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
      //update_d(active.vtx(i), tots[i]);
    } //else tots[i] = 0;
  }
  if (induced != nullptr) { delete induced; }
}
else {
  parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active.size(), [&](size_t i, HybridSpace_lw* induced) {
    if (G.get_vertex(active.vtx(i)).getOutDegree() != 0) {
      auto ignore_f = [&](const uintE& u, const uintE& v) {
        if (still_active[u] == 2 || still_active[v] == 2) return false;
        if (still_active[u] == 1 && still_active[v] == 0) return true;
        if (still_active[u] == 0 && still_active[v] == 1) return false;
        return rank[u] < rank[v];
        //return still_active[u] != 2 && (still_active[u] != 1 || u > active.vtx(i));
      }; // false if u is dead, false if u is in active and u < active.vtx(i), true otherwise
      induced->setup(G, k, active.vtx(i), ignore_f);
      auto update_d = [&](uintE vtx, size_t count) { pbbslib::xadd(&(D_update[eltsPerCacheLine*vtx]), (long) count); };
      induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
      //update_d(active.vtx(i), tots[i]);
    } //else tots[i] = 0;
  }, 1, false);
}
  for (size_t j=0; j < active.size(); j++) { still_active[active.vtx(j)] = 2; }

  // filter D_update for nonzero elements
  // subtract these from D and then we can rebucket these elements
  auto D_delayed_f = [&](size_t i) { return std::make_tuple(i, D_update[eltsPerCacheLine*i]); };
  auto D_delayed = pbbslib::make_sequence<std::tuple<uintE, long>>(G.n, D_delayed_f);
  auto D_filter_f = [&](const std::tuple<uintE, long>& tup) { return std::get<1>(tup) > 0; } ;
  size_t filter_size = pbbs::filter_out(D_delayed, D_filter.slice(), D_filter_f);

  /*size_t filter_size = 0;
  for (size_t l=0; l < G.n; l++) {
    if (D_update[l] > 0) {
      D_filter[filter_size] = std::make_tuple(l, D_update[l]);
      assert (cliques[eltsPerCacheLine*l] >= D_update[l]);
      cliques[eltsPerCacheLine*l] -= D_update[l];
      D_update[l] = 0;
      filter_size++;
    }
  }*/

  parallel_for(0, filter_size, [&] (size_t i) {
    const uintE v = std::get<0>(D_filter[i]);
    assert (v < G.n);
    D_update[eltsPerCacheLine*v] = 0;
    assert (cliques[eltsPerCacheLine*v] >= std::get<1>(D_filter[i]));
    cliques[eltsPerCacheLine*v] -= std::get<1>(D_filter[i]);
    uintE deg = D[v];
    if (deg > cur_bkt) {
      long new_deg = std::max(cliques[eltsPerCacheLine*v], (long) cur_bkt);
      D[v] = new_deg;
      long bkt = b.get_bucket(deg, new_deg);
      // store (v, bkt) in an array now, pass it to apply_f below instead of what's there right now -- maybe just store it in D_filter?
      D_filter[i] = std::make_tuple(v, bkt);
    } else D_filter[i] = std::make_tuple(UINT_E_MAX, LONG_MAX);
  });

  auto apply_f = [&](size_t i) -> Maybe<std::tuple<uintE, uintE>> {
    const uintE v = std::get<0>(D_filter[i]);
    const uintE bkt = std::get<1>(D_filter[i]);
    if (v != UINT_E_MAX) return wrap(v, bkt);
    return Maybe<std::tuple<uintE, uintE> >();
    //ret.exists = std::get<0>(D[v]);
    //return ret;
  };
  b.update_buckets(apply_f, filter_size);

    // so update buckets can take fn rep
    // so what we should do is use apply_f for 0 to filter_size, and then create

    // here, we use apply_f for 0 to filter_size, get a tuple<uintE,long>* to update buckets
    // use another filter for that
    // then put it through update_buckets

    //auto moved = edgeMapData<uintE>(
    //    G, active, kcore_fetch_add<W>(ER.begin(), D.begin(), cur_bkt));
    //vertexMap(moved, apply_f);

    //if (moved.dense()) {
    //  b.update_buckets(moved.get_fn_repr(), n);
    //} else {
    //  b.update_buckets(moved.get_fn_repr(), moved.size());
    //}
    //moved.del();
    active.del();

    rounds++;
  }
   std::cout << "rho = " << rounds << std::endl;

  b.del();
  free(still_active);

  return D;
}