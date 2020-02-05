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
  pbbslib::sample_sort_inplace(o.slice(), [&](const uintE u, const uintE v) {
    return G.get_vertex(u).getOutDegree() < G.get_vertex(v).getOutDegree();
  });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { r[o[i]] = i; });
  t.stop();
  debug(t.reportTotal("Rank time"););
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
inline size_t KClique(Graph& GA, size_t k, long order_type, double epsilon,
long space_type, bool label, bool filter, bool use_base, uintE* per_vert) {
  std::cout << "### Starting clique counting" << std::endl;
  const size_t eltsPerCacheLine = 64/sizeof(long);
  using W = typename Graph::weight_type;
  assert (k >= 1);
  if (k == 1) return GA.n;
  else if (k == 2) return GA.m;

  timer t_rank; t_rank.start();
  sequence<uintE> rank = get_ordering(GA, order_type, epsilon);
  double tt_rank = t_rank.stop();
  std::cout << "### Rank Running Time: " << tt_rank << std::endl;

  timer t_filter; t_filter.start();
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return (rank[u] < rank[v]) && GA.get_vertex(u).getOutDegree() >= k-1 && GA.get_vertex(v).getOutDegree() >= k-1;
  };
  auto DG = filter ? filter_graph(GA, pack_predicate) : relabel_graph(GA, rank.begin(), pack_predicate); //TODO see if relabel is really needed or not
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
    count = induced_split::CountCliques(DG, k-1, base_f, use_base, label, 300);
  }
  } else {
    auto base_f = [&](uintE vtx, size_t count) {
      pbbs::write_add(&(per_vert[eltsPerCacheLine*vtx]), count);
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
    count = induced_split::CountCliques(DG, k-1, base_f, use_base, label, 300);
  }
  }

  double tt = t.stop();
  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << k << " cliques = " << count << "\n";
  return count;
}








template <class Graph>
sequence<uintE> Peel(Graph& G, size_t k, uintE* cliques, bool label=true, size_t num_buckets=128) {
  const size_t eltsPerCacheLine = 64/sizeof(long);
  auto D = sequence<uintE>(G.n, [&](size_t i) { return cliques[eltsPerCacheLine*i]; });
  //auto ER = sequence<uintE>(G.n, [&](size_t i) { return 0; });
  auto D_update = sequence<uintE>(G.n, [&](size_t i) { return 0; });
  auto D_filter = sequence<std::tuple<uintE, uintE>>(G.n);
  auto b = make_vertex_buckets(G.n, D, increasing, num_buckets);

  char* still_active = (char*) calloc(G.n, sizeof(char));
  for (size_t i=0; i < G.n; i++) {still_active[i] = 0;}

  size_t rounds = 0;
  size_t finished = 0;
  uintE cur_bkt = 0;
  uintE max_bkt = 0;
  // Peel each bucket
   while (finished != G.n) {
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    auto active = vertexSubset(G.n, bkt.identifiers);
    finished += active.size();
    cur_bkt = bkt.id;
    max_bkt = std::max(cur_bkt, (uintE) bkt.id);
    //active.toSparse();

  for (size_t j=0; j < active.size(); j++) { still_active[active.vtx(j)] = 1; }

// here, update D[i] if necessary
// for each vert in active, just do the same kickoff, but we drop neighbors if they're earlier in the active set
// also drop if already peeled -- check using D
  

  //sequence<size_t> tots = sequence<size_t>::no_init(active.size());
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active
  auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, G.n, label, true); };
  auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del(); 
  
  
  parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active.size(), [&](size_t i, HybridSpace_lw* induced) {
    if (G.get_vertex(active.vtx(i)).getOutDegree() != 0) {
      auto ignore_f = [&](const uintE& u) { return still_active[u] != 2 && (still_active[u] != 1 || u > active.vtx(i)); }; // false if u is dead, false if u is in active and u < active.vtx(i), true otherwise
      induced->setup(G, k, active.vtx(i), ignore_f);
      auto update_d = [&](uintE vtx, size_t count) { assert(ignore_f(vtx)); pbbs::write_add(&(D_update[vtx]), count); };
      induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
      //update_d(active.vtx(i), tots[i]);
    } //else tots[i] = 0;
  }, 1, false);

  for (size_t j=0; j < active.size(); j++) { still_active[active.vtx(j)] = 2; }

  // filter D_update for nonzero elements
  // subtract these from D and then we can rebucket these elements
  //auto D_delayed_f = [&](size_t i) { return std::make_tuple(i, D_update[i]); };
  //auto D_delayed = pbbs::delayed_sequence<std::tuple<uintE, uintE>, decltype(D_delayed_f)>(G.n, D_delayed_f);
  //auto D_filter_f = [&](std::tuple<uintE,uintE> tup) { return std::get<1>(tup) > 0; } ;
  //size_t filter_size = pbbs::filter_out(D_delayed, D_filter, D_filter_f);

  size_t filter_size = 0;
  for (size_t l=0; l < G.n; l++) {
    if (D_update[l] > 0) {
      D_filter[filter_size] = std::make_tuple(l, D_update[l]);
      assert (cliques[eltsPerCacheLine*l] >= D_update[l]);
      cliques[eltsPerCacheLine*l] -= D_update[l];
      D_update[l] = 0;
      filter_size++;
    }
  }

  parallel_for(0, filter_size, [&] (size_t i) {
    const uintE v = std::get<0>(D_filter[i]);
    assert (v < G.n);
    //D_update[v] = 0;
    //assert (cliques[v] >= std::get<1>(D_filter[i]));
    //cliques[v] -= std::get<1>(D_filter[i]);
    uintE deg = D[v];
    if (deg > cur_bkt) {
      uintE new_deg = std::max(cliques[eltsPerCacheLine*v], cur_bkt);
      D[v] = new_deg;
      uintE bkt = b.get_bucket(deg, new_deg);
      // store (v, bkt) in an array now, pass it to apply_f below instead of what's there right now -- maybe just store it in D_filter?
      D_filter[i] = std::make_tuple(v, bkt);
    } else D_filter[i] = std::make_tuple(UINT_E_MAX, UINT_E_MAX);
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

  b.del();
  free(still_active);

  return D;
}