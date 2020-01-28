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
    auto base_f = [&](uintE* base) {};
  if (space_type == 2) {
    count = induced_intersection::CountCliques(DG, k-1);
  }
  else if (space_type == 3) {
    count = induced_neighborhood::CountCliques(DG, k-1);
  }
  else if (space_type == 5) {
    count = induced_hybrid::CountCliques(DG, k-1, base_f, use_base, label);
  }
  } else {
    auto base_f = [&](uintE* base) {
      for (size_t i=0; i < k; i++) {
        pbbs::write_add(&(per_vert[eltsPerCacheLine*base[i]]), 1);
      }
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
  auto d_slice = D.slice();
  auto b = make_vertex_buckets(G.n, d_slice, increasing, num_buckets);

  char* still_active = (char*) calloc(G.n, sizeof(char));

  size_t rounds = 0;
  // Peel each bucket
  while (true) {
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    auto active = vertexSubset(G.n, bkt.identifiers);
    size_t cur_bkt = bkt.id;
    if (cur_bkt == b.null_bkt) {
      break;
    }
    active.toSparse();

  for (size_t j=0; j < active.size(); j++) { still_active[active.vtx(j)] = 1; }

// here, update D[i] if necessary
// for each vert in active, just do the same kickoff, but we drop neighbors if they're earlier in the active set
// also drop if already peeled -- check using D
  auto update_d = [&](uintE* base) {
    for (size_t i=0; i <= k; i++) {
      auto tmp = D[base[i]];
      pbbs::write_add(&(D[base[i]]), -1);
      assert (tmp < D[base[i]]);
    }
  };

  sequence<size_t> tots = sequence<size_t>::no_init(active.size());
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active
  auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, G.n, label, true); };
  auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del(); 
  parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active.size(), [&](size_t i, HybridSpace_lw* induced) {
    if (G.get_vertex(active.vtx(i)).getOutDegree() != 0) {
      auto ignore_f = [&](const uintE& u) { return still_active[u] != 2 && (still_active[u] != 1 || u > active.vtx(i)); }; // false if u is dead, false if u is in active and u < active.vtx(i), true otherwise
      induced->setup(G, k, active.vtx(i), ignore_f);
      tots[i] = induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
    } else tots[i] = 0;
  }, 1, false);

    
    auto f = [&](size_t i) -> Maybe<std::tuple<uintE, uintE>> {
      const uintE v = active.vtx(i);
      const uintE v_bkt = D[v];
      uintE bkt = UINT_E_MAX;
      if (!(v_bkt == UINT_E_MAX))
        bkt = b.get_bucket(v_bkt);
      return Maybe<std::tuple<uintE, uintE>>(std::make_tuple(v, bkt));
    };
    b.update_buckets(f, active.size());

    for (size_t j=0; j < active.size(); j++) { still_active[active.vtx(j)] = 2; }

    active.del();
    rounds++;
  }

  b.del();
  free(still_active);

  return D;
}