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
#include "gbbs/bucket.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"
#include "gbbs/pbbslib/dyn_arr.h"
#include "pbbslib/assert.h"
#include "pbbslib/list_allocator.h"
#include "pbbslib/integer_sort.h"

// Ordering files
#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"

// Clique files
#include "intersect.h"
#include "induced_intersection.h"
#include "induced_neighborhood.h"
#include "induced_hybrid.h"
#include "induced_split.h"
#include "relabel.h"
#include "peel.h"

#define SIMD_STATE 4

namespace gbbs {

template <class Graph>
inline uintE* degreeOrderNodes(Graph& G, size_t n) {
  uintE* r = pbbslib::new_array_no_init<uintE>(n); // to hold degree rank per vertex id

  sequence<uintE> o(n); // to hold vertex ids in degree order
  par_for(0, n, pbbslib::kSequentialForThreshold, [&](size_t i){ o[i] = i; });

  pbbs::integer_sort_inplace(o.slice(), [&] (size_t p) {
    return G.get_vertex(p).getOutDegree();
  });

  par_for(0, n, pbbslib::kSequentialForThreshold, 
          [&](size_t i){ r[o[i]] = i; });
  return r;
}

template <class Graph>
pbbs::sequence<uintE> get_ordering(Graph& GA, long order_type, double epsilon = 0.1) {
  const size_t n = GA.n;
  if (order_type == 0) return goodrichpszona_degen::DegeneracyOrder_intsort(GA, epsilon);
  else if (order_type == 1) return barenboimelkin_degen::DegeneracyOrder(GA, epsilon);
  else if (order_type == 2) {
    auto rank = sequence<uintE>(n, [&](size_t i) { return i; });
    auto kcore = KCore(GA);
    auto get_core = [&](uintE p) -> uintE { return kcore[p]; };
    pbbs::integer_sort_inplace(rank.slice(), get_core);
    return rank;
  }
  else if (order_type == 3) return pbbslib::make_sequence(degreeOrderNodes(GA, n), n);
  else if (order_type == 4) {
    auto rank = sequence<uintE>(n, [&](size_t i) { return i; });
    return rank;
  } else ABORT("Unexpected directed type: " << order_type);
}

template <class Graph>
inline size_t TriClique_count(Graph& DG, bool use_base, size_t* per_vert) {
  const size_t n = DG.n;
  size_t count = 0;
  auto counts = sequence<size_t>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) { counts[i] = 0; });

  if (!use_base) { // if counting in total
    auto base_f = [&](uintE a, uintE b, uintE ngh) {};
    count = CountDirectedBalanced(DG, counts.begin(), base_f);
  } else { // if counting per vertex
    auto base_f = [&](uintE a, uintE b, uintE ngh) {
      // Add triangle count to each vertex in triangle
      per_vert[(a+worker_id()*DG.n)]++;
      per_vert[(b+worker_id()*DG.n)]++;
      per_vert[(ngh+worker_id()*DG.n)]++;
    };
    count = CountDirectedBalanced(DG, counts.begin(), base_f);
  }
  return count;
}

template <class Graph>
inline size_t Clique_count(Graph& DG, size_t k, long space_type, bool label, bool use_base, long recursive_level,
  size_t* per_vert) {
  const size_t n = DG.n;
  size_t count = 0;
  // All vertices should be allowed in cliques
  auto use_f = [&](const uintE& src, const uintE& u) { return true; };

  if (!use_base) { // if counting in total
    auto base_f = [&](uintE vtx, size_t _count) {};
    // Everything except space_type = 5 should not be used -- for testing on larger graphs
    if (space_type == 2) count = induced_intersection::CountCliques(DG, k-1, use_f, base_f, use_base);
    else if (space_type == 3) count = induced_neighborhood::CountCliques(DG, k-1);
    else if (space_type == 5) count = induced_hybrid::CountCliques(DG, k-1, base_f, use_base, label, recursive_level);
    else if (space_type == 6) count = induced_split::CountCliques(DG, k-1, base_f, use_base, label, recursive_level);
  } else { // if counting per vertex
    auto base_f = [&](uintE vtx, size_t _count) { per_vert[(vtx+worker_id()*n)] += _count; };
    // Everything except space_type = 5 should not be used -- for testing on larger graphs
    if (space_type == 2) count = induced_intersection::CountCliques(DG, k-1, use_f, base_f, use_base);
    else if (space_type == 3) count = induced_neighborhood::CountCliques(DG, k-1);
    else if (space_type == 5) count = induced_hybrid::CountCliques(DG, k-1, base_f, use_base, label, recursive_level);
    else if (space_type == 6) count = induced_split::CountCliques(DG, k-1, base_f, use_base, label, recursive_level);
  }
  return count;
}


template <class Graph>
inline size_t Clique(Graph& GA, size_t k, long order_type, double epsilon, long space_type, bool label, bool filter,
  bool use_base, long recursive_level, bool approx_peel, double approx_eps) {
  if (k < 3) ABORT("k must be >= 3: " <<  k);

  using W = typename Graph::weight_type;
  const size_t n = GA.n;

  // Store per vertex counts if peeling
  size_t* per_vert = use_base ? (size_t*) calloc(n*num_workers(), sizeof(size_t)) : nullptr;

  // Obtain vertex ordering
  timer t_rank; t_rank.start();
  sequence<uintE> rank = get_ordering(GA, order_type, epsilon);
  double tt_rank = t_rank.stop();
  std::cout << "### Rank Running Time: " << tt_rank << std::endl;

  // Direct the graph based on ordering
  timer t_filter; t_filter.start();
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return (rank[u] < rank[v]) && GA.get_vertex(u).getOutDegree() >= k-1 && GA.get_vertex(v).getOutDegree() >= k-1;
  };
  auto DG = filter ? filterGraph(GA, pack_predicate) : relabel_graph(GA, rank.begin(), pack_predicate);
  double tt_filter = t_filter.stop();
  std::cout << "### Filter Graph Running Time: " << tt_filter << std::endl;


  timer t; t.start();
  size_t count = 0;
  // Clique counting
  if (k == 3) count = TriClique_count(DG, use_base, per_vert);
  else count = Clique_count(DG, k, space_type, label, use_base, recursive_level, per_vert);

  double tt = t.stop();
  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << k << " cliques = " << count << "\n";

  // Return if only counting in total
  if (!use_base) {
    DG.del();
    return count;
  }

  // Collate per vertex counts
  for (size_t j=1; j < static_cast<size_t>(num_workers()); j++) {
    parallel_for(0, n, [&](size_t l) {
      per_vert[l] += per_vert[(l + j*n)];
    });
  }

  // TODO:
  // If the graph was relabeled, move counts back to correct vertices (deprecated)
  // Note that if we want to allow relabeling + peeling, then we should keep an undirected
  // version of the graph that's also relabeled (or have a way to convert from relabeled 
  // vertices to un-relabeled vertices, and use this to query an undirected version
  // of the graph)
  if (!filter) {
    size_t* inverse_per_vert = use_base && !filter ? (size_t*) malloc(n*sizeof(size_t)) : nullptr;
    parallel_for(0, n, [&] (size_t i) { inverse_per_vert[i] = per_vert[rank[i]]; });
    free(per_vert);
    per_vert = inverse_per_vert;
  }

  auto per_vert_seq = pbbslib::make_sequence<size_t>(n, [&] (size_t i) { return per_vert[i]; });
  auto max_per_vert = pbbslib::reduce_max(per_vert_seq);
  if (!approx_peel) {
  // Exact vertex peeling
    if (max_per_vert >= std::numeric_limits<uintE>::max()) Peel<size_t>(GA, DG, k-1, per_vert, label, rank);
    else Peel<uintE>(GA, DG, k-1, per_vert, label, rank);
  } else {
  // Approximate vertex peeling
    if (max_per_vert >= std::numeric_limits<uintE>::max())
      ApproxPeel(GA, DG, k-1, per_vert, count, label, rank, approx_eps);
    else ApproxPeel(GA, DG, k-1, per_vert, count, label, rank, approx_eps);
  }

  // Cleanup
  free(per_vert);
  DG.del();

  return count;
}

}  // namespace gbbs
