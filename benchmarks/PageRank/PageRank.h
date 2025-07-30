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
//
// This file provides an implementation of the power-iteration kernel for
// computing PageRank. The implementation works on both undirected and directed
// graphs.
//
// The code handles vertices with zero weighted out-degree by giving them
// (implicit) out-edges to all source nodes (by default, all nodes are
// considered source nodes). The implementation sums up the mass of the zero
// weighted out-degree vertices in each iteration and spreads it evenly across
// all source nodes.
//
// Note that we provide some alternate implementations of the same algorithm in
// other files in this package (PageRank_edgeMapReduce and PageRank_delta). The
// difference for the edgeMapReduce implementration are in how the matrix-vector
// product is computed (in terms of cache-locality). The PageRank_delta
// implementation implements an approximate version of PageRank that is faster
// in many cases but is not yet carefully tested. The difference between
// implementations (1) and (2) are some optimizations used to speed up how the
// matrix-vector product works. We should carefully benchmark the two
// implementations again, but from a few years ago (~2020), the PageRank code
// was consistently faster than PageRank_edgeMap by 20--30% on the WDC2012
// graph.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <type_traits>
#include <utility>

#include "absl/container/flat_hash_set.h"
#include "absl/types/span.h"
#include "gbbs/bridge.h"
#include "gbbs/edge_map_data.h"
#include "gbbs/flags.h"
#include "gbbs/helpers/assert.h"
#include "gbbs/macros.h"
#include "gbbs/vertex_subset.h"
#include "parlay/monoid.h"
#include "parlay/sequence.h"

namespace gbbs {

namespace pagerank_utils {

inline void ValidatePageRankParameters(const double eps,
                                       const double damping_factor) {
  ASSERT(eps >= 0.0);
  ASSERT(0.0 <= damping_factor && damping_factor < 1.0);
}

// Returns a boolean indicating whether `Graph` is weighted.
template <class Graph>
constexpr bool IsWeighted() {
  return !std::is_same_v<typename Graph::weight_type, gbbs::empty>;
}

}  // namespace pagerank_utils

// edgeMap struct used for distributing `_p_curr[s]`, for each vertex `s` in the
// considered subset of nodes (which is assumed to contain only nodes with
// non-zero weighted out-degree) to the `_p_next` values of its out-neighbors.
// The amount distributed to each out-neighbor is proportional to the weight of
// the edge from `s` to the out-neighbor. For unweighted graphs, we treat all
// edges as having weight 1.
template <class Graph>
struct PR_F {
  using W = typename Graph::weight_type;
  const double* p_curr;
  double* p_next;
  const Graph& graph;
  const double* weighted_out_degrees;  // Used only if `Graph` is weighted.

  // `_p_curr` and `_p_next` must have size equal to `G.n`. The same applies to
  // `weighted_out_degrees` if `Graph` is weighted; otherwise,
  // `weighted_out_degrees` is ignored.
  PR_F(const double* _p_curr, double* _p_next, const Graph& graph,
       const double* weighted_out_degrees)
      : p_curr(_p_curr),
        p_next(_p_next),
        graph(graph),
        weighted_out_degrees(weighted_out_degrees) {}

  inline bool update(
      const uintE& s, const uintE& d,
      const W& wgh) {  // update function applies PageRank equation
    p_next[d] += ContributionFromInNeighbor(s, wgh);
    return true;
  }

  inline bool updateAtomic(const uintE& s, const uintE& d,
                           const W& wgh) {  // atomic Update
    gbbs::fetch_and_add(&p_next[d], ContributionFromInNeighbor(s, wgh));
    return true;
  }

  inline bool cond(intT d) { return cond_true(d); }

 private:
  // Returns the amount of mass that node `s` sends over an out-edge with weight
  // `wgh`.
  inline double ContributionFromInNeighbor(const uintE& s, const W& wgh) {
    // Note: the (unweighted or weighted) out-degree of `s` is assumed be
    // non-zero, so the division below is safe.
    if constexpr (pagerank_utils::IsWeighted<Graph>()) {
      return p_curr[s] * (wgh / weighted_out_degrees[s]);
    } else {
      return p_curr[s] / graph.get_vertex(s).out_degree();
    }
  }
};

// Power iteration implementation of PageRank that uses Ligra's edgeMap
// functionality to perform sparse matrix-vector (SpMV) products. The dense
// implementation of edgeMap should be called in every iteration, which will
// require the in-edges of the graph being materialized (this expectation is met
// if using an undirected graph, or if the directed graph has both
// in-/out-edges materialized.
//
// The last argument should only be set to false if the graph you are using does
// not have its in-edges materialized.
//
// The convergence threshold `eps` must be nonnegative. The algorithm stops when
// the L1 distance between the PageRank values (seen as a vector with one entry
// per node) in two consecutive iterations becomes smaller than `eps * (number
// of nodes in G)`, or when `max_iters` iterations have been executed (whichever
// comes first).
// `damping_factor` must be in the range [0, 1).
//
// If the source vector has non-zero length, the algorithm:
// (1) sets the teleportation probabilities (i.e., the personalization vector)
// to 1/|sources| for each source, and 0 for any non-source, and
// (2) sets the initial PageRank scores to 1/|sources| for each source, and 0
// for any non-source.
//
// If the source vector is empty, the algorithm sets the teleportation
// probabilities and the initial PageRank values to 1/n for each node.
template <class Graph>
sequence<double> PageRank_edgeMap(const Graph& G, double eps = 0.000001,
                                  absl::Span<const uintE> sources = {},
                                  double damping_factor = 0.85,
                                  size_t max_iters = 100,
                                  bool has_in_edges = true) {
  pagerank_utils::ValidatePageRankParameters(eps, damping_factor);
  const uintE n = G.n;
  if (n == 0) return sequence<double>();

  double one_over_num_sources = 1 / (double)n;
  // Current PageRank values.
  sequence<double> p_curr;

  absl::flat_hash_set<uintE> sources_set(sources.begin(), sources.end());
  if (!sources.empty()) {
    // If the source set is non-empty, we need to set the initial mass over only
    // the sources.
    one_over_num_sources = 1 / (double)sources.size();
    p_curr = parlay::sequence<double>(n);
    parlay::parallel_for(0, sources.size(), [&](size_t i) {
      p_curr[sources[i]] = one_over_num_sources;
    });
  } else {
    p_curr = sequence<double>(n, one_over_num_sources);
  }

  // Tentative PageRank values for the next iteration.
  auto p_next = sequence<double>(n, 0.0);

  // Compute the weighted out-degrees if the graph is weighted.
  sequence<double> weighted_out_degrees;
  if constexpr (pagerank_utils::IsWeighted<Graph>()) {
    auto get_edge_weight = [](uintE unused_source_id, uintE unused_target_id,
                              const typename Graph::weight_type weight) {
      return double{weight};
    };
    weighted_out_degrees =
        sequence<double>::from_function(n, [&](const size_t i) {
          return G.get_vertex(i).out_neighbors().reduce(get_edge_weight,
                                                        parlay::plus<double>());
        });
  }

  auto is_non_dangling_node = sequence<bool>::from_function(n, [&](uintE v) {
    if constexpr (pagerank_utils::IsWeighted<Graph>()) {
      return weighted_out_degrees.at(v) != 0.0;
    } else {
      return G.get_vertex(v).out_degree() != 0;
    }
  });
  // Nodes with zero (weighted) out-degree.
  parlay::sequence<uintE> dangling_nodes =
      parlay::pack_index<uintE>(parlay::delayed_seq<bool>(
          n, [&](size_t i) { return !is_non_dangling_node.at(i); }));
  // Nodes that will propagate weight to their out-neighbors at each iteration.
  vertexSubset Frontier(n, n, std::move(is_non_dangling_node));

  // If the graph does not have in-edges materialized, use the dense-forward
  // setting.
  gbbs::flags flags = (!has_in_edges) ? gbbs::dense_forward : 0;
  // edgeMap does not require returning an output vertex subset.
  flags |= gbbs::no_output;

  // Each iteration uses `p_curr` to compute `p_next`, then swaps `p_curr` and
  // `p_next`, and sets `p_next` to zero in preparation for the next iteration.
  for (size_t iter = 0; iter < max_iters; ++iter) {
    gbbs_debug(timer t; t.start(););

    // Sum of the PageRank values of the dangling nodes.
    double dangling_sum = parlay::reduce(parlay::delayed_map(
        dangling_nodes, [&](uintE v) { return p_curr[v]; }));

    // SpMV
    edgeMap(G, Frontier,
            PR_F<Graph>(p_curr.begin(), p_next.begin(), G,
                        weighted_out_degrees.data()),
            /*threshold=*/0, flags);
    if (sources.empty()) {
      // Updates p_next by assuming each node has an equal probability being
      // teleported to, and also takes into account the contribution from nodes
      // with no outgoing edges.
      double added_constant = (1 - damping_factor) * one_over_num_sources;
      parlay::parallel_for(0, n, [&](size_t i) {
        p_next[i] += dangling_sum * one_over_num_sources;
        p_next[i] = damping_factor * p_next[i] + added_constant;
      });
    } else {
      // When we have sources, the teleportation probabilities are now uniform
      // across the source set. If desired, we could update this logic in the
      // future to use a user-specified distribution.
      // TODO(laxmand): update this case so that we use a bit-vector for the
      // sources since they are read-only.
      one_over_num_sources = 1 / (double)sources.size();
      double added_constant = (1 - damping_factor) * one_over_num_sources;
      parlay::parallel_for(0, n, [&](size_t i) {
        if (sources_set.contains(i)) {
          p_next[i] += dangling_sum * one_over_num_sources;
          p_next[i] = damping_factor * p_next[i] + added_constant;
        } else {
          p_next[i] = damping_factor * p_next[i];
        }
      });
    }

    // Check convergence: compute L1-norm between p_curr and p_next.
    auto differences = parlay::delayed_seq<double>(
        n, [&](size_t i) { return fabs(p_curr[i] - p_next[i]); });
    double L1_norm = parlay::reduce(differences);

    // Swap p_curr and p_next. The final vector returned will be p_curr.
    std::swap(p_curr, p_next);
    if (L1_norm < eps * n) {
      break;
    }

    gbbs_debug(std::cout << "L1_norm = " << L1_norm << std::endl;);
    // Reset p_curr
    parallel_for(0, n, [&](size_t i) { p_next[i] = 0.0; });

    gbbs_debug(t.stop(); t.next("iteration time"););
  }
  gbbs_debug(auto max_pr = parlay::reduce_max(p_curr);
             std::cout << "max_pr = " << max_pr << std::endl;);
  return p_curr;
}

}  // namespace gbbs
