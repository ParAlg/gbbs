#ifndef THIRD_PARTY_GBBS_BENCHMARKS_PAGERANK_PAGERANK_EDGEMAPREDUCE_H_
#define THIRD_PARTY_GBBS_BENCHMARKS_PAGERANK_PAGERANK_EDGEMAPREDUCE_H_

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
// This file provides an alternate implementation of computing PageRank. Please
// see PageRank.h for comments that also apply to this file. The difference 
// between this implementation and the one in PageRank.h is some optimizations
// used to speed up how the matrix-vector product works. We should carefully
// benchmark the two implementations again, but from a few years ago (~2020),
// the PageRank code was consistently faster than PageRank_edgeMap by 20--30% on
// the WDC2012 graph.
//
// TODOs(laxmand):
// - There are unit tests for the first two implementations, but unit tests need
//   to be added for PageRankDelta.
// - PageRankDelta needs to be updated to handle dangling edges.
// - PageRankDelta needs to be updated to handle sources.
// - Add and test support for weighted graphs.
// - Benchmark PageRank_edgeMap and PageRank_edgeMapReduce and update the
//   performance numbers above.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "gbbs/bridge.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/flags.h"
#include "gbbs/helpers/assert.h"
#include "gbbs/macros.h"
#include "gbbs/vertex_subset.h"
#include "parlay/monoid.h"
#include "parlay/sequence.h"

namespace gbbs {

namespace pagerank_edgemapreduce_utils {

inline void ValidatePageRankParameters(const double eps,
                                       const double damping_factor) {
  ASSERT(eps >= 0.0);
  ASSERT(0.0 <= damping_factor && damping_factor < 1.0);
}  
  
}

// This version of PageRank uses edgeMapReduce_dense, an implementation
// which reduces over the in-neighbors of every vertex and aggregates the
// incoming contributions to each vertex in parallel.
//
// The key difference between PageRank_edgeMapReduce and PageRank_edgeMap
// (above) is that PageRank_edgeMap will *sequentially* aggregate the incoming
// contributions to a vertex, whereas PageRank_edgeMapReduce will do this
// reduction in parallel.
template <class Graph>
sequence<double> PageRank_edgeMapReduce(Graph& G, double eps = 0.000001,
                                        std::vector<uintE> sources = {},
                                        double damping_factor = 0.85,
                                        size_t max_iters = 100) {
  pagerank_edgemapreduce_utils::ValidatePageRankParameters(eps, damping_factor);
  using W = typename Graph::weight_type;
  const uintE n = G.n;
  if (n == 0) return sequence<double>();
  double added_constant;
  if (sources.empty()) {
    added_constant = (1 - damping_factor) * (1 / static_cast<double>(n));
  } else {
    added_constant =
        (1 - damping_factor) * (1 / static_cast<double>(sources.size()));
  }
  absl::flat_hash_set<uintE> sources_set(sources.begin(), sources.end());

  double one_over_n = 1 / (double)n;
  sequence<double> p_curr;
  sequence<double> p_div;
  double one_over_num_sources = 0.0;
  if (!sources.empty()) {
    one_over_num_sources = 1 / (double)sources.size();
    p_curr = parlay::sequence<double>(n);
    p_div = sequence<double>(n);
    parlay::parallel_for(0, sources.size(), [&](size_t i) {
      uintE source = sources[i];
      p_curr[source] = one_over_num_sources;
      p_div[source] =
          one_over_n /
          std::max(double{1},
                   static_cast<double>(G.get_vertex(source).out_degree()));
    });
  } else {
    p_curr = sequence<double>(n, one_over_n);
    p_div = sequence<double>::from_function(n, [&](size_t i) -> double {
      return one_over_n /
             std::max(double{1},
                      static_cast<double>(G.get_vertex(i).out_degree()));
    });
  }

  auto p_next = sequence<double>(n, static_cast<double>(0));
  auto frontier = sequence<bool>(n, true);
  auto p_div_next = sequence<double>(n);

  // read from special array of just degrees
  auto degrees = sequence<uintE>::from_function(
      n, [&](size_t i) { return G.get_vertex(i).out_degree(); });

  parlay::sequence<uintE> dangling_nodes = parlay::pack_index<uintE>(
      parlay::delayed_seq<bool>(n, [&](size_t i) { return degrees[i] == 0; }));

  double dangling_sum{0};

  vertexSubset Frontier(n, n, std::move(frontier));
  auto EM = EdgeMap<double, Graph>(
      G, std::make_tuple(UINT_E_MAX, static_cast<double>(0)),
      (size_t)G.m / 1000);

  auto cond_f = [&](const uintE& v) { return true; };
  auto map_f = [&](const uintE& d, const uintE& s, const W& wgh) -> double {
    return p_div[s];
  };
  auto reduce_f = [&](double l, double r) { return l + r; };
  auto apply_f_no_source = [&](std::tuple<uintE, double> k)
      -> std::optional<std::tuple<uintE, double>> {
    uintE u = std::get<0>(k);
    double contribution = std::get<1>(k);
    contribution += dangling_sum * one_over_n;
    p_next[u] = damping_factor * contribution + added_constant;
    p_div_next[u] = (p_next[u] / static_cast<double>(degrees[u]));
    return std::nullopt;
  };
  auto apply_f_source = [&](std::tuple<uintE, double> k)
      -> std::optional<std::tuple<uintE, double>> {
    uintE u = std::get<0>(k);
    double contribution = std::get<1>(k);
    if (sources_set.contains(u)) {
      contribution += dangling_sum * one_over_num_sources;
      p_next[u] = damping_factor * contribution + added_constant;
    } else {
      p_next[u] = damping_factor * contribution;
    }
    p_div_next[u] = p_next[u] / static_cast<double>(degrees[u]);
    return std::nullopt;
  };

  for (size_t iter = 0; iter < max_iters; ++iter) {
    dangling_sum = parlay::reduce(parlay::delayed_map(
        dangling_nodes, [&](uintE v) { return p_curr[v]; }));

    timer t;
    t.start();
    // SpMV
    timer tt;
    tt.start();
    // Ensure we map over the in-edges here.
    if (sources.empty()) {
      EM.template edgeMapReduce_dense<double, double>(
          Frontier, cond_f, map_f, reduce_f, apply_f_no_source, 0.0,
          no_output | in_edges);
    } else {
      EM.template edgeMapReduce_dense<double, double>(
          Frontier, cond_f, map_f, reduce_f, apply_f_source, 0.0,
          no_output | in_edges);
    }
    tt.stop();
    tt.next("em time");

    // Check convergence: compute L1-norm between p_curr and p_next
    auto differences = parlay::delayed_seq<double>(n, [&](size_t i) {
      auto d = p_curr[i];
      p_curr[i] = 0;
      return fabs(d - p_next[i]);
    });
    double L1_norm = parlay::reduce(differences, parlay::plus<double>());
    std::swap(p_curr, p_next);
    // Reset p_curr and p_div.
    std::swap(p_div, p_div_next);
    if (L1_norm < eps) {
      break;
    }
    gbbs_debug(std::cout << "L1_norm = " << L1_norm << std::endl;);

    t.stop();
    t.next("iteration time");
  }
  gbbs_debug(auto max_pr = parlay::reduce_max(p_curr);
             std::cout << "max_pr = " << max_pr << std::endl;);
  return p_curr;
}


  
}

#endif  // THIRD_PARTY_GBBS_BENCHMARKS_PAGERANK_PAGERANK_EDGEMAPREDUCE_H_
