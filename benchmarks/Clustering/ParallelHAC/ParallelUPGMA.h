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

#include <queue>
#include <unordered_set>
#include <vector>

#include "gbbs/gbbs.h"
#include "ClusteredGraph.h"

namespace gbbs {
namespace clustering {

// Algorithm (high-level):
// 1. Bucket the edges into O(log_{1+epsilon}(W_max / W_min)) edge classes (buckets).
// 2. For the i-th bucket of edges:  // process edges from highest to smallest bucket
//    a. random-mate step, performing (1+epsilon)-close merges.
//    b. subset of edges are no longer in this bucket---move to a lower bucket
//    c. repeat until this bucket is fully contracted.

// For average-linkage:
//   - simply store float edge weights in the graph?
//   - calculate actual weight by checking the sizes of both endpoints.
//   - can leave the weight in the current graph (edge-set) until one of the two
//     endpoints grows larger than (1+eps) of its previous size, at which point we
//     may need to filter it out (if it is no longer in the current bucket).
//   - Big ++ if we can avoid the std::pair<uintE, float/double> weights to
//     represent edges in avg-linkage.
//

// Questions:
//   - How to represent dendrogram output?
//   - Meaningful (large-scale) weighted graph inputs?
//   - What kind of synthetic weights can we claim are interesting? Would be
//     good to run on XXB---XXXB edge graphs to show the real superiority of our
//     approach.

// Need to motivate parallel fork-join setting. Discuss definitions.
//
//
// p_i is the probability of being chosen


// Bucket the edges using a bucketing structure?
// Another option is to do implicit bucketing.
// If we bucket by a bucketing structure, we need to hash the edges or something
// similar so that we can index them in [0, O(m)) as identifiers.
//
// So the edges would be stored both in the hash table, and in the trees.


template <class ClusteredGraph, class W>
void ProcessEdgesComplete(ClusteredGraph& CG, sequence<std::tuple<uintE, uintE, W>>&& edges) {
  using edge = std::tuple<uintE, uintE, W>;


}


template <class Weights,
// provides get_weight : () -> Weights::weight_type which is the
// datatype that is stored for each edge incident to a _cluster_. This
// could involve more than simply storing the underlying weight, or
// could internally be a representation like gbbs::empty.
template <class WW> class w_vertex, class IW>  // the weight type of the underlying graph
auto ParallelUPGMA(symmetric_graph<w_vertex, IW>& G, Weights& weights, double epsilon = 1.0) {
  using clustered_graph =
      gbbs::clustering::clustered_graph<Weights, IW, w_vertex>;
  timer bt; bt.start();
  auto CG = clustered_graph(G, weights);
  bt.next("build clustered graph time");

  using W = typename Weights::weight_type;
  using edge = std::tuple<uintE, uintE, W>;

  // For log_{1+eps}(...) rounds:
  //   - extract active edges in a sequence
  //   - process active edges to obtain a sub-clustering at this level
  // - emit overall clustering

  auto degrees = sequence<size_t>::from_function(CG.n, [&] (size_t i) { return CG.degree(i); });
  size_t sum_degs = scan_inplace(make_slice(degrees));
  std::cout << "sum_degs = " << sum_degs << std::endl;

  // Map-reduce to figure out the number of edges with weight > threshold.
  auto wgh_degrees = sequence<size_t>::from_function(CG.n, [&] (size_t i) {
    // should only process active clusters
    auto pred_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      // check if wgh > threshold for this round. replace first cond with this.
      // symmetry break.
      return true && (u < v);
    };
    return CG.clusters[i].countNeighbors(i, pred_f);
  });
  size_t sum_wgh_degrees = scan_inplace(make_slice(wgh_degrees));
  std::cout << "sum_wgh_degs = " << sum_wgh_degrees << std::endl;

  // Extract the edges.
  auto edges = sequence<edge>::uninitialized(sum_wgh_degrees);
  parallel_for(0, CG.n, [&] (size_t i) {
    size_t offset = wgh_degrees[i];
    size_t delta = ((i == CG.n-1) ? sum_wgh_degrees : wgh_degrees[i+1]) - offset;
    if (delta > 0) {
      size_t k = 0;
      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        edges[offset + k] = {u, v, wgh};
        k++;
      };
      CG.clusters[i].iterate(i, map_f);
    }
  });

  // All edges passing the threshold for this round now in edges. Now process
  // using random-mate-based algorithm.
  //
  // The stuff above is generic for any linkage, but what comes next is specific
  // to the given linkage function.

  ProcessEdgesComplete(CG, std::move(edges));
}



//  size_t n = G.n;
//  size_t round = 0;
//  while ((round < 8) && ((1 << round) < n)) {
//    std::cout << "starting round:" << round << std::endl;
//    size_t stride = 1 << round;
//    size_t end_ind = n/stride;
//    std::cout << "Merging " << end_ind << " many vertices." << std::endl;
//    parallel_for(0, end_ind, [&] (size_t off) {
//      size_t i = off * stride;
//      if (off % 2 == 0) {
//        size_t j = (off+1) * stride;
//        CG.unite(i, j);
//      }
//    });
//    bt.next("Merge time");
//    round++;
//  }


}  // namespace clustering
}  // namespace gbbs
