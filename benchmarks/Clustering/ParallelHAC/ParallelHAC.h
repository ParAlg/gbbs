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

#include "UnweightedAverageLinkage.h"

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


//template <class ClusteredGraph, class W>
//void ProcessEdgesCompleteOld(ClusteredGraph& CG, sequence<std::tuple<uintE, uintE, W>>&& edges) {
//  // using edge = std::tuple<uintE, uintE, W>;
//
//  auto V = sequence<vtx_status>::uninitialized(CG.n);
//
//  // Do we need to sort the edges / construct the graph in CSR format?
//  // - need to sample a random neighbor of a vertex (seems weird to implement
//  //   without sorted edges and offsets (CSR))
//  // - can simulate by using writeMin with random hashes (newly generated
//  //   per-round). Do a write
//
//  // 1. Compute the set of vertices that are active (incident to any edges).
//  //    This set shrinks over the course of this algorithm.
//  auto BoolSeq = sequence<bool>(CG.n, false);
//  parallel_for(0, edges.size(), [&] (size_t i) {
//    auto [u, v, wgh] = edges[i];
//    assert(u < CG.n);
//    assert(v < CG.n);
//    if (!BoolSeq[u]) BoolSeq[u] = true;
//    if (!BoolSeq[v]) BoolSeq[v] = true;
//  });
//  auto I = pack_index(BoolSeq);
//
//  // 2. Select red and blue vertices.
//  // TODO: use a better source of randomness (fresh per-round).
//  auto is_red = pbbslib::make_delayed<bool>(I.size(), [&] (size_t i) {
//    uintE u = I[i];
//    return (bool)(pbbslib::hash32_3(u) % 2); });
//  auto Split = parlay::internal::split_two(I, is_red);
//  auto& S = std::get<0>(Split);
//  auto R = S.cut(0, std::get<1>(Split));
//  auto B = S.cut(std::get<1>(Split), S.size());
//  std::cout << "I size = " << I.size() << " R size = " << R.size() << " B size = " << B.size() << std::endl;
//
//  // Blue vertices select a random red neighbor to try and propose to.
//  // Red vertices compute an MIS on their desired neighbors. Conflict edges are
//  // edges with weight < threshold.
//  //
//  // Idea: just pick local minima (if there are conflict edges).
//
//  // Per vertex array stores:
//  // - red/blue identification
//
//
//  // Need to subsequently filter out edges between (red, blue) vertices where a
//  // bad neighbor of the blue vertex merged with the red vertex. How?
//  // Idea (again use local minima): if we have a conflict edge to our chosen red
//  // vertex, and are not the local minima in the conflict nghood, delete the
//  // edge (mark it).
//
//  // After finding all successful merges, filter the edges array. Need to remove
//  // new self-loops, and edges that are deleted due to conflict. Self-loop can
//  // be done using a map.
//  // Keep iterating util edges.size() == 0.
//
//  // Note that the set of active vertices is also shrinking.
//}
//
//template <class ClusteredGraph, class W>
//void ProcessEdgesComplete(ClusteredGraph& CG, sequence<std::tuple<uintE, uintE, W>>&& edges) {
//  using edge = std::tuple<uintE, uintE, W>;
//  auto V = sequence<vtx_status>::uninitialized(CG.n);
//  auto edges_2 = sequence<edge>::uninitialized(edges.size());
//
////  while (edges.size() > 0) {
////
////
////  }
//
//  // 1. Compute the set of vertices that are active (incident to any edges).
//  //    This set shrinks over the course of this algorithm.
//  auto BoolSeq = sequence<bool>(CG.n, false);
//  parallel_for(0, edges.size(), [&] (size_t i) {
//    auto [u, v, wgh] = edges[i];
//    assert(u < CG.n);
//    assert(v < CG.n);
//    if (!BoolSeq[u]) BoolSeq[u] = true;
//    if (!BoolSeq[v]) BoolSeq[v] = true;
//  });
//  auto I = pack_index(BoolSeq);
//
//  // 2. Select red and blue vertices.
//  // TODO: use a better source of randomness (fresh per-round).
//  auto is_red = pbbslib::make_delayed<bool>(I.size(), [&] (size_t i) {
//    uintE u = I[i];
//    return (bool)(pbbslib::hash32_3(u) % 2); });
//  auto Split = parlay::internal::split_two(I, is_red);
//  auto& S = std::get<0>(Split);
//  auto R = S.cut(0, std::get<1>(Split));
//  auto B = S.cut(std::get<1>(Split), S.size());
//  std::cout << "I size = " << I.size() << " R size = " << R.size() << " B size = " << B.size() << std::endl;
//
//  // Blue vertices select a random red neighbor to try and propose to.
//  // Red vertices compute an MIS on their desired neighbors. Conflict edges are
//  // edges with weight < threshold.
//  //
//  // Idea: just pick local minima (if there are conflict edges).
//
//  // Per vertex array stores:
//  // - red/blue identification
//
//
//  // Need to subsequently filter out edges between (red, blue) vertices where a
//  // bad neighbor of the blue vertex merged with the red vertex. How?
//  // Idea (again use local minima): if we have a conflict edge to our chosen red
//  // vertex, and are not the local minima in the conflict nghood, delete the
//  // edge (mark it).
//
//  // After finding all successful merges, filter the edges array. Need to remove
//  // new self-loops, and edges that are deleted due to conflict. Self-loop can
//  // be done using a map.
//  // Keep iterating util edges.size() == 0.
//}


template <class Weights,
// provides get_weight -> Weights::weight_type which is the
// datatype that is stored for each edge incident to a _cluster_. This
// could involve more than simply storing the underlying weight, or
// could internally be a representation like gbbs::empty.
template <class WW> class w_vertex, class IW>  // the weight type of the underlying graph
auto ParallelUPGMA(symmetric_graph<w_vertex, IW>& G, Weights& weights, double epsilon = 0.5) {
  timer tt; tt.start();
  using clustered_graph =
      gbbs::clustering::clustered_graph<Weights, IW, w_vertex>;
  timer bt; bt.start();
  auto CG = clustered_graph(G, weights);
  bt.next("build clustered graph time");

  using Sim = double;
  using W = typename Weights::weight_type;
  using edge = std::tuple<uintE, uintE, W>;

  Sim max_weight = 0;
  Sim min_weight = std::numeric_limits<Sim>::max();
  parallel_for(0, CG.n, [&] (size_t i) {
    auto f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      auto actual_weight = Weights::get_weight(wgh, u, v, CG);
      if (actual_weight > max_weight) { pbbslib::write_max(&max_weight, actual_weight); }
      if (actual_weight < min_weight) { pbbslib::write_min(&min_weight, actual_weight); }
    };
    CG.clusters[i].iterate(i, f);
  });
  assert(max_weight >= min_weight);
  std::cout << "Max weight = " << max_weight << " Min weight = " << min_weight << std::endl;
  auto orig_max_weight = max_weight;

  double one_plus_eps = 1 + epsilon;
  long rounds = max((size_t)ceil(log(max_weight / min_weight) / log(one_plus_eps)), (size_t)1);

  parlay::random rnd;
  auto Colors = sequence<size_t>::uninitialized(G.n);

  while (rounds > 0) {
    Sim lower_threshold = max_weight / one_plus_eps;

    std::cout << "Round = " << rounds << std::endl;
    ProcessGraphUnweightedAverage</*AggressiveMerge=*/true, Weights>(CG, lower_threshold, max_weight, rnd);

    rnd = rnd.next();
    max_weight /= one_plus_eps;
    rounds--;
  }

  auto lower_threshold = max_weight / one_plus_eps;
  if (lower_threshold > 0) {
    // Final round.

    std::cout << "Final round." << std::endl;
    ProcessGraphUnweightedAverage</*AggressiveMerge=*/true, Weights>(CG, (Sim)0, orig_max_weight, rnd);
  }

//  for (size_t i=0; i<G.n; i++) {
//    if (CG.clusters[i].active && CG.clusters[i].neighbor_size() > 0) {
//      std::cout << "i = " << i << " is active. Degree = " <<
//        CG.clusters[i].neighbor_size() << " Cluster size = " <<
//        CG.clusters[i].cluster_size() << std::endl;
//      CG.clusters[i].print_edges();
//    }
//  }
  tt.stop(); tt.reportTotal("total time");
  return CG.get_dendrogram();
}


//    std::cout << "Round = " << rounds << ". Extracting edges with weight between " << lower_threshold << " and " << max_weight << std::endl;
//    timer rt; rt.start();
//    // Map-reduce to figure out the number of edges with weight > threshold.
//    auto wgh_degrees = sequence<size_t>::from_function(CG.n, [&] (size_t i) {
//      // TODO: should only process active clusters
//      auto pred_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
//        assert(wgh.get_weight() <= max_weight);
//        return (wgh.get_weight() >= lower_threshold) && (wgh.get_weight() <= max_weight) && (u < v);   // predicate
//      };
//      return CG.clusters[i].countNeighbors(i, pred_f);
//    });
//    size_t sum_wgh_degrees = scan_inplace(make_slice(wgh_degrees));
//
//    // Extract the edges.
//    auto edges = sequence<edge>::uninitialized(sum_wgh_degrees);
//    parallel_for(0, CG.n, [&] (size_t i) {
//      size_t offset = wgh_degrees[i];
//      size_t delta = ((i == CG.n-1) ? sum_wgh_degrees : wgh_degrees[i+1]) - offset;
//      if (delta > 0) {
//        size_t k = 0;
//        auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
//          if ((wgh.get_weight() >= lower_threshold) && (u < v)) {  // predicate
//            edges[offset + k] = {u, v, wgh};
//            k++;
//          }
//        };
//        CG.clusters[i].iterate(i, map_f);
//      }
//    });
//
//    std::cout << "Extracted: " << edges.size() << " edges." << std::endl;
//    rt.next("Extract time");
//
//    // All edges passing the threshold for this round now in edges. Now process
//    // using random-mate-based algorithm.
//    //
//    // The stuff above is generic for any linkage, but what comes next is specific
//    // to the given linkage function.
//    //ProcessEdgesComplete(CG, std::move(edges));
//    ProcessEdgesUnweightedAverage(CG, std::move(edges), Colors, rnd);
//    rt.next("Process time");


// For log_{1+eps}(...) rounds:
//   - extract active edges in a sequence
//   - process active edges to obtain a sub-clustering at this level
// - emit overall clustering

//auto degrees = sequence<size_t>::from_function(CG.n, [&] (size_t i) { return CG.degree(i); });
//size_t sum_degs = scan_inplace(make_slice(degrees));
//std::cout << "sum_degs = " << sum_degs << std::endl;


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
