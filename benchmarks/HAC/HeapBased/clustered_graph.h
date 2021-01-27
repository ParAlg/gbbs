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

#include "benchmarks/Connectivity/UnionFind/union_find_rules.h"

namespace gbbs {
namespace clustering {

template <class Weights, class IW, template <class W> class w_vertex>
struct clustered_graph {

  using orig_vertex = w_vertex<IW>;
  using Graph = symmetric_graph<w_vertex, IW>;

  using W = typename Weights::weight_type;
  using edge = std::pair<uintE, W>;

  struct clustered_vertex {
    uintE staleness;  // tracks the last cluster update size

    // An augmented map storing our neighbors. The augmented map stores our
    // neighbors, with an augmented value of the *weight*, using
    aug_map<entry> neighbors;

    // if |ordered_edges| > 2*|neighbors|, clear and reinsert neighbors
    // (amortized constant work)
    // Todo: What did this comment mean?

    // constexpr bool get_edge_priority(const edge& lhs, const edge& rhs) {
    //   return Weights::less(lhs.second, rhs.second);
    // };
    // std::priority_queue<edge, std::vector<edge>, get_edge_priority> ordered_edges;

    clustered_vertex(orig_vertex& vertex, const Weights& weights) {
      cluster_size = vtx.out_degree();
      staleness = cluster_size;

      auto map_f = [&] (const uintE& u, const uintE& v, const IW& wgh) {
        W true_weight = weights.get_weight(u, v, wgh);
        neighbors.insert({v, true_weight});

        ordered_edges.insert({v, true_weight});
      };
      vertex.out_neighbors().map(map_f, /* parallel = */false);
    }

    uintE size() {
      return neighbors.size();
    }

  };

  Graph& G;
  size_t n;

  parlay::sequence<uintE> components;
  parlay::sequence<clustered_vertex> clusters;

  // union : (clustered_vertex* vtx_1, clustered_vertex* vtx_2)

  // get_min_wgh(cluster_id id) {
  //
  // }


  // extract dendrogram



};


}  // namespace clustering
}  // namespace gbbs
