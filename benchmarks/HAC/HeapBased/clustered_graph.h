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
#include "pam/pam.h"

namespace gbbs {
namespace clustering {

template <class Weights, class IW, template <class W> class w_vertex>
struct clustered_graph {

  using orig_vertex = w_vertex<IW>;
  using Graph = symmetric_graph<w_vertex, IW>;

  using W = typename Weights::weight_type;
  using edge = std::pair<uintE, W>;

  struct neighbor_entry {
    using key_t = uintE;  // neighbor_id
    using val_t = W;      // weight
    using aug_t = W;      // aggregated weight
    static inline bool comp(key_t a, key_t b) { return a < b; }
    static aug_t get_empty() { return Weights::id(); }  // check
    static aug_t from_entry(key_t k, val_t v) { return v; }  // v
    static aug_t combine(aug_t a, aug_t b) { return Weights::combine(a, b); }  // max
  };

  using neighbor_map = aug_map<neighbor_entry>;

  struct clustered_vertex {

    clustered_vertex() : staleness(0) {}

    clustered_vertex(uintE vtx_id, orig_vertex& vertex, const Weights& weights) {
      auto cluster_size = vertex.out_degree();
      staleness = cluster_size;

      auto edges = sequence<edge>::uninitialized(cluster_size);

      size_t i = 0;
      auto map_f = [&] (const uintE& u, const uintE& v, const IW& wgh) {
        W true_weight = weights.get_weight(u, v, wgh);
        edges[i++] = std::make_pair(v, true_weight);
        if (vtx_id == 0) {
          std::cout << "ngh = " << v << " " << true_weight << std::endl;
        }
      };
      vertex.out_neighbors().map(map_f, /* parallel = */false);

      neighbors = neighbor_map(edges);
      if (vtx_id == 0) {
        std::cout << "aug_val = " << neighbors.aug_val() << std::endl;
      }
    }

    std::optional<edge> highest_priority_edge() {
      if (size() == 0) return {};
      W m = neighbors.aug_val();
      auto f = [m](W v) { return v < m; };
      auto entry = *neighbors.aug_select(f);
      assert(entry.second == m);
      return entry;
    }

    uintE size() {
      return neighbors.size();
    }

    // Tracks the last cluster update size.
    uintE staleness;
    // An augmented map storing our neighbors + weights.
    neighbor_map neighbors;
  };

  Graph& G;
  Weights& weights;
  size_t n;

  // parlay::sequence<uintE> components;  // necessary?
  parlay::sequence<clustered_vertex> clusters;

  // union : (clustered_vertex* vtx_1, clustered_vertex* vtx_2)

  // get_min_wgh(cluster_id id) {
  //
  // }

  clustered_graph(Graph& G, Weights& weights) : G(G), weights(weights) {
    n = G.n;
    clusters = parlay::sequence<clustered_vertex>(n);

    parallel_for(0, n, [&] (size_t i) {
      auto orig = G.get_vertex(i);
      clusters[i] = clustered_vertex(i, orig, weights);
    });
    std::cout << "Built all vertices" << std::endl;
    for (size_t i=0; i<=100; i++) {
      std::cout << clusters[i].size() << " " << G.get_vertex(i).out_degree() << std::endl;
    }
    parallel_for(0, n, [&] (size_t i) {
      assert(clusters[i].size() == G.get_vertex(i).out_degree());
    });
  }

  // extract dendrogram



};


}  // namespace clustering
}  // namespace gbbs
