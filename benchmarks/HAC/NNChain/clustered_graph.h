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
    static aug_t combine(aug_t a, aug_t b) { return Weights::augmented_combine(a, b); }  // used to select min/max edges based on similarity/dissimilarity clustering.
  };

  using neighbor_map = aug_map<neighbor_entry>;

  struct clustered_vertex {

    clustered_vertex() : staleness(0), active(0) {}

    clustered_vertex(uintE vtx_id, orig_vertex& vertex, const Weights& weights) {
      auto cluster_size = vertex.out_degree();
      staleness = cluster_size;
      active = true;

      auto edges = sequence<edge>::uninitialized(cluster_size);

      size_t i = 0;
      auto map_f = [&] (const uintE& u, const uintE& v, const IW& wgh) {
        W true_weight = weights.get_weight(u, v, wgh);
        edges[i++] = std::make_pair(v, true_weight);
      };
      vertex.out_neighbors().map(map_f, /* parallel = */false);

      neighbors = neighbor_map(edges);
    }

    std::optional<edge> highest_priority_edge() {
      if (size() == 0) return {};
      W m = neighbors.aug_val();
      edge entry;
      // if constexpr (Weights::similarity_clustering) {
      //   auto f_sim = [m](W v) { return v < m; };
      //   entry = *neighbors.aug_select(f_sim);
      // } else { // dissimilarity clustering
      //   auto f_dissim = [m](W v) { return v > m; };
      //   entry = *neighbors.aug_select(f_dissim);
      // }
      entry = *neighbors.aug_eq(m);
      if (entry.second != m) {
        // auto entries = neighbor_map::entries(neighbors);
        // std::cout << "Num entries = " << entries.size() << std::endl;
        // for (size_t i=0; i<entries.size(); i++) {
        //   std::cout << entries[i].first << " " << entries[i].second.bundle_size << " " << entries[i].second.total_weight << std::endl;
        // }
      }
      assert(entry.second == m);
      return entry;
    }

    uintE size() {
      return neighbors.size();
    }

    bool is_active() {
      return active;
    }

    // Tracks the last cluster update size.
    uintE staleness;
    // Active == false iff this cluster has not yet been clustered
    bool active;
    // An augmented map storing our neighbors + weights.
    neighbor_map neighbors;
  };

  Graph& G;
  Weights& weights;
  size_t n;

  parlay::sequence<clustered_vertex> clusters;

  // Returns whether this cluster is still active, or whether it has been merged
  // into a _larger_ cluster.
  bool is_active(uintE id) {
    return clusters[id].is_active();
  }

  uintE unite(uintE a, uintE b) {
    assert(is_active(a));
    assert(is_active(b));
    // Identify smaller/larger clusters (will merge smaller -> larger).
    uintE d_a = clusters[a].size();
    uintE d_b = clusters[b].size();
    uintE smaller, larger;
    if (d_a < d_b) {
      smaller = a; larger = b;
    } else {
      larger = a; smaller = b;
    }

    // Deactivate smaller.
    clusters[smaller].active = false;

    // Merge smaller and larger's neighbors.
    auto smaller_ngh = std::move(clusters[smaller].neighbors);
    auto larger_ngh = std::move(clusters[larger].neighbors);

    // Some sanity asserts, we are merging an edge incident to both after all.
    assert(smaller_ngh.size() > 0);
    assert(larger_ngh.size() > 0);

    // Remove larger's id from smaller, and vice versa.
    assert(smaller_ngh.contains(larger));
    assert(larger_ngh.contains(smaller));
    auto small_pre_merge = neighbor_map::remove(std::move(smaller_ngh), larger);
    auto large_pre_merge = neighbor_map::remove(std::move(larger_ngh), smaller);

    auto smaller_keys = neighbor_map::keys(small_pre_merge);

    auto merged = neighbor_map::map_union(
        std::move(small_pre_merge),
        std::move(large_pre_merge),
        Weights::linkage);
    clusters[larger].neighbors = std::move(merged);

    // Map over _all_ of smaller's edges, and update its neighbors to point to
    // larger. If the neighbor, w, also has an edge to larger (a
    // smaller-larger-w triangle), then update the weight of this edge.

    // TODO: check if parallel_for really helps here.
    //for (size_t i=0; i<smaller_keys.size(); i++) {
    parallel_for(0, smaller_keys.size(), [&] (size_t i) {
      uintE w = smaller_keys[i];
      assert(clusters[w].neighbors.contains(smaller));  // Sanity.

      auto w_zero = std::move(clusters[w].neighbors);
      auto found_weight = *(w_zero.find(smaller));
      auto w_one = neighbor_map::remove(std::move(w_zero), smaller);

      // Insert larger, merging using Weights::linkage if it already exists in
      // the tree.
      auto larger_ent = std::make_pair(larger, found_weight);
      w_one.insert(larger_ent, Weights::linkage);

      // Move the neighbors back.
      clusters[w].neighbors = std::move(w_one);
    });
    return larger;
  }


  clustered_graph(Graph& G, Weights& weights) : G(G), weights(weights) {
    n = G.n;
    clusters = parlay::sequence<clustered_vertex>(n);

    parallel_for(0, n, [&] (size_t i) {
      auto orig = G.get_vertex(i);
      clusters[i] = clustered_vertex(i, orig, weights);
    });
    std::cout << "Built all vertices" << std::endl;
    // parallel_for(0, n, [&] (size_t i) {
    //   assert(clusters[i].size() == G.get_vertex(i).out_degree());
    // });
  }

  // extract dendrogram



};


}  // namespace clustering
}  // namespace gbbs
