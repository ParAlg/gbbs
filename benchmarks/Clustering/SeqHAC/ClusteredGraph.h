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
    static aug_t get_empty() { return Weights::id(); }       // check
    static aug_t from_entry(key_t k, val_t v) { return v; }  // v
    static aug_t combine(aug_t a, aug_t b) {
      return Weights::augmented_combine(a, b);
    }  // used to select min/max edges based on similarity/dissimilarity
       // clustering.
  };

  using neighbor_map = aug_map<neighbor_entry>;

  struct clustered_vertex {
    clustered_vertex() : staleness(0), active(0) {}

    clustered_vertex(uintE vtx_id, orig_vertex& vertex,
                     const Weights& weights) {
      auto cluster_size = vertex.out_degree();
      staleness = cluster_size;
      active = true;
      current_id = vtx_id;

      auto edges = sequence<edge>::uninitialized(cluster_size);

      size_t i = 0;
      auto map_f = [&](const uintE& u, const uintE& v, const IW& wgh) {
        W true_weight = weights.get_weight(u, v, wgh);
        edges[i++] = std::make_pair(v, true_weight);
      };
      vertex.out_neighbors().map(map_f, /* parallel = */ false);

      neighbors = neighbor_map(edges);
    }

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
    std::optional<edge> highest_priority_edge() {
      if (size() == 0) return std::optional<edge>();
      W m = neighbors.aug_val();
      edge entry;
      entry = *neighbors.aug_eq(m);
      assert(entry.second == m);
      return entry;
    }
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

    uintE size() { return neighbors.size(); }

    bool is_active() { return active; }

    uintE get_current_id() { return current_id; }

    void set_current_id(uintE id) { current_id = id; }

    // Tracks the last cluster update size.
    uintE staleness;
    // The "current" id of this cluster, updated upon a merge that keeps this
    // cluster active.
    uintE current_id;
    // Active == false iff this cluster has not yet been clustered
    bool active;
    // An augmented map storing our neighbors + weights.
    neighbor_map neighbors;
  };

  Graph& G;
  Weights& weights;
  uintE n;
  uintE last_cluster_id;
  uintE num_merges_performed;

  parlay::sequence<clustered_vertex> clusters;
  parlay::sequence<std::pair<uintE, W>> dendrogram;

  // Returns whether this cluster is still active, or whether it has been merged
  // into a _larger_ cluster.
  bool is_active(uintE id) { return clusters[id].is_active(); }

  uintE new_cluster_id() {
    uintE ret = last_cluster_id;
    last_cluster_id++;
    return ret;
  }

  uintE unite(uintE a, uintE b, W wgh) {
    assert(is_active(a));
    assert(is_active(b));
    // Identify smaller/larger clusters (will merge smaller -> larger).
    uintE d_a = clusters[a].size();
    uintE d_b = clusters[b].size();
    uintE smaller, larger;
    if (d_a < d_b) {
      smaller = a;
      larger = b;
    } else {
      larger = a;
      smaller = b;
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

    auto merged =
        neighbor_map::map_union(std::move(small_pre_merge),
                                std::move(large_pre_merge), Weights::linkage);
    clusters[larger].neighbors = std::move(merged);

    // Save that clusters a and b are merged.
    uintE current_a = clusters[a].get_current_id();
    uintE current_b = clusters[b].get_current_id();
    uintE new_id = new_cluster_id();  // increments next_id
    num_merges_performed++;

    dendrogram[current_a] = {new_id, wgh};
    dendrogram[current_b] = {new_id, wgh};

    // Update the current id of the remaining vertex.
    clusters[larger].current_id = new_id;

    // Map over _all_ of smaller's edges, and update its neighbors to point to
    // larger. If the neighbor, w, also has an edge to larger (a
    // smaller-larger-w triangle), then update the weight of this edge.

    for (size_t i = 0; i < smaller_keys.size(); i++) {
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
    }
    return larger;
  }

  clustered_graph(Graph& G, Weights& weights) : G(G), weights(weights) {
    n = G.n;
    last_cluster_id = n;
    num_merges_performed = 0;
    clusters = parlay::sequence<clustered_vertex>(n);
    dendrogram = parlay::sequence<std::pair<uintE, W>>(
        2 * n - 1, std::make_pair(UINT_E_MAX, W()));

    parallel_for(0, n, [&](size_t i) {
      auto orig = G.get_vertex(i);
      clusters[i] = clustered_vertex(i, orig, weights);
    });
    debug(std::cout << "Built all vertices" << std::endl;);
  }

  // extract dendrogram
  sequence<std::pair<uintE, W>> get_dendrogram() {
    debug(std::cout << "num_merges_performed = " << num_merges_performed
                    << std::endl;);
    debug(std::cout << "n = " << n << std::endl;);

    if (num_merges_performed < n - 1) {
      size_t last_clust = last_cluster_id;
      auto ids = parlay::delayed_seq<uintE>(last_clust + 1, [&](size_t i) {
        if (dendrogram[i].first == UINT_E_MAX) return (uintE)i;
        return UINT_E_MAX;
      });
      auto bad =
          parlay::filter(ids, [&](const uintE& e) { return e != UINT_E_MAX; });

      debug(std::cout << "num bad = " << bad.size() << std::endl;);

      std::queue<uintE> bad_queue;
      for (size_t i = 0; i < bad.size(); i++) {
        bad_queue.push(bad[i]);
      }

      while (bad_queue.size() > 1) {
        uintE fst = bad_queue.front();
        bad_queue.pop();
        uintE snd = bad_queue.front();
        bad_queue.pop();

        uintE new_id = new_cluster_id();  // increments next_id
        dendrogram[fst] = {new_id, Weights::id()};
        dendrogram[snd] = {new_id, Weights::id()};

        debug(std::cout << "Merged components for: " << fst << " " << snd
                        << " dend_size = " << dendrogram.size() << std::endl;);

        bad_queue.push(new_id);
      }
    }

    return std::move(dendrogram);
  }
};

}  // namespace clustering
}  // namespace gbbs
