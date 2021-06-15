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
namespace approx_average_linkage {

template <class Weights, class IW, template <class W> class w_vertex>
struct clustered_graph {

  using orig_vertex = w_vertex<IW>;
  using Graph = symmetric_graph<w_vertex, IW>;

  using W = typename Weights::weight_type;
  using internal_edge = std::pair<uintE, std::pair<uintE, W>>;
  using edge = std::pair<uintE, W>;

  struct neighbor_entry {
    using key_t = uintE;                // neighbor_id
    using val_t = std::pair<uintE, W>;  // (id * weight)
    using aug_t = W;                    // aggregated weight
    static inline bool comp(key_t a, key_t b) { return a < b; }
    static aug_t get_empty() { return Weights::id(); }
    static aug_t from_entry(key_t k, val_t v) { return v.second; }  // (get weight)
    // used to select min/max edges based on similarity/dissimilarity clustering.
    static aug_t combine(aug_t a, aug_t b) { return Weights::augmented_combine(a, b); }
  };

  using neighbor_map = aug_map<neighbor_entry>;

  struct clustered_vertex {

    clustered_vertex() {}

    clustered_vertex(uintE vtx_id, orig_vertex& vertex, const Weights& weights) {
      auto cluster_size = vertex.out_degree();
      num_in_cluster = 1;  // initially just this vertex
      staleness = 1;
      active = true;
      current_id = vtx_id;

      auto edges = sequence<internal_edge>::uninitialized(cluster_size);

      size_t i = 0;
      auto map_f = [&] (const uintE& u, const uintE& v, const IW& wgh) {
        W true_weight = Weights::get_weight(u, v, wgh);
        edges[i++] = std::make_pair(v, std::make_pair(v, true_weight));
      };
      vertex.out_neighbors().map(map_f, /* parallel = */false);

      neighbors = neighbor_map(edges);

      cur_best_edge = best_edge();
    }

    std::optional<edge> best_edge() {
      if (neighbor_size() == 0) return {};
      W m = neighbors.aug_val();
      internal_edge entry;
      entry = *neighbors.aug_eq(m);
      assert(entry.second.second == m);
      return entry.second;
    }

    std::optional<W> best_edge_weight() {
      if (neighbor_size() == 0) return {};
      return neighbors.aug_val();
    }

    void update_best_edge() {
      cur_best_edge = best_edge();
    }

    uintE neighbor_size() {
      return neighbors.size();
    }

    uintE size() {
      return num_in_cluster;
    }

    bool is_active() {
      return active;
    }

    uintE get_current_id() {
      return current_id;
    }

    void set_current_id(uintE id) {
      current_id = id;
    }

    bool is_stale(double epsilon) {
      debug(std::cout << "staleness check: staleness = " << staleness << " eps = " << epsilon << " lhs = " << ((staleness * (1 + epsilon))) << " rhs = " << size() << std::endl;);
      return ((staleness * (1 + epsilon)) < ((double)size()));
    }

    // Tracks the best edge.
    std::optional<edge> cur_best_edge;
    // Tracks the last cluster update size.
    uintE staleness;
    // The "current" id of this cluster, updated upon a merge that keeps this cluster active.
    uintE current_id;
    // Number of vertices contained in this cluster.
    uintE num_in_cluster;
    // Active == false iff this cluster has not yet been clustered
    bool active;
    // An augmented map storing our neighbors + weights.
    neighbor_map neighbors;
  };

  Graph& G;
  Weights& weights;
  double epsilon;
  uintE n;
  uintE last_cluster_id;
  uintE num_merges_performed;

  parlay::sequence<clustered_vertex> clusters;
  parlay::sequence<std::pair<uintE, W>> dendrogram;

  // Returns whether this cluster is still active, or whether it has been merged
  // into a _larger_ cluster.
  bool is_active(uintE id) {
    return clusters[id].is_active();
  }

  uintE new_cluster_id() {
    uintE ret = last_cluster_id;
    last_cluster_id++;
    return ret;
  }

  template <class UH>
  uintE unite(uintE a, uintE b, W wgh, const UH& update_heap) {
    assert(is_active(a));
    assert(is_active(b));

    // Identify smaller/larger clusters (will merge smaller -> larger).
    uintE size_a = clusters[a].neighbor_size();
    uintE size_b = clusters[b].neighbor_size();
    uintE smaller, larger;
    if (size_a < size_b) {
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

    size_t new_cluster_size = clusters[larger].size() + clusters[smaller].size();
    debug(std::cout << "New cluster size = " << new_cluster_size << std::endl;);

    // Map the small map's weights to use the new cluster size. This
    // is important in the case when we merge and don't find something
    // in the intersection. Note for the future that a lazy / delayed map
    // would be OK here.
    auto small_map_f = [&] (const auto& entry) {
      return Weights::RebuildWeight(clusters, entry.first, entry.second, new_cluster_size);
    };
    auto updated_small = neighbor_map::map(std::move(small_pre_merge), small_map_f);
    auto smaller_keys = neighbor_map::keys(updated_small);

    // Merge smaller (updated) and larger, using linkage for elements
    // found in the intersection.
    auto linkage = Weights::GetLinkage(clusters, new_cluster_size);
    auto merged = neighbor_map::map_union(
        std::move(updated_small),
        std::move(large_pre_merge),
        linkage);
    clusters[larger].neighbors = std::move(merged);

    // Save that clusters a and b are merged.
    uintE current_a = clusters[a].get_current_id();
    uintE current_b = clusters[b].get_current_id();
    uintE new_id = new_cluster_id();  // increments next_id
    num_merges_performed++;

    dendrogram[current_a] = {new_id, wgh};
    dendrogram[current_b] = {new_id, wgh};

    // Update the current id of the remaining cluster.
    clusters[larger].current_id = new_id;
    // Update the size of the remaining cluster.
    clusters[larger].num_in_cluster = new_cluster_size;

    // Map over _all_ of smaller's edges, and update its neighbors to point to
    // larger. If the neighbor, w, also has an edge to larger (a
    // smaller-larger-w triangle), then update the weight of this edge.
    for (size_t i=0; i<smaller_keys.size(); i++) {
      uintE w = smaller_keys[i];
      assert(clusters[w].neighbors.contains(smaller));  // Sanity.

      auto w_zero = std::move(clusters[w].neighbors);
      auto found_value = *(w_zero.find(smaller));  // value

      auto w_one = neighbor_map::remove(std::move(w_zero), smaller);

      // Insert larger, merging if it already exists in the tree.
      found_value.first = larger;
      auto new_value = Weights::RebuildWeight(clusters, w, found_value, new_cluster_size);
      auto larger_ent = std::make_pair(larger, new_value);

      auto ngh_linkage = Weights::GetLinkage(clusters, clusters[w].size());
      w_one.insert(larger_ent, ngh_linkage);

      // Move the neighbors back.
      clusters[w].neighbors = std::move(w_one);

      update_heap(w, larger, new_value.second.get_weight());
    }

    // Staleness check.
    if (clusters[larger].is_stale(epsilon)) {
      debug(std::cout << "LARGER = " << larger << " is STALE" << std::endl;);
      // Update our own edges.
      auto edges = std::move(clusters[larger].neighbors);
      auto map_f = [&] (const auto& entry) {
        return Weights::RebuildWeight(clusters, entry.first, entry.second, new_cluster_size);
      };
      auto updated_edges = neighbor_map::map(edges, map_f);
      clusters[larger].neighbors = std::move(updated_edges);

      // Map over the edges, and update on our neighbors endpoints.
      auto update_ngh_f = [&] (const auto& entry) {
        uintE ngh_id = entry.first;
        auto val = entry.second;
        uintE val_id = val.first;
        assert(ngh_id == val_id);
        assert(clusters[ngh_id].active);
        assert(ngh_id != larger);

        val.first = larger; // place our id
        //auto updated_val = Weights::UpdateWeight(clusters, val, new_cluster_size);  // update weight
        auto updated_val = Weights::RebuildWeight(clusters, ngh_id, val, new_cluster_size);  // update weight
        auto new_entry = std::make_pair(larger, updated_val);
        assert(updated_val.first == larger);

        // Now update our neighbor.
        assert(clusters[ngh_id].neighbors.contains(larger));
        clusters[ngh_id].neighbors.insert(new_entry);

        if (larger == 1029) {
          std::cout << "Inserted into ngh_id = " << ngh_id << " wgh = " << updated_val.second.get_weight() << " cut_weight = " << updated_val.second.total_weight <<  std::endl;
        }

        update_heap(ngh_id, larger, updated_val.second.get_weight());
      };
      neighbor_map::map_void(clusters[larger].neighbors, update_ngh_f);

      // Update staleness.
      clusters[larger].staleness = clusters[larger].size();
      debug(std::cout << "Finished update." << std::endl;);
    } else {
      debug(std::cout << "NOT STALE" << std::endl;);
    }

    update_heap(larger, larger, std::numeric_limits<double>::max());

    clusters[larger].update_best_edge();
    return larger;
  }


  clustered_graph(Graph& G, Weights& weights, double epsilon) : G(G), weights(weights), epsilon(epsilon) {
    n = G.n;
    last_cluster_id = n;
    num_merges_performed = 0;
    clusters = parlay::sequence<clustered_vertex>(n);
    dendrogram = parlay::sequence<std::pair<uintE, W>>(2*n - 2, std::make_pair(UINT_E_MAX, W()));

    parallel_for(0, n, [&] (size_t i) {
      auto orig = G.get_vertex(i);
      clusters[i] = clustered_vertex(i, orig, weights);
    });
    std::cout << "Built all vertices" << std::endl;
  }

  // extract dendrogram
  sequence<std::pair<uintE, W>> get_dendrogram() {

    std::cout << "num_merges_performed = " << num_merges_performed << std::endl;
    std::cout << "n = " << n << std::endl;

    if (num_merges_performed < n-1) {
      size_t last_clust = last_cluster_id;
      auto ids = parlay::delayed_seq<uintE>(last_clust + 1, [&] (size_t i) {
        if (dendrogram[i].first == UINT_E_MAX) return (uintE)i;
        return UINT_E_MAX;
      });
      auto bad = parlay::filter(ids, [&] (const uintE& e) { return e != UINT_E_MAX; });

      std::cout << "num bad = " << bad.size() << std::endl;

      std::queue<uintE> bad_queue;
      for (size_t i=0; i<bad.size(); i++) {
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

        std::cout << "Merged components for: " << fst << " " << snd << std::endl;

        bad_queue.push(new_id);
      }
    }

    auto root = (2*n-2);

    for (size_t i=0; i<(2*n - 2); i++) {
      std::cout << "Checking i = " << i << std::endl;
      auto cluster_id = i;
      double wgh = std::numeric_limits<double>::max();
      while (true) {
        auto parent = dendrogram[cluster_id].first;
        auto merge_wgh = dendrogram[cluster_id].second.get_weight();
        std::cout << "id = " << cluster_id << " parent = " << parent << " wgh = " << merge_wgh << std::endl;
        assert(wgh >= merge_wgh);
        wgh = merge_wgh;
        if (cluster_id == parent || parent == root) break;
        cluster_id = parent;
      }
      std::cout << "i = " << i << " is good." << std::endl;
    }


    return std::move(dendrogram);
  }


};


}  // namespace clustering
}  // namespace gbbs
