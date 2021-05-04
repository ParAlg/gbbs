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

  using Graph = symmetric_graph<w_vertex, IW>;
  using Graph_vertex = w_vertex<IW>;

  using W = typename Weights::weight_type;
  using edge = std::pair<uintE, W>;

  struct neighbor_entry {
    using key_t = uintE;  // neighbor_id
    using val_t = W;      // weight
    static inline bool comp(key_t a, key_t b) { return a < b; }
  };

  using neighbor_map = pam_map<neighbor_entry>;

  struct clustered_vertex {

    clustered_vertex() {}

    clustered_vertex(uintE vtx_id, Graph_vertex& vertex, const Weights& weights) {
      auto cluster_size = vertex.out_degree();

      auto edges = sequence<edge>::uninitialized(cluster_size);

      auto map_f = [&] (const uintE& u, const uintE& v, const IW& wgh, size_t index) {
        W true_weight = weights.get_weight(u, v, wgh);
        edges[index] = {v, true_weight};
      };

      vertex.out_neighbors().map_with_index(map_f);

      neighbors = neighbor_map(edges);
      auto combine_w = [&] (W l, W r) { return l; };
      neighbors = neighbor_map(edges, combine_w);
    }

    uintE size() {
      return neighbors.size();
    }

    struct Add {
      using T = size_t;
      static T identity() { return 0;}
      static T add(T a, T b) { return a + b;}
    };

    // F is a predicate from neighbor -> bool
    template <class F>
    uintE countNeighbors(uintE id, F& f) {
      auto pred = [&] (edge e) -> size_t {
        return f(id, e.first, e.second);
      };
      neighbor_map::map_reduce(neighbors, pred, Add());
    }

    template <class F>
    void iterate(uintE id, F& f) {
      auto iter = [&] (edge e) -> void {
        f(id, e.first, e.second);
      };
      neighbor_map::foreach_seq(neighbors, iter);
    }

    // A map storing our neighbors + weights.
    neighbor_map neighbors;
  };

  Graph& G;
  Weights& weights;

  uintE n;
  uintE last_cluster_id;
  uintE num_merges_performed;

  sequence<clustered_vertex> clusters;
  sequence<std::pair<uintE, W>> dendrogram;

  uintE degree(uintE v) {
    return clusters[v].size();
  }

  uintE new_cluster_id() {
    uintE ret = last_cluster_id;
    last_cluster_id++;
    return ret;
  }

  uintE unite(uintE a, uintE b) {
    // Identify smaller/larger clusters (will merge smaller -> larger).
    uintE d_a = clusters[a].size();
    uintE d_b = clusters[b].size();
    uintE smaller, larger;
    if (d_a < d_b) {
      smaller = a; larger = b;
    } else {
      larger = a; smaller = b;
    }

    // Merge smaller and larger's neighbors.
    auto smaller_ngh = std::move(clusters[smaller].neighbors);
    auto larger_ngh = std::move(clusters[larger].neighbors);

    // Some sanity asserts, we are merging an edge incident to both after all.
    assert(smaller_ngh.size() > 0);
    assert(larger_ngh.size() > 0);

    // Remove larger's id from smaller, and vice versa.
    auto smaller_keys = neighbor_map::keys(smaller_ngh);

    auto merged = neighbor_map::map_union(
        std::move(smaller_ngh),
        std::move(larger_ngh),
        Weights::linkage);
    clusters[larger].neighbors = std::move(merged);

    return larger;
  }

  clustered_graph(Graph& G, Weights& weights) : G(G), weights(weights) {
    n = G.n;
    last_cluster_id = n;
    num_merges_performed = 0;
    clusters = sequence<clustered_vertex>(n);
    dendrogram = sequence<std::pair<uintE, W>>(2*n - 1, std::make_pair(UINT_E_MAX, W()));

    parallel_for(0, n, [&] (size_t i) {
      auto orig = G.get_vertex(i);
      clusters[i] = clustered_vertex(i, orig, weights);
    });
    debug(std::cout << "Built all vertices" << std::endl;);
  }

};


}  // namespace clustering
}  // namespace gbbs
