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


    clustered_vertex(uintE vtx_id, Graph_vertex& vertex, const Weights& weights, edge* edges) {
      auto cluster_size = vertex.out_degree();
      auto combine_w = [&] (W l, W r) { return l; };
      neighbors = neighbor_map(edges, edges + cluster_size, combine_w);
    }

    clustered_vertex(uintE vtx_id, Graph_vertex& vertex, const Weights& weights) {
      auto cluster_size = vertex.out_degree();

      auto edges = sequence<edge>::uninitialized(cluster_size);

      auto map_f = [&] (const uintE& u, const uintE& v, const IW& wgh, size_t index) {
        W true_weight = weights.get_weight(u, v, wgh);
        edges[index] = {v, true_weight};
      };

      vertex.out_neighbors().map_with_index(map_f);

      auto combine_w = [&] (W l, W r) { return l; };
      neighbors = neighbor_map(edges, combine_w);
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
      return neighbor_map::map_reduce(neighbors, pred, Add());
    }

    template <class F>
    void iterate(uintE id, F& f) {
      auto iter = [&] (edge e) -> void {
        f(id, e.first, e.second);
      };
      neighbor_map::foreach_seq(neighbors, iter);
    }

    uintE neighbor_size() {
      return neighbors.size();
    }

    // Number of vertices contained in this cluster.
    uintE cluster_size() {
      return num_in_cluster;
    }

    // Tracks the last cluster update size.
    uintE staleness;
    // Initially equal to staleness, but cas'd in parallel rounds.
    uintE cas_size;
    // The "current" id of this cluster, updated upon a merge that keeps this cluster active.
    uintE current_id;
    // Number of vertices contained in this cluster.
    uintE num_in_cluster;
    // Active = false iff this cluster is no longer active.
    bool active;
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
    return clusters[v].neighbor_size();
  }

  // Need to implement a unite_merge operation.
  // input: sequence of (u, v) pairs representing that v will merge to u
  //
  // outline:
  // - sort by 2nd component
  // - perform merges from u_1,...,u_k to v.
  // -
  void unite_merge(sequence<std::pair<uintE, uintE>>&& merge_seq) {
    //auto sorted = parlay::sample_sort(make_slice(merge_seq), [&] (const auto& pair) { return pair.second; });

    auto sorted = parlay::sort(make_slice(merge_seq));

    // For each instance, find largest component.
    auto all_starts = parlay::delayed_seq<uintE>(sorted.size(), [&] (size_t i) {
      if ((i == 0) || sorted[i].first != sorted[i-1].first) {
        return (uintE)i;
      }
      return UINT_E_MAX;
    });
    auto starts = parlay::filter(all_starts, [&] (uintE v) { return v != UINT_E_MAX; });

    std::cout << "Number of merge targets = " << starts.size() << std::endl;

    parallel_for(0, starts.size(), [&] (size_t i) {
      uintE our_id = sorted[start].first;
      auto our_size = CG.clusters[our_id].cluster_size();

      size_t start = starts[i];
      size_t end = (i == starts.size() - 1) ? sorted.size() : starts[i+1];
      auto sizes_and_id = parlay::delayed_seq<std::pair<uintE, uintE>>(end - start, [&] (size_t i) {
        uintE vtx_id = sorted[start + i].second;
        return {CG.clusters[vtx_id].cluster_size(), vtx_id};
      });
      std::pair<uintE, uintE> id = std::make_pair((uintE)0, (uintE)UINT_E_MAX);
      auto mon = make_monoid([&] (const auto& l, const auto&, r) { return (l.first < r.first) ? r : l;});
      auto [largest_size, largest_id] = parlay::reduce(sizes_and_id, mon);

      if (our_size < largest_size) {  // relabel
        for (size_t i=start; i<end; i++) {
          sorted[i].first = largest_id;
          if (sorted[i].second == largest_id) {
            sorted[i].second = our_id;
          }
        }

        our_id = largest_id;
        our_size = largest_size;
      }


    });
  }

  clustered_graph(Graph& G, Weights& weights) : G(G), weights(weights) {
    n = G.n;
    last_cluster_id = n;
    num_merges_performed = 0;
    clusters = sequence<clustered_vertex>(n);
    dendrogram = sequence<std::pair<uintE, W>>(2*n - 1, std::make_pair(UINT_E_MAX, W()));

    neighbor_map::reserve(G.m);

    auto edges = sequence<edge>::uninitialized(G.m);
    auto offsets = sequence<size_t>::from_function(G.n, [&] (size_t i) {
      return G.get_vertex(i).out_degree();
    });
    pbbslib::scan_inplace(make_slice(offsets));
    parallel_for(0, n, [&] (size_t i) {
      size_t off = offsets[i];
      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh, size_t j) {
        edges[off + j] = {v, wgh};
      };
      G.get_vertex(i).out_neighbors().map_with_index(map_f);
    });

    parallel_for(0, n, [&] (size_t i) {
      auto orig = G.get_vertex(i);
      size_t off = offsets[i];
      clusters[i] = clustered_vertex(i, orig, weights, edges.begin() + off);
    }, 1);
    std::cout << "Built all vertices" << std::endl;
  }

};


}  // namespace clustering
}  // namespace gbbs
