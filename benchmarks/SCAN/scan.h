#pragma once

#include <cmath>

#include <unordered_set>
#include <utility>
#include <vector>

#include "ligra/macros.h"
#include "pbbslib/parallel.h"
#include "benchmarks/SCAN/symmetric_2d_table.h"

// Index for an undirected graph from which clustering the graph with SCAN is
// quick.
//
// Based off of SCAN index presented in "Efficient Structural Graph Clustering:
// An Index-Based Approach" by Wen et al.
class ScanIndex {
 public:
  template <class Graph>
  explicit ScanIndex(Graph& graph);

 private:
  // Gets the structural similarity between vertices `u` and `v`.
  float GetSimilarity(uintE u, uintE v);
  // Sets the structural similarity between vertices `u` and `v` to be
  // `similarity`.
  void SetSimilarity(uintE u, uintE v, float similarity);

  const uintE num_vertices_;
  // Stores structural similarities between each pair of vertices.
  // TODO(tom.tseng): This is a poor way to store these similarities.
  // A similarity are only non-zero if there's an edge between the pair of
  // vertices. Replace this with a symmetric hash map.
  Symmetric2DTable similarities_;
};

template <class Graph>
ScanIndex::ScanIndex(Graph& graph)
  : num_vertices_{static_cast<uintE>(graph.n)}
  , similarities_{graph.n} {
  using Vertex = typename Graph::vertex;
  using Weight = typename Graph::weight_type;

  // Compute structural similarities between each pair of adjacent vertices.

  std::vector<std::unordered_set<uintE>> adjacency_list{graph.n};
  parallel_for(0, graph.n, [&](const size_t vertex_id) {
    Vertex vertex{graph.get_vertex(vertex_id)};
    auto& neighbors{adjacency_list[vertex_id]};
    neighbors.reserve(vertex.getOutDegree());

    const auto update_adjacency_list = [&neighbors](
        const uintE source_vertex,
        const uintE neighbor_vertex,
        const Weight weight) {
      neighbors.insert(neighbor_vertex);
    };
    const bool kParallel{false};
    vertex.mapOutNgh(vertex_id, update_adjacency_list, kParallel);
  });

  graph.map_edges([&](
        const uintE u,
        const uintE v,
        const Weight weight) {
      const auto& u_neighbors{adjacency_list[u]};
      const auto& v_neighbors{adjacency_list[v]};
      const bool u_neighbors_is_smaller{
        u_neighbors.size() < v_neighbors.size()
      };
      const auto& smaller_neighbor_list{
        u_neighbors_is_smaller? u_neighbors : v_neighbors
      };
      const auto& larger_neighbor_list{
        u_neighbors_is_smaller? v_neighbors : u_neighbors
      };

      uintE num_shared_neighbors = 0;
      for (const uintE neighbor : smaller_neighbor_list) {
        if (larger_neighbor_list.count(neighbor) > 0) {
          num_shared_neighbors++;
        }
      }

      similarities_.SetEntry(u, v,
          num_shared_neighbors /
              (sqrt(u_neighbors.size()) * sqrt(v_neighbors.size())));
  });
}
