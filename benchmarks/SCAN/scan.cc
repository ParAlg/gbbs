#define NOTMAIN

#include "benchmarks/SCAN/scan.h"

#include <cmath>
#include <atomic>
#include <unordered_set>
#include <utility>
#include <vector>

#include "ligra/macros.h"
#include "pbbslib/parallel.h"

template <class Graph>
ScanIndex::ScanIndex(Graph* graph)
  : similarities_{
      graph->m,
      std::make_pair(UndirectedEdge{UINT_E_MAX, UINT_E_MAX}, 0.0),
      HashUndirectedEdge} {
  using Vertex = typename Graph::vertex;
  using Weight = typename Graph::weight_type;

  // Compute structural similarities between each pair of adjacent vertices.

  std::vector<sparse_table<
    uintE, pbbslib::empty, std::function<decltype(pbbslib::hash64_2)>>>
    adjacency_list{graph->n};
  parallel_for(0, graph->n, [&](const size_t vertex_id) {
    Vertex vertex{graph->get_vertex(vertex_id)};
    auto& neighbors{adjacency_list[vertex_id]};
    neighbors = make_sparse_table<
      uintE, pbbslib::empty, std::function<decltype(pbbslib::hash64_2)>>(
        // Adding 1 avoids having small tables completely full
        vertex.getOutDegree() + 1,
        {UINT_E_MAX, pbbslib::empty{}},
        pbbslib::hash64_2);

    const auto update_adjacency_list = [&neighbors](
        const uintE source_vertex,
        const uintE neighbor_vertex,
        const Weight weight) {
      neighbors.insert(std::make_pair(neighbor_vertex, pbbslib::empty{}));
    };
    const bool kParallel{false};
    vertex.mapOutNgh(vertex_id, update_adjacency_list, kParallel);
  });

  graph->map_edges([&](
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

      std::atomic<uintE> num_shared_neighbors{0};
      smaller_neighbor_list.map(
        [&](const std::pair<uintE, pbbslib::empty>& kv) {
          if (larger_neighbor_list.contains(kv.first)) {
              num_shared_neighbors++;
          }
      });

      similarities_.insert({UndirectedEdge{u, v},
          num_shared_neighbors /
              (sqrt(u_neighbors.size()) * sqrt(v_neighbors.size()))});
  });
}

template
ScanIndex::ScanIndex(symmetric_graph<symmetric_vertex, pbbslib::empty>*);
