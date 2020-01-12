#define NOTMAIN

#include "benchmarks/SCAN/scan.h"

#include <cmath>
#include <atomic>
#include <unordered_set>
#include <utility>
#include <vector>

#include "ligra/macros.h"
#include "pbbslib/parallel.h"

namespace scan {

namespace internal {

// Compute structural similarities (as defined by SCAN) between each pair of
// adjacent vertices.
//
// The structural similarity between two vertices u and v is
//   (size of intersection of closed neighborhoods of u and v) /
//   (geometric mean of size of closed neighborhoods of u and of v)
// where the closed neighborhood of a vertex x consists of all neighbors of x
// along with x itself.
template <class Graph>
StructuralSimilarities ComputeStructuralSimilarities(Graph* graph) {
  using Vertex = typename Graph::vertex;
  using Weight = typename Graph::weight_type;

  StructuralSimilarities similarities{
    graph->m,
    std::make_pair(UndirectedEdge{UINT_E_MAX, UINT_E_MAX}, 0.0),
    std::hash<UndirectedEdge>{}};

  std::vector<sparse_table<
    uintE, pbbslib::empty, std::function<decltype(pbbslib::hash64_2)>>>
    adjacency_list{graph->n};
  parallel_for(0, graph->n, [&graph, &adjacency_list](const size_t vertex_id) {
    Vertex vertex{graph->get_vertex(vertex_id)};
    auto* neighbors{&adjacency_list[vertex_id]};
    *neighbors = make_sparse_table<
      uintE, pbbslib::empty, std::function<decltype(pbbslib::hash64_2)>>(
        // Adding 1 avoids having small tables completely full
        vertex.getOutDegree() + 1,
        {UINT_E_MAX, pbbslib::empty{}},
        pbbslib::hash64_2);

    const auto update_adjacency_list = [&neighbors](
        const uintE source_vertex,
        const uintE neighbor_vertex,
        const Weight weight) {
      neighbors->insert(std::make_pair(neighbor_vertex, pbbslib::empty{}));
    };
    vertex.mapOutNgh(vertex_id, update_adjacency_list);
  });

  graph->map_edges([&graph, &adjacency_list, &similarities](
        const uintE u_id,
        const uintE v_id,
        const Weight) {
      // Only perform this computation once for each undirected edge
      if (u_id < v_id) {
        Vertex u{graph->get_vertex(u_id)};
        Vertex v{graph->get_vertex(v_id)};
        const auto& u_neighbors{adjacency_list[u_id]};
        const auto& v_neighbors{adjacency_list[v_id]};

        const bool u_degree_is_smaller{u.getOutDegree() < v.getOutDegree()};
        const uintE smaller_degree_vertex_id{u_degree_is_smaller ? u_id : v_id};
        Vertex* smaller_degree_vertex{u_degree_is_smaller ? &u : &v};
        const auto& larger_degree_vertex_neighbors{
          u_degree_is_smaller ? v_neighbors : u_neighbors
        };

        std::atomic<uintE> num_shared_neighbors{0};
        const auto count_shared_neighbors{
          [&larger_degree_vertex_neighbors, &num_shared_neighbors]
          (const uintE, const uintE neighbor, const Weight) {
            if (larger_degree_vertex_neighbors.contains(neighbor)) {
                num_shared_neighbors++;
            }
          }};
        smaller_degree_vertex->mapOutNgh(
            smaller_degree_vertex_id,
            count_shared_neighbors);

        // The neighborhoods we've computed are open neighborhoods -- since
        // structural similarity uses closed neighborhoods, we need to adjust the
        // number and denominator a little.
        similarities.insert({UndirectedEdge{u_id, v_id},
            (num_shared_neighbors + 2) /
                (sqrt(graph->get_vertex(u_id).getOutDegree() + 1) *
                 sqrt(graph->get_vertex(v_id).getOutDegree() + 1))});
      }
  });

  return similarities;
}

template
StructuralSimilarities ComputeStructuralSimilarities(
    symmetric_graph<symmetric_vertex, pbbslib::empty>*);

}  // namespace internal

template <class Graph>
ScanIndex::ScanIndex(Graph* graph)
  : similarities_{internal::ComputeStructuralSimilarities(graph)} {}

template
ScanIndex::ScanIndex(symmetric_graph<symmetric_vertex, pbbslib::empty>*);

}  // namespace scan
