// Helper functions for SCAN logic that shouldn't clutter the main SCAN
// header file.
#pragma once

#include <utility>

#include "ligra/graph.h"
#include "ligra/macros.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/undirected_edge.h"
#include "ligra/vertex_subset.h"
#include "pbbslib/seq.h"
#include "pbbslib/utilities.h"

namespace naive_scan {

using Clustering = pbbs::sequence<pbbs::sequence<uintE>>;

namespace internal {  // internal declarations

using StructuralSimilarities =
  sparse_table<UndirectedEdge, float, std::hash<UndirectedEdge>>;
using NoWeight = pbbslib::empty;

class CoreBFSEdgeMapFunctions {
 public:
  CoreBFSEdgeMapFunctions(
      const StructuralSimilarities& similarities,
      const Clustering& current_clustering,
      float epsilon);

  bool update(uintE, uintE, NoWeight) const;
  bool updateAtomic(uintE, uintE, NoWeight) const;
  bool cond(uintE) const;

 private:
  const StructuralSimilarities& similarities_;
  const Clustering& current_clustering_;
  float epsilon_;
};

// Compute structural similarities (as defined by SCAN) between each pair of
// adjacent vertices.
//
// The structural similarity between two vertices u and v is
//   (size of intersection of closed neighborhoods of u and v) /
//   (geometric mean of size of closed neighborhoods of u and of v)
// where the closed neighborhood of a vertex x consists of all neighbors of x
// along with x itself.
template <template <typename WeightType> class VertexType>
StructuralSimilarities
ComputeStructuralSimilarities(symmetric_graph<VertexType, NoWeight>* graph) {
  using Vertex = VertexType<NoWeight>;
  using VertexSet =
    sparse_table<uintE, pbbslib::empty, decltype(&pbbslib::hash64_2)>;

  StructuralSimilarities similarities{
    graph->m,
    std::make_pair(UndirectedEdge{UINT_E_MAX, UINT_E_MAX}, 0.0),
    std::hash<UndirectedEdge>{}};

  pbbs::sequence<VertexSet> adjacency_list{
    pbbs::sequence<VertexSet>::no_init(graph->n)};

  parallel_for(0, graph->n, [&graph, &adjacency_list](const size_t vertex_id) {
    Vertex vertex{graph->get_vertex(vertex_id)};
    auto* const neighbors{&adjacency_list[vertex_id]};
    *neighbors = VertexSet{
      vertex.getOutDegree(),
      {UINT_E_MAX, internal::NoWeight{}},
      pbbslib::hash64_2};
    const auto update_adjacency_list{[&neighbors](
        uintE, const uintE neighbor, NoWeight) {
      neighbors->insert({neighbor, pbbslib::empty{}});
    }};
    vertex.mapOutNgh(vertex_id, update_adjacency_list);
  });

  graph->map_edges([&graph, &adjacency_list, &similarities](
        const uintE u_id, const uintE v_id, NoWeight) {
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

        const auto is_shared_neighbor{
          [&](uintE, const uintE neighbor, NoWeight) {
            return larger_degree_vertex_neighbors.contains(neighbor);
          }};
        const auto addition_monoid{pbbslib::addm<size_t>()};
        const size_t num_shared_neighbors{
          smaller_degree_vertex->template reduceOutNgh<size_t>(
              smaller_degree_vertex_id,
              is_shared_neighbor,
              addition_monoid)};

        // The neighborhoods we've computed are open neighborhoods -- since
        // structural similarity uses closed neighborhoods, we need to adjust
        // the numbers slightly, hence the `+ 1` and `+ 2` terms.
        similarities.insert(
            {UndirectedEdge{u_id, v_id},
            (num_shared_neighbors + 2) /
                (sqrt(u.getOutDegree() + 1) * sqrt(v.getOutDegree() + 1))});
      }
  });
  return similarities;
}

// Remove duplicate vertices from a vertex subset.
void RemoveDuplicates(vertexSubset* vertex_subset);

}  // namespace internal

}  // namespace naive_scan
