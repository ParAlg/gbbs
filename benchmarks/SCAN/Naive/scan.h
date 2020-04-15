// A simple implementation of SCAN.
//
// This implementation doesn't attempt to be fast but is useful for checking
// correctness of other SCAN implementations.
#pragma once

#include <functional>
#include <tuple>
#include <utility>

#include "benchmarks/SCAN/Naive/scan_helpers.h"
#include "ligra/graph.h"
#include "ligra/ligra.h"
#include "ligra/macros.h"
#include "ligra/vertex_subset.h"
#include "pbbslib/seq.h"
#include "pbbslib/stlalgs.h"
#include "pbbslib/utilities.h"

namespace naive_scan {

// We represent a clustering on an n-vertex graph by an n-length sequence S such
// that S[i] is the IDs of clusters that vertex i is in. Some vertices might not
// be a member of any cluster, and some "border" vertices may be a member of
// many clusters.
using Clustering = pbbs::sequence<pbbs::sequence<uintE>>;

// Compute a SCAN clustering of a graph using SCAN parameters mu and epsilon.
template <template <typename> class VertexTemplate>
Clustering Cluster(
    symmetric_graph<VertexTemplate, pbbslib::empty>* graph,
    const uint64_t mu,
    const float epsilon) {
  const size_t num_vertices{graph->n};
  Clustering clustering(graph->n, pbbs::sequence<uintE>{});

  internal::StructuralSimilarities
    similarities{internal::ComputeStructuralSimilarities(graph)};
  const auto neighbor_has_epsilon_similarity{
    [&](const uintE u, const uintE v, internal::NoWeight) {
      constexpr float defaultSimilarity{-1.0};
      return similarities.find({u, v}, defaultSimilarity) >= epsilon;
    }};
  const pbbs::sequence<bool> core_bitmap(
      num_vertices,
      [&](const size_t i) {
        return
          graph->get_vertex(i)
            .countOutNgh(i, neighbor_has_epsilon_similarity) >= mu;
      });
  const auto is_core{[&](const uintE v) { return core_bitmap[v]; }};
  // Cluster all the cores by running BFS on the more-than-epsilon-similar edges
  // of the core vertices.
  for (uintE root = 0; root < num_vertices; root++) {
    if (!clustering[root].empty()) {
      continue;
    }
    const auto update_clustering_for_cores{[&](const uintE v) {
      clustering[v] = pbbs::sequence<uintE>(1, root);
    }};

    vertexSubset frontier(num_vertices, root);
    while (!frontier.isEmpty()) {
      vertexSubset core_frontier =
        frontier.isDense
        ? vertexFilter(frontier, is_core)
        : vertexFilter2(frontier, is_core);
      vertexMap(core_frontier, update_clustering_for_cores);

      vertexSubset next_frontier{edgeMap(
          *graph,
          core_frontier,
          internal::CoreBFSEdgeMapFunctions{
            similarities, clustering, epsilon})};
      frontier.del();
      core_frontier.del();
      frontier = std::move(next_frontier);
    }
    frontier.del();
  }

  const auto is_core_neighbor{
    [&](const std::tuple<uintE, internal::NoWeight> neighbor) {
      return core_bitmap[std::get<0>(neighbor)];
    }};
  const auto get_core_cluster{
    [&](const std::tuple<uintE, internal::NoWeight> neighbor) {
      return clustering[std::get<0>(neighbor)][0];
    }};
  // Cluster all the non-cores by attaching them to the clusters of neighboring
  // cores.
  par_for(0, num_vertices, [&](const size_t vertex_id) {
    if (clustering[vertex_id].empty()) {
      auto vertex{graph->get_vertex(vertex_id)};
      const uintE degree{vertex.getOutDegree()};
      std::tuple<uintE, internal::NoWeight>*
        neighbors{vertex.getOutNeighbors()};
      const sequence<std::tuple<uintE, internal::NoWeight>> core_neighbors{
        filter(
            pbbs::range<std::tuple<uintE, internal::NoWeight>*>{
              neighbors, neighbors + degree},
            is_core_neighbor)};
      const sequence<uintE> neighboring_clusters{
        pbbs::map<uintE>(core_neighbors, get_core_cluster)};
      clustering[vertex_id] = pbbs::remove_duplicates_ordered(
          neighboring_clusters, std::less<uintE>{});
    }
  });

  return clustering;
}

}  // namespace naive_scan
