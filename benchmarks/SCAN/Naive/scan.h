// A simple implementation of SCAN.
//
// This implementation doesn't attempt to be fast but is useful for checking
// correctness of other SCAN implementations.
#pragma once

#include <functional>
#include <limits>
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
  const auto neighbor_is_epsilon_similar{
    [&](const uintE u, const uintE v, internal::NoWeight) {
      constexpr float
        defaultSimilarity{std::numeric_limits<float>::signaling_NaN()};
      return similarities.find({u, v}, defaultSimilarity) >= epsilon;
    }};
  const pbbs::sequence<bool> core_bitmap(
      num_vertices,
      [&](const size_t i) {
        // `+ 1` to account for open vs. closed neighborhood
        return
          graph->get_vertex(i)
            .countOutNgh(i, neighbor_is_epsilon_similar) + 1 >= mu;
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
      internal::RemoveDuplicates(&core_frontier);
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

  // Cluster all the non-cores by attaching them to the clusters of
  // epsilon-neighbor cores.
  par_for(0, num_vertices, [&](const uintE vertex_id) {
    if (clustering[vertex_id].empty()) {
      auto vertex{graph->get_vertex(vertex_id)};
      const uintE degree{vertex.getOutDegree()};

      sequence<uintE> neighboring_clusters(degree);
      size_t index{0};
      const auto get_neighboring_clusters{
        [&](const uintE v_id, const uintE ngh_id, internal::NoWeight) {
          constexpr float
            kDefaultSimilarity{std::numeric_limits<float>::signaling_NaN()};
          if (core_bitmap[ngh_id] &&
              similarities.find({v_id, ngh_id}, kDefaultSimilarity)
                >= epsilon) {
            neighboring_clusters[index++] = clustering[ngh_id][0];
          }
        }};
      constexpr bool kParallel{false};
      vertex.mapOutNgh(vertex_id, get_neighboring_clusters, kParallel);
      clustering[vertex_id] = pbbs::remove_duplicates_ordered(
          neighboring_clusters.slice(0, index), std::less<uintE>{});
    }
  });

  return clustering;
}

}  // namespace naive_scan
