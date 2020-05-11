#pragma once

#include <string>

#include "benchmarks/SCAN/IndexBased/scan_helpers.h"
#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "pbbslib/seq.h"
#include "pbbslib/utilities.h"

namespace indexed_scan {

using Clustering = pbbs::sequence<uintE>;

// Value in `Clustering` for a vertex that does not belong to any cluster.
constexpr uintE kUnclustered{UINT_E_MAX};

// Type of an unclustered vertex.
enum class UnclusteredType {
  kHub,  // Vertex is adjacent to two or more clusters.
  kOutlier  // Vertex is adjacent to at most one cluster.
};
std::ostream& operator<<(std::ostream&, UnclusteredType);

// Index for an undirected graph from which clustering the graph with SCAN is
// quick, though index construction may be expensive.
class Index {
 public:
  // Constructor.
  //
  // The neighbor lists for each vertex in the graph must be sorted by ascending
  // neighbor ID.
  template <template <typename> class VertexTemplate>
  explicit Index(
      symmetric_graph<VertexTemplate, pbbslib::empty>* graph)
    : num_vertices_{graph->n}
    , neighbor_order_{graph}
    , core_order_{neighbor_order_} {}

  // Compute a SCAN clustering of the indexed graph using SCAN parameters
  // mu and epsilon.
  //
  // Those who are familiar with SCAN may know that some "border" vertices of
  // clusters can belong to multiple clusters at once. This implementation picks
  // an arbitrary choice of a single cluster assignment for those vertices.
  //
  // Arguments:
  //   epsilon
  //     A threshold value on the "similarity" between adjacent vertices based
  //     on how much they share neighbors. Increasing this makes finer-grained,
  //     smaller clusters.
  //   mu
  //     How many neighbors a vertex needs to be epsilon-similar to in order to
  //     be considered a "core" vertex from which a cluster is grown.
  //     Increasing this increases the minimum cluster size and but also makes
  //     large clusters less likely to appear.
  //
  // Returns:
  //   `graph->n`-length sequence S where S[i] is the cluster ID of vertex i or
  //   is `kUnclustered` if vertex i does not belong to any cluster. The cluster
  //   IDs will be in the range `[0, graph->n)` but will not necessarily be
  //   contiguous.
  Clustering Cluster(uint64_t mu, float epsilon) const;

 private:
  const size_t num_vertices_;
  const internal::NeighborOrder neighbor_order_;
  const internal::CoreOrder core_order_;
};

// Relabels `clustering` so that every cluster ID is in the range
// [0, <number_of_clusters>) and returns the number of clusters.
size_t CleanClustering(Clustering* clustering);

// Converts clustering to a readable string.
//
// We don't write a `operator<<` overload because `Clustering` is an alias for a
// more general type that isn't necessarily printed in the same way.
std::string ClusteringToString(const Clustering& clustering);

// Determines what type of unclustered vertex the input vertex is.
//
// Arguments:
//   clustering
//     A clustering of the graph `graph`.
//  vertex
//    The vertex of interest.
//  vertex_id
//     Vertex ID of the vertex of interest. `clustering[vertex_id]` should be
//     `kUnclustered`.
template <class Vertex>
UnclusteredType DetermineUnclusteredType(
    const Clustering& clustering, Vertex vertex, uintE vertex_id) {
  bool is_hub{false};
  // `candidate_cluster` holds a cluster ID that vertex i is adjacent to, or
  // `kUnclustered` before such a cluster is found.
  uintE candidate_cluster{kUnclustered};
  const auto check_neighbor{[&](
      const uintE v_id, const uintE neighbor_id, internal::NoWeight) {
    const uintE neighbor_cluster{clustering[neighbor_id]};
    // If `candidate_cluster` is at its default value of `kUnclustered`, assign
    // `neighbor_cluster` to it. Otherwise, if it has a value differing from
    // `neighbor_cluster`, then we know vertex i is adjacent to multiple
    // clusters and is a hub.
    if (neighbor_cluster != kUnclustered &&
        !(candidate_cluster == UINT_E_MAX &&
          pbbs::atomic_compare_and_swap(
                  &candidate_cluster, UINT_E_MAX, neighbor_cluster)) &&
        candidate_cluster != neighbor_cluster && !is_hub) {
      is_hub = true;
    }
  }};
  vertex.mapOutNgh(vertex_id, check_neighbor);
  return is_hub ? UnclusteredType::kHub : UnclusteredType::kOutlier;
}

}  // namespace indexed_scan
