// Miscellaneous utility functions for SCAN. May be worth moving out of the
// `SCAN/IndexBased` folder if we add other SCAN implementations.
//
// Here, a SCAN clustering is a partition of the vertices of a graph.
#pragma once

#include <iostream>
#include <string>

#include "ligra/graph.h"
#include "ligra/macros.h"
#include "pbbslib/seq.h"

namespace scan {

// A clustering is given by assigning each vertex a cluster ID.
using Clustering = pbbs::sequence<uintE>;

// Value in `Clustering` for a vertex that does not belong to any cluster.
constexpr uintE kUnclustered{UINT_E_MAX};

// Type of an unclustered vertex.
enum class UnclusteredType {
  kHub,  // Vertex is adjacent to two or more clusters.
  kOutlier  // Vertex is adjacent to at most one cluster.
};
std::ostream& operator<<(std::ostream&, UnclusteredType);

// Converts clustering to a readable string.
//
// We don't write a `operator<<` overload because `Clustering` is an alias for a
// more general type that isn't necessarily printed in the same way.
std::string ClusteringToString(const Clustering& clustering);

// Given a clustering where each cluster ID is in the range `[0,
// clustering->size())`, relabels the clustering so that every cluster ID is in
// the range [0, <number of clusters>) and returns the number of clusters.
size_t CompactClustering(Clustering* clustering);

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
      const uintE v_id, const uintE neighbor_id, pbbslib::empty) {
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

}  // namespace scan
