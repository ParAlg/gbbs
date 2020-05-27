// Miscellaneous utility functions for SCAN. May be worth moving out of the
// `SCAN/IndexBased` folder if we add other SCAN implementations.
//
// Here, a SCAN clustering is a partition of the vertices of a graph.
#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <utility>

#include "ligra/bridge.h"
#include "ligra/graph.h"
#include "ligra/macros.h"
#include "pbbslib/integer_sort.h"
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

// Quality measure of clustering based on the difference of the edge density of
// each cluster compared to the expected edge density of the cluster given a
// random graph with the same degree distribution. The modularity falls in the
// range [-1/2, 1].
//
// Further reading:
// - Wikipedia page on "Modularity (networks)"
// - Section 3.3.2 of "Community detection in graphs" by Santo Fortunato (2009).
template <template <typename> class VertexTemplate>
double Modularity(
    symmetric_graph<VertexTemplate, pbbslib::empty>* graph,
    const Clustering& clustering,
    const uintE max_cluster_id) {
  const size_t num_edges{graph->m};  // two times the number of undirected edges
  if (num_edges == 0) {
    return 0.0;
  }
  const size_t num_vertices{graph->n};
  const size_t num_clusters{max_cluster_id + 1};

  // Fraction of edges that fall within a cluster.
  const double intracluster_edge_proportion{
    static_cast<double>(pbbslib::reduce_add(pbbs::delayed_seq<uintT>(
      num_vertices,
      [&](const size_t vertex_id) {
        const uintE cluster_id{clustering[vertex_id]};
        const auto is_same_cluster{
          [&](const uintE v_id, const uintE ngh_id, pbbslib::empty) {
            return clustering[ngh_id] == cluster_id;
          }};
        return graph->get_vertex(vertex_id)
          .countOutNgh(vertex_id, is_same_cluster);
      }))) / num_edges};

  constexpr auto get_first{
    [](const std::pair<uintE, uintE> p) { return p.first; }};
  // <cluster ID, vertex id> pairs, bucketed by cluster ID
  const pbbs::sequence<std::pair<uintE, uintE>> vertices_by_cluster{
    pbbs::integer_sort(
      pbbs::delayed_seq<std::pair<uintE, uintE>>(
        num_vertices,
        [&](const size_t i) {
          return std::make_pair(clustering[i], i);
        }),
      get_first)};
  const pbbs::sequence<uintE> cluster_offsets{
    pbbs::get_counts<uintE>(vertices_by_cluster, get_first, num_clusters)
  };
  // Fraction of edges that fall within a cluster for a random graph with the
  // same degree distribution.
  const double null_intracluster_proportion{
    pbbslib::reduce_add(pbbs::delayed_seq<double>(
      num_clusters,
      [&](const size_t cluster_id) {
        const uintE cluster_start{cluster_offsets[cluster_id]};
        const uintT cluster_degree{pbbslib::reduce_add(pbbs::delayed_seq<uintT>(
          cluster_offsets[cluster_id + 1] - cluster_start,
          [&](const size_t i) {
            return graph->get_vertex(
              vertices_by_cluster[cluster_start + i].second).degree;
          }))};
        return std::pow(static_cast<double>(cluster_degree) / num_edges, 2.0);
      }))};

  // TODO(tomtseng): deal with kUnclustered, document in comment.
  // either remove unclustered vertices, or put them in a separate cluster.
  // see how DBScan comparison papers deal with this

  return intracluster_edge_proportion - null_intracluster_proportion;
}

}  // namespace scan
