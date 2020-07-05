// Miscellaneous utility functions for SCAN. May be worth moving out of the
// `SCAN/IndexBased` folder if we add other SCAN implementations.
//
// Here, a SCAN clustering is a partition of the vertices of a graph.
#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <utility>

#include "gbbs/bridge.h"
#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "pbbslib/monoid.h"
#include "pbbslib/seq.h"

namespace gbbs {
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
    const Clustering& clustering, Vertex vertex, uintE vertex_id);

// Quality measure of clustering based on the difference of the edge density of
// each cluster compared to the expected edge density of the cluster given a
// random graph with the same degree distribution.
//
// The modularity is at most 1. It may be negative if a clustering is "worse
// than random."
//
// This implementation treats each unclustered vertex as if it is in its own
// cluster.
//
// Further reading:
// - Wikipedia page on "Modularity (networks)"
// - Section 3.3.2 of "Community detection in graphs" by Santo Fortunato (2009).
template <template <typename> class VertexTemplate>
double Modularity(
    symmetric_graph<VertexTemplate, pbbslib::empty>* graph,
    const Clustering& clustering);

//////////////
// Internal //
//////////////

namespace internal {


template <class A, class B, class C, class D, class E>
const pbbs::sequence<uintT> CollectReduce(A, B, C, D, E)  {
  // TODO implement this --- this used to be collect_reduce from
  // collect_reduce.h, but that's buggy
  assert(false);  // TODO
  return {};
}

}  // namespace internal

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

template <template <typename> class VertexTemplate>
double Modularity(
    symmetric_graph<VertexTemplate, pbbslib::empty>* graph,
    const Clustering& clustering) {
  const size_t num_edges{graph->m};  // two times the number of undirected edges
  if (num_edges == 0) {
    return 0.0;
  }
  const size_t num_vertices{graph->n};
  const size_t num_clusters{
    1 +
    pbbslib::reduce_max(pbbs::delayed_seq<uintE>(
      num_vertices,
      [&](const size_t vertex_id) {
        const uintE cluster_id{clustering[vertex_id]};
        return cluster_id == kUnclustered ? 0 : cluster_id;
      }))};

  // Fraction of edges that fall within a cluster.
  const double intracluster_edge_proportion{
    static_cast<double>(pbbslib::reduce_add(pbbs::delayed_seq<uintT>(
      num_vertices,
      [&](const size_t vertex_id) -> uintT {
        const uintE cluster_id{clustering[vertex_id]};
        if (cluster_id == kUnclustered) {
          return 0;
        }
        const auto is_same_cluster{
          [&](const uintE v_id, const uintE ngh_id, pbbslib::empty) {
            return clustering[ngh_id] == cluster_id;
          }};
        return graph->get_vertex(vertex_id)
          .countOutNgh(vertex_id, is_same_cluster);
      }))) / num_edges};

  const auto degrees_split_result{
    pbbs::split_two(
      pbbs::delayed_seq<std::pair<uintE, uintT>>(
        num_vertices,
        [&](const size_t i) {
          return std::make_pair(clustering[i], graph->get_vertex(i).degree);
        }),
      pbbs::delayed_seq<bool>(
        num_vertices,
        [&](const size_t i) { return clustering[i] != kUnclustered; }))};
  // <cluster id, degree> of each vertex, with unclustered vertices at the start
  // of the list
  const auto& degrees{degrees_split_result.first};
  const auto& num_unclustered_vertices{degrees_split_result.second};
  // degrees_by_cluster[i] == sum of degrees over vertices in cluster i
  const pbbs::sequence<uintT> degrees_by_cluster{
    internal::CollectReduce(
      degrees.slice(num_unclustered_vertices, degrees.size()),
      [&](const std::pair<uintE, uintT> p) { return p.first; },
      [&](const std::pair<uintE, uintT> p) { return p.second; },
      pbbs::addm<uintT>{},
      num_clusters)};
  // Approximately the fraction of edges that fall within a cluster for a random
  // graph with the same degree distribution:
  //   sum((sum(degree) for each vertex in cluster) / (2 * <number of edges>)
  //       for each cluster in graph)
  const double null_intracluster_proportion{
     pbbslib::reduce_add(pbbs::delayed_seq<double>(
       num_clusters,
       [&](const size_t cluster_id) {
         return std::pow(
             static_cast<double>(degrees_by_cluster[cluster_id]) / num_edges,
             2.0);
       })) +
     // special case for unclustered vertices
     pbbslib::reduce_add(pbbs::delayed_seq<double>(
       num_unclustered_vertices,
       [&](const size_t i) {
         return std::pow(
             static_cast<double>(degrees[i].second) / num_edges, 2.0);
       }))};

  return intracluster_edge_proportion - null_intracluster_proportion;
}

}  // namespace scan
}  // namespace gbbs
