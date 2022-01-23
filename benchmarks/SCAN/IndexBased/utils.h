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

namespace gbbs {
namespace scan {

// A clustering is given by assigning each vertex a cluster ID.
using Clustering = sequence<uintE>;

// Value in `Clustering` for a vertex that does not belong to any cluster.
constexpr uintE kUnclustered{UINT_E_MAX};

// Type of an unclustered vertex.
enum class UnclusteredType {
  kHub,     // Vertex is adjacent to two or more clusters.
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
UnclusteredType DetermineUnclusteredType(const Clustering& clustering,
                                         Vertex vertex, uintE vertex_id);

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
template <template <typename> class VertexTemplate, typename Weight>
double Modularity(symmetric_graph<VertexTemplate, Weight>* graph,
                  const Clustering& clustering);

//////////////
// Internal //
//////////////

namespace internal {

// For each collection of values that have the same key in the input sequence,
// reduces over those values. Returns a sequence `R` such that `R[i]` is the
// reduction result over values with key `i`.
//
// This function's interface matches the interface of `collect_reduce`. We
// don't use `collect_reduce` because the implementation is broken at the
// time of writing this comment.
//
// Arguments
// ---------
// seq: Sequence<Element>
//   Input sequence of key-value pairs
// get_key: Element -> size_t
//   Function that gets the key of an element from `seq`.
// get_value: Element -> Value
//   Function that gets the value of an element from `seq`
// reduce: Monoid<Value>
//   Monoid with a reduction function over values
// num_keys:
//   Must be such that keys returned by `get_key on elements of `seq` are in the
//   range [0, `num_keys`).
template <class Seq, class Key_fn, class Value_fn, class Monoid>
sequence<typename std::remove_reference_t<Monoid>::T> CollectReduce(
    const Seq& seq, Key_fn&& get_key, Value_fn&& get_value, Monoid&& reduce_fn,
    size_t num_keys) {
  using Value = typename std::remove_reference_t<Monoid>::T;
  if (seq.empty()) {
    return sequence<Value>(num_keys, reduce_fn.identity);
  }
  sequence<size_t> bucketed_indices = sequence<size_t>::from_function(
      seq.size(), [](const size_t i) { return i; });
  const auto index_to_key{[&](const size_t i) { return get_key(seq[i]); }};
  integer_sort_inplace(make_slice(bucketed_indices), index_to_key);
  sequence<size_t> key_offsets{
      parlay::get_counts(make_slice(bucketed_indices), index_to_key, num_keys)};
  parlay::scan_inplace(key_offsets);
  sequence<Value> result =
      sequence<Value>::from_function(num_keys, [&](const size_t i) {
        const size_t values_start{key_offsets[i]};
        const size_t values_end{
            i + 1 == key_offsets.size() ? seq.size() : key_offsets[i + 1]};
        const auto values{parlay::delayed_seq<Value>(
            values_end - values_start, [&](const size_t j) {
              return get_value(seq[bucketed_indices[values_start + j]]);
            })};
        return parlay::reduce(values, std::forward<Monoid>(reduce_fn));
      });
  return result;
}

}  // namespace internal

template <class Vertex>
UnclusteredType DetermineUnclusteredType(const Clustering& clustering,
                                         Vertex vertex, uintE vertex_id) {
  bool is_hub{false};
  // `candidate_cluster` holds a cluster ID that vertex i is adjacent to, or
  // `kUnclustered` before such a cluster is found.
  uintE candidate_cluster{kUnclustered};
  const auto check_neighbor{
      [&](const uintE v_id, const uintE neighbor_id, gbbs::empty) {
        const uintE neighbor_cluster{clustering[neighbor_id]};
        // If `candidate_cluster` is at its default value of `kUnclustered`,
        // assign
        // `neighbor_cluster` to it. Otherwise, if it has a value differing from
        // `neighbor_cluster`, then we know vertex i is adjacent to multiple
        // clusters and is a hub.
        if (neighbor_cluster != kUnclustered &&
            !(candidate_cluster == UINT_E_MAX &&
              gbbs::atomic_compare_and_swap(&candidate_cluster, UINT_E_MAX,
                                            neighbor_cluster)) &&
            candidate_cluster != neighbor_cluster && !is_hub) {
          is_hub = true;
        }
      }};
  vertex.out_neighbors().map(check_neighbor);
  return is_hub ? UnclusteredType::kHub : UnclusteredType::kOutlier;
}

template <template <typename> class VertexTemplate, typename Weight>
double Modularity(symmetric_graph<VertexTemplate, Weight>* graph,
                  const Clustering& clustering) {
  const size_t num_edges{graph->m};  // two times the number of undirected edges
  if (num_edges == 0) {
    return 0.0;
  }
  const size_t num_vertices{graph->n};
  const size_t num_clusters{parlay::reduce_max(
      parlay::delayed_seq<uintE>(num_vertices, [&](const size_t vertex_id) {
        const uintE cluster_id{clustering[vertex_id]};
        return cluster_id == kUnclustered ? 0 : cluster_id + 1;
      }))};

  if
    constexpr(std::is_same<Weight, gbbs::empty>::value) {
      // unweighted case

      // Fraction of edges that fall within a cluster.
      const double intracluster_edge_proportion{
          static_cast<double>(parlay::reduce(parlay::delayed_seq<uintT>(
              num_vertices,
              [&](const size_t vertex_id) -> uintT {
                const uintE cluster_id{clustering[vertex_id]};
                if (cluster_id == kUnclustered) {
                  return 0;
                }
                const auto is_same_cluster{
                    [&](const uintE v_id, const uintE ngh_id, gbbs::empty) {
                      return clustering[ngh_id] == cluster_id;
                    }};
                return graph->get_vertex(vertex_id).out_neighbors().count(
                    is_same_cluster);
              }))) /
          num_edges};

      const auto degrees_split_result{parlay::split_two(
          parlay::delayed_seq<std::pair<uintE, uintT>>(
              num_vertices,
              [&](const size_t i) {
                return std::make_pair(clustering[i],
                                      graph->get_vertex(i).degree);
              }),
          parlay::delayed_seq<bool>(num_vertices, [&](const size_t i) {
            return clustering[i] != kUnclustered;
          }))};
      // <cluster id, degree> of each vertex, with unclustered vertices at the
      // start of the list
      const auto& clusters_and_degrees{degrees_split_result.first};
      const auto& num_unclustered_vertices{degrees_split_result.second};
      // degrees_by_cluster[i] == sum of degrees over vertices in cluster i
      const sequence<uintT> degrees_by_cluster{internal::CollectReduce(
          clusters_and_degrees.cut(num_unclustered_vertices,
                                   clusters_and_degrees.size()),
          [&](const std::pair<uintE, uintT> p) { return p.first; },
          [&](const std::pair<uintE, uintT> p) { return p.second; },
          parlay::addm<uintT>{}, num_clusters)};
      // Approximately the fraction of edges that fall within a cluster for a
      // random graph with the same degree distribution:
      //   sum((sum(degree) for each vertex in cluster) / (2 * <number of
      //   edges>)
      //       for each cluster in graph)
      const double null_intracluster_proportion{
          parlay::reduce(parlay::delayed_seq<double>(
              num_clusters,
              [&](const size_t cluster_id) {
                return std::pow(
                    static_cast<double>(degrees_by_cluster[cluster_id]) /
                        num_edges,
                    2.0);
              })) +
          // special case for unclustered vertices
          parlay::reduce(parlay::delayed_seq<double>(
              num_unclustered_vertices, [&](const size_t i) {
                return std::pow(
                    static_cast<double>(clusters_and_degrees[i].second) /
                        num_edges,
                    2.0);
              }))};

      return intracluster_edge_proportion - null_intracluster_proportion;
    }
  else {  // weighted case
    constexpr auto get_weight{
        [](uintE, uintE, const Weight weight) { return weight; }};
    const parlay::addm<double> add_weights{};
    // weighted_degrees[i] = sum of weights of incident edges on vertex i
    const auto weighted_degrees{
        sequence<double>::from_function(graph->n, [&](const size_t i) {
          return graph->get_vertex(i).out_neighbors().reduce(get_weight,
                                                             add_weights);
        })};
    // 2 * sum of all edge weights
    const double total_weight{parlay::reduce(weighted_degrees)};

    // Fraction of edge weight that falls within a cluster.
    const double intracluster_edge_proportion{
        parlay::reduce(parlay::delayed_seq<double>(
            num_vertices,
            [&](const size_t vertex_id) {
              const uintE cluster_id{clustering[vertex_id]};
              if (cluster_id == kUnclustered) {
                return 0.0;
              }
              const auto get_intracluster_weight{
                  [&](const uintE v_id, const uintE ngh_id, const Weight w) {
                    return clustering[ngh_id] == cluster_id ? w : 0.0;
                  }};
              return graph->get_vertex(vertex_id).out_neighbors().reduce(
                  get_intracluster_weight, add_weights);
            })) /
        total_weight};

    const auto degrees_split_result{parlay::split_two(
        parlay::delayed_seq<std::pair<uintE, double>>(
            num_vertices,
            [&](const size_t i) {
              return std::make_pair(clustering[i], weighted_degrees[i]);
            }),
        parlay::delayed_seq<bool>(num_vertices, [&](const size_t i) {
          return clustering[i] != kUnclustered;
        }))};
    // <cluster id, degree> of each vertex, with unclustered vertices at the
    // start of the list
    const auto& clusters_and_degrees{degrees_split_result.first};
    const auto& num_unclustered_vertices{degrees_split_result.second};
    // degrees_by_cluster[i] == sum of weighted degrees over vertices in cluster
    // i
    const sequence<double> degrees_by_cluster{internal::CollectReduce(
        clusters_and_degrees.cut(num_unclustered_vertices,
                                 clusters_and_degrees.size()),
        [&](const std::pair<uintE, double> p) { return p.first; },
        [&](const std::pair<uintE, double> p) { return p.second; }, add_weights,
        num_clusters)};
    // Approximately the fraction edge weight that falls within a cluster for a
    // random graph with the same degree distribution.
    const double null_intracluster_proportion{
        parlay::reduce(parlay::delayed_seq<double>(
            num_clusters,
            [&](const size_t cluster_id) {
              return std::pow(degrees_by_cluster[cluster_id] / total_weight,
                              2.0);
            })) +
        // special case for unclustered vertices
        parlay::reduce(parlay::delayed_seq<double>(
            num_unclustered_vertices, [&](const size_t i) {
              return std::pow(clusters_and_degrees[i].second / total_weight,
                              2.0);
            }))};

    return intracluster_edge_proportion - null_intracluster_proportion;
  }
}

}  // namespace scan
}  // namespace gbbs
