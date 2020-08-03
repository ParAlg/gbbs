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
#include "pbbslib/integer_sort.h"
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
template <class Graph>
double Modularity(
    Graph* graph,
    const Clustering& clustering);

//////////////
// Internal //
//////////////

namespace internal {

// For each collection of values that have the same key in the input sequence,
// reduces over those values. Returns a sequence `R` such that `R[i]` is the
// reduction result over values with key `i`.
//
// This function's interface matches the interface of `pbbs::collect_reduce`. We
// don't use `pbbs::collect_reduce` because the implementation is broken at the
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
pbbs::sequence<typename Monoid::T> CollectReduce(
    const Seq& seq,
    Key_fn&& get_key,
    Value_fn&& get_value,
    Monoid&& reduce_fn,
    size_t num_keys) {
  using Value = typename Monoid::T;
  pbbs::sequence<size_t> bucketed_indices{
    seq.size(),
    [](const size_t i) { return i; }};
  const auto index_to_key{[&](const size_t i) { return get_key(seq[i]); }};
  integer_sort_inplace(bucketed_indices.slice(), index_to_key);
  pbbs::sequence<size_t> key_offsets{
    pbbs::get_counts(bucketed_indices, index_to_key, num_keys)};
  pbbslib::scan_add_inplace(key_offsets);
  pbbs::sequence<Value> result{
    num_keys,
    [&](const size_t i) {
      const size_t values_start{key_offsets[i]};
      const size_t values_end{
        i + 1 == key_offsets.size() ? seq.size() : key_offsets[i + 1]};
      const auto values{pbbs::delayed_seq<Value>(
          values_end - values_start,
          [&](const size_t j) {
            return get_value(seq[bucketed_indices[values_start + j]]);
          })};
      return pbbslib::reduce(values, reduce_fn);
    }};
  return result;
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

template <class Graph>
double Modularity(
    Graph* graph,
    const Clustering& clustering) {
  const size_t num_edges{graph->m};  // two times the number of undirected edges
  if (num_edges == 0) {
    return 0.0;
  }

  double wm = 0.0;
  pbbs::sequence<double> degs(graph->n, 0.0);
  for (size_t i = 0; i < graph->n; i++) {
    auto v = graph->get_vertex(i);
    const size_t d = v.getOutDegree();
    auto* nghs = v.getOutNeighbors();
    for (size_t j = 0; j < d; j++) {
      degs[i] += std::get<1>(nghs[j]);
    }
    wm += degs[i];
  }
  double pos = 0.0;
  double null = 0.0;
  for (size_t i = 0; i < graph->n; i++) {
    auto v = graph->get_vertex(i);
    const size_t d = v.getOutDegree();
    auto* nghs = v.getOutNeighbors();
    const uintE vc = clustering[i];
    size_t k = 0;
    for (size_t j = 0; j < graph->n; j++) {
      if (vc == clustering[std::get<0>(nghs[j])]) {
        while (k + 1 < d && std::get<0>(nghs[k]) < j) {
          ++k;
        }
        if (std::get<0>(nghs[k]) == j) {
          pos += std::get<1>(nghs[k]);
        }
        null += degs[i] * degs[j];
      }
    }
  }

  return (pos - (null / wm)) / wm;
}

}  // namespace scan
}  // namespace gbbs
