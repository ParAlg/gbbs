// Utility functions for working with graphs that are useful for writing unit
// tests.
#pragma once

#include <utility>
#include <unordered_set>
#include <vector>

#include "gtest/gtest.h"
#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "gbbs/undirected_edge.h"
#include "gbbs/vertex.h"
#include "pbbslib/utilities.h"

namespace gbbs {
namespace graph_test {

namespace internal {  // Internal declarations

template <typename Weight>
void CheckUnweightedNeighbors(
    size_t actual_degree,
    const std::tuple<uintE, Weight>* actual_neighbors,
    const std::vector<uintE>& expected_neighbors);

}  // namespace internal

// Make an undirected, unweighted graph from a list of edges.
symmetric_graph<symmetric_vertex, gbbs::empty> MakeUnweightedSymmetricGraph(
    const uintE num_vertices,
    const std::unordered_set<UndirectedEdge>& edges);

// Check that vertex has `expected_neighbors` as its out-neighbors. Does not
// check edge weights. Ordering matters.
template <class Vertex>
void CheckUnweightedOutNeighbors(
    Vertex& vertex,
    const std::vector<uintE>& expected_neighbors) {
  internal::CheckUnweightedNeighbors(
      vertex.out_degree(), vertex.out_neighbors().neighbors, expected_neighbors);
}

// Check that vertex has `expected_neighbors` as its in-neighbors. Does not
// check edge weights. Ordering matters.
template <class Vertex>
void CheckUnweightedInNeighbors(
    Vertex& vertex,
    const std::vector<uintE>& expected_neighbors) {
  internal::CheckUnweightedNeighbors(
      vertex.in_degree(), vertex.in_neighbors().neighbors, expected_neighbors);
}

namespace internal {  // Internal definitions

// Check that the first `actual_degree` entries of `actual_neighbors` equal
// the entire vector `expected_neighbors`.
template <typename Weight>
void CheckUnweightedNeighbors(
    const size_t actual_degree,
    const std::tuple<uintE, Weight>* const actual_neighbors,
    const std::vector<uintE>& expected_neighbors) {
  EXPECT_EQ(actual_degree, expected_neighbors.size());

  std::vector<uintE> neighbors_vector;
  neighbors_vector.reserve(actual_degree);
  for (size_t i = 0; i < actual_degree; i++) {
    neighbors_vector.emplace_back(std::get<0>(actual_neighbors[i]));
  }

  EXPECT_EQ(neighbors_vector, expected_neighbors);
}

}  // namespace internal

}  // namespace graph_test
}  // namespace gbbs
