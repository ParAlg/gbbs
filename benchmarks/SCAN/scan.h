#pragma once

#include "ligra/graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "benchmarks/SCAN/undirected_edge.h"

namespace scan {

namespace internal {

using StructuralSimilarities =
  sparse_table<UndirectedEdge, float, std::hash<UndirectedEdge>>;

template <class Graph>
StructuralSimilarities ComputeStructuralSimilarities(Graph* graph);

}  // namespace internal

// Index for an undirected graph from which clustering the graph with SCAN is
// quick.
//
// Based off of SCAN index presented in "Efficient Structural Graph Clustering:
// An Index-Based Approach" by Wen et al.
class ScanIndex {
 public:
  template <class Graph>
  explicit ScanIndex(Graph* graph);

 private:
  // Stores structural similarities between each pair of adjacent vertices.
  const internal::StructuralSimilarities similarities_;
};

}  // namespace scan
