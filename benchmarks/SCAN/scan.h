#pragma once

#include "ligra/graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "pbbslib/seq.h"
#include "benchmarks/SCAN/undirected_edge.h"

namespace scan {

namespace internal {

using StructuralSimilarities =
  sparse_table<UndirectedEdge, float, std::hash<UndirectedEdge>>;

struct NeighborSimilarity {
  uintE neighbor;
  float similarity;
};

bool operator==(const NeighborSimilarity&, const NeighborSimilarity&);

using NeighborOrder = pbbs::sequence<pbbs::sequence<NeighborSimilarity>>;

template <class Graph>
StructuralSimilarities ComputeStructuralSimilarities(Graph* graph);

template <class Graph>
NeighborOrder
ComputeNeighborOrder(Graph* graph, const StructuralSimilarities& similarities);

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

  // TODO(tom.tseng): Make this class actually store stuff
  // TODO(tom.tseng): add functions for outputting SCAN clusterings of graph
};

}  // namespace scan
