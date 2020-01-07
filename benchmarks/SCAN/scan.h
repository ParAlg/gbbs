#pragma once

#include "ligra/graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "benchmarks/SCAN/undirected_edge.h"

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
  sparse_table<
    UndirectedEdge, float, std::function<decltype(HashUndirectedEdge)>>
    similarities_;
};

template
ScanIndex::ScanIndex(symmetric_graph<symmetric_vertex, pbbslib::empty>*);
