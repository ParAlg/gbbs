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
  // Vertex ID.
  uintE neighbor;
  // Similarity of neighbor vertex to some original reference vertex.
  float similarity;
};

bool operator==(const NeighborSimilarity&, const NeighborSimilarity&);

std::ostream& operator<<(std::ostream& os, const NeighborSimilarity&);

using NeighborOrder = pbbs::sequence<pbbs::sequence<NeighborSimilarity>>;

struct CoreThreshold {
  uintE vertex_id;
  // Maximum value of the SCAN parameter epsilon for which `vertex_id` is a core
  // vertex (given some fixed reference value for SCAN parameter mu).
  float threshold;
};

bool operator==(const CoreThreshold&, const CoreThreshold&);

std::ostream& operator<<(std::ostream& os, const CoreThreshold&);

using CoreOrder = pbbs::sequence<pbbs::sequence<CoreThreshold>>;

template <class Graph>
StructuralSimilarities ComputeStructuralSimilarities(Graph* graph);

template <class Graph>
NeighborOrder
ComputeNeighborOrder(Graph* graph, const StructuralSimilarities& similarities);

CoreOrder ComputeCoreOrder(const NeighborOrder& neighbor_order);

}  // namespace internal

// Clustering resulting from running SCAN.
struct Clustering {
  // A list of clusters; `clusters[i]` is a list of vertices in the i-th
  // cluster.
  // Clusters are not necessarily disjoint, and they do not necessarily cover
  // the entire graph.
  pbbs::sequence<pbbs::sequence<uintE>> clusters;
  // Vertices not in any cluster but is adjacent to at least two clusters.
  pbbs::sequence<uintE> hubs;
  // Vertices not in any cluster and adjacent to at most one cluster.
  pbbs::sequence<uintE> outliers;
};

// Index for an undirected graph from which clustering the graph with SCAN is
// quick.
//
// Based off of SCAN index presented in "Efficient Structural Graph Clustering:
// An Index-Based Approach" by Wen et al.
class ScanIndex {
 public:
  template <class Graph>
  explicit ScanIndex(Graph* graph);

  // Compute a SCAN clustering of the indexed graph using SCAN parameters
  // epsilon and mu.
  //
  // Explanation of parameters:
  // - epsilon: a threshold value on the "similarity" between adjacent vertices
  //   based on how much they share neighbors. Increasing this makes
  //   finer-grained, smaller clusters.
  // - mu: how many neighbors a vertex needs to be epsilon-similar to in order
  //   to be considered a "core" vertex from which a cluster is grown.
  //   Increasing this increases the minimum cluster size and but also makes
  //   large clusters less likely to appear.
  Clustering Cluster(float epsilon, uint64_t mu) const;

 private:
  const size_t num_vertices;
  const internal::NeighborOrder neighbor_order;
  const internal::CoreOrder core_order;
};

}  // namespace scan
