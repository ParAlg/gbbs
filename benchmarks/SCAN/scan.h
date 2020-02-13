#pragma once

#include <variant>

#include "ligra/graph.h"
#include "ligra/vertex.h"
#include "ligra/pbbslib/sparse_table.h"
#include "pbbslib/seq.h"
#include "benchmarks/SCAN/undirected_edge.h"

namespace scan {

namespace internal {

struct NeighborSimilarity {
  // Vertex ID.
  uintE neighbor;
  // Similarity of neighbor vertex to some original reference vertex.
  float similarity;
};

using NeighborOrder = pbbs::sequence<pbbs::sequence<NeighborSimilarity>>;

struct CoreThreshold {
  uintE vertex_id;
  // Maximum value of the SCAN parameter epsilon for which `vertex_id` is a core
  // vertex (given some fixed reference value for SCAN parameter mu).
  float threshold;
};

class CoreOrder {
 public:
  explicit CoreOrder(const NeighborOrder& neighbor_order);

  uint64_t MaxMu() const;
  pbbs::sequence<uintE> GetCores(uint64_t mu, float epsilon) const;

 private:
  const size_t num_vertices_;
  pbbs::sequence<pbbs::sequence<CoreThreshold>> order_{};
};

}  // namespace internal

// Represents a vertex that is in at least one SCAN cluster.
struct ClusterMember {
  // IDs of the clusters that the vertex is in. Vertices at the peripheries of
  // clusters may belong to multiple clusters.
  pbbs::sequence<uintE> clusters;
};
bool operator==(const ClusterMember&, const ClusterMember&);

// Represents a vertex that is not in a cluster but is adjacent to at least two
// clusters.
struct Hub{};
bool operator==(const Hub&, const Hub&);

// Represents a vertex that is not in a cluster and is adjacent to at most one
// cluster.
struct Outlier{};
bool operator==(const Outlier&, const Outlier&);

// Possible results for a vertex in a SCAN clustering.
using VertexType = std::variant<ClusterMember, Hub, Outlier>;

// Clustering resulting from running SCAN.
struct Clustering {
  // Number of clusters.
  uintE num_clusters;
  // clusters_by_vertex[i] is the SCAN clustering result for vertex i. Cluster
  // IDs are in the range [0, num_clusters).
  pbbs::sequence<VertexType> clusters_by_vertex;
};
bool operator==(const Clustering&, const Clustering&);

// Index for an undirected graph from which clustering the graph with SCAN is
// quick.
//
// Based off of SCAN index presented in "Efficient Structural Graph Clustering:
// An Index-Based Approach" by Wen et al.
class ScanIndex {
 public:
  explicit ScanIndex(symmetric_graph<symmetric_vertex, pbbslib::empty>* graph);

  // Compute a SCAN clustering of the indexed graph using SCAN parameters
  // mu and epsilon.
  //
  // Returns an <number of vertices in graph>-length sequence S in which S[i] is
  // the clustering status of vertex i.
  //
  // Explanation of parameters:
  // - epsilon: a threshold value on the "similarity" between adjacent vertices
  //   based on how much they share neighbors. Increasing this makes
  //   finer-grained, smaller clusters.
  // - mu: how many neighbors a vertex needs to be epsilon-similar to in order
  //   to be considered a "core" vertex from which a cluster is grown.
  //   Increasing this increases the minimum cluster size and but also makes
  //   large clusters less likely to appear.
  Clustering Cluster(uint64_t mu, float epsilon) const;

 private:
  const size_t num_vertices_;
  const internal::NeighborOrder neighbor_order_;
  const internal::CoreOrder core_order_;
};

namespace internal {

using StructuralSimilarities =
  sparse_table<UndirectedEdge, float, std::hash<UndirectedEdge>>;

bool operator==(const NeighborSimilarity&, const NeighborSimilarity&);
std::ostream& operator<<(std::ostream& os, const NeighborSimilarity&);

bool operator==(const CoreThreshold&, const CoreThreshold&);
std::ostream& operator<<(std::ostream& os, const CoreThreshold&);

StructuralSimilarities ComputeStructuralSimilarities(
    symmetric_graph<symmetric_vertex, pbbslib::empty>* graph);

NeighborOrder ComputeNeighborOrder(
    symmetric_graph<symmetric_vertex, pbbslib::empty>* graph,
    const StructuralSimilarities& similarities);

pbbs::sequence<pbbs::sequence<CoreThreshold>> ComputeCoreOrder(
    const NeighborOrder& neighbor_order);

}  // namespace internal

}  // namespace scan
