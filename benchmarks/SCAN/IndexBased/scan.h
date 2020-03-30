#pragma once

#include <variant>

#include "benchmarks/SCAN/IndexBased/scan_helpers.h"
#include "ligra/graph.h"
#include "pbbslib/seq.h"

namespace indexed_scan {

// Represents a vertex that is in at least one SCAN cluster.
struct ClusterMember {
  // IDs of the clusters that the vertex is in. Vertices at the peripheries of
  // clusters may belong to multiple clusters.
  pbbs::sequence<uintE> clusters;
};
bool operator==(const ClusterMember&, const ClusterMember&);
std::ostream& operator<<(std::ostream&, const ClusterMember&);

// Represents a vertex that is not in a cluster but is adjacent to at least two
// clusters.
struct Hub{};
bool operator==(const Hub&, const Hub&);
std::ostream& operator<<(std::ostream&, const Hub&);

// Represents a vertex that is not in a cluster and is adjacent to at most one
// cluster.
struct Outlier{};
bool operator==(const Outlier&, const Outlier&);
std::ostream& operator<<(std::ostream&, const Outlier&);

// Possible results for a vertex in a SCAN clustering.
using VertexType = std::variant<ClusterMember, Hub, Outlier>;
std::ostream& operator<<(std::ostream&, const VertexType&);

// Clustering resulting from running SCAN.
struct Clustering {
  // Number of clusters.
  uintE num_clusters;
  // clusters_by_vertex[i] is the SCAN clustering result for vertex i. Cluster
  // IDs are in the range [0, num_clusters).
  pbbs::sequence<VertexType> clusters_by_vertex;
};
bool operator==(const Clustering&, const Clustering&);
std::ostream& operator<<(std::ostream&, const Clustering&);

// Index for an undirected graph from which clustering the graph with SCAN is
// quick.
//
// Based off of SCAN index presented in "Efficient Structural Graph Clustering:
// An Index-Based Approach" by Wen et al.
class Index {
 public:
  template <template <typename> class VertexTemplate>
  explicit Index(
      symmetric_graph<VertexTemplate, pbbslib::empty>* graph)
    : num_vertices_{graph->n}
    , neighbor_order_{
        internal::ComputeNeighborOrder(
            graph,
            internal::ComputeStructuralSimilarities(graph))}
    , core_order_{neighbor_order_} {}

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

}  // namespace indexed_scan
