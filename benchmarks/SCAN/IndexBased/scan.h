#pragma once

#include "benchmarks/SCAN/IndexBased/scan_helpers.h"
#include "cereal/access.hpp"
#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "pbbslib/seq.h"
#include "pbbslib/utilities.h"
#include "benchmarks/SCAN/IndexBased/utils.h"

namespace indexed_scan {

using scan::Clustering;
using scan::kUnclustered;

// Index for an undirected graph from which clustering the graph with SCAN is
// quick, though index construction may be expensive.
//
// The Index may be serialized and deserialized as a cereal
// (https://uscilab.github.io/cereal/) archive.
class Index {
 public:
  // Constructor.
  //
  // The neighbor lists for each vertex in the graph must be sorted by ascending
  // neighbor ID.
  template <template <typename> class VertexTemplate>
  explicit Index(
      symmetric_graph<VertexTemplate, pbbslib::empty>* graph)
    : num_vertices_{graph->n}
    , neighbor_order_{graph}
    , core_order_{neighbor_order_} {}

  // Empty index constructor.
  Index();

  // Compute a SCAN clustering of the indexed graph using SCAN parameters
  // mu and epsilon.
  //
  // Those who are familiar with SCAN may know that some "border" vertices of
  // clusters can belong to multiple clusters at once. This implementation picks
  // an arbitrary choice of a single cluster assignment for those vertices.
  //
  // Arguments:
  //   epsilon
  //     A threshold value on the "similarity" between adjacent vertices based
  //     on how much they share neighbors. Increasing this makes finer-grained,
  //     smaller clusters.
  //   mu
  //     How many neighbors a vertex needs to be epsilon-similar to in order to
  //     be considered a "core" vertex from which a cluster is grown.
  //     Increasing this increases the minimum cluster size and but also makes
  //     large clusters less likely to appear.
  //
  // Returns:
  //   `graph->n`-length sequence S where S[i] is the cluster ID of vertex i or
  //   is `kUnclustered` if vertex i does not belong to any cluster. The cluster
  //   IDs will be in the range `[0, graph->n)` but will not necessarily be
  //   contiguous.
  Clustering Cluster(uint64_t mu, float epsilon) const;

 private:
  friend class cereal::access;
  friend bool operator==(const Index&, const Index&);

  // For cereal serialization.
  template<class CerealArchive>
  void serialize(CerealArchive& archive) {
    archive(num_vertices_, neighbor_order_, core_order_);
  }

  size_t num_vertices_;
  internal::NeighborOrder neighbor_order_;
  internal::CoreOrder core_order_;
};

}  // namespace indexed_scan
