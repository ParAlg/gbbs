#pragma once

#include "benchmarks/SCAN/IndexBased/scan_helpers.h"
#include "benchmarks/SCAN/IndexBased/similarity_measure.h"
#include "benchmarks/SCAN/IndexBased/utils.h"
#include "gbbs/graph.h"
#include "gbbs/macros.h"

namespace gbbs {
namespace indexed_scan {

using scan::Clustering;
using scan::kUnclustered;

// Index for an undirected graph from which clustering the graph with SCAN is
// quick, though index construction may be expensive.
class Index {
 public:
  // Constructor.
  //
  // Arguments:
  //   graph
  //     The graph on which to construct the index. The neighbor lists for each
  //     vertex in the graph must be sorted by ascending neighbor ID.
  //   similarity_measure: similarity measure from `similarity_measure.h`
  //     Determines how to compute the similarity between two adjacency
  //     vertices. The traditional choice for SCAN is `scan::CosineSimilarity`.
  template <template <typename> class VertexTemplate, typename Weight,
            class SimilarityMeasure = scan::CosineSimilarity>
  explicit Index(
      symmetric_graph<VertexTemplate, Weight>* graph,
      const SimilarityMeasure& similarity_measure = scan::CosineSimilarity{})
      : num_vertices_{graph->n},
        neighbor_order_{graph, similarity_measure},
        core_order_{neighbor_order_} {}

  Index();

  // Compute a SCAN clustering of the indexed graph using SCAN parameters
  // mu and epsilon.
  //
  // Those who are familiar with SCAN may know that some "border" vertices of
  // clusters can belong to multiple clusters at once. This implementation, by
  // default, non-deterministically picks an arbitrary choice of a single
  // cluster assignment for those vertices.
  //
  // Arguments:
  //   epsilon
  //     A threshold value on the similarity between adjacent vertices based
  //     on how much they share neighbors. Increasing this makes finer-grained,
  //     smaller clusters.
  //   mu
  //     How many neighbors a vertex needs to be epsilon-similar to in order to
  //     be considered a "core" vertex from which a cluster is grown.
  //     Increasing this increases the minimum cluster size and but also makes
  //     large clusters less likely to appear.
  //   get_deterministic_result
  //     If this is set to true, the output result will be deterministic, but
  //     but the computation time of this function may be longer.
  //
  // Returns:
  //   `graph->n`-length sequence S where S[i] is the cluster ID of vertex i or
  //   is `kUnclustered` if vertex i does not belong to any cluster. The cluster
  //   IDs will be in the range `[0, graph->n)` but will not necessarily be
  //   contiguous.
  Clustering Cluster(uint64_t mu, float epsilon,
                     bool get_deterministic_result = false) const;

  // Computes clusterings for several epsilon values.
  // Instead of returning all the clusterings, which may be memory intensive,
  // this function runs f on each clustering as follows:
  //   f(<clustering with parameters (mu, epsilons[i])>, i) for each i.
  // f is called sequentially, but the order in which f is called on each i is
  // arbitrary.
  //
  // TODO(tomtseng) change std::function argument to a templated argument
  void Cluster(uint64_t mu, const sequence<float>& epsilons,
               const std::function<void(Clustering&&, size_t)> f,
               bool get_deterministic_result = false) const;

 private:
  size_t num_vertices_;
  internal::NeighborOrder neighbor_order_;
  internal::CoreOrder core_order_;
};

}  // namespace indexed_scan
}  // namespace gbbs
