#pragma once

#include "benchmarks/SCAN/IndexBased/scan_helpers.h"
#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "pbbslib/seq.h"
#include "pbbslib/utilities.h"
#include "benchmarks/SCAN/IndexBased/utils.h"

namespace indexed_scan {

using scan::Clustering;
using scan::kUnclustered;

// TODO add comment
class CosineSimilaritiesFunctor {
 public:
  CosineSimilaritiesFunctor() = default;

  template <template <typename> class VertexTemplate>
  pbbs::sequence<internal::EdgeSimilarity>
  operator()(symmetric_graph<VertexTemplate, pbbslib::empty>* graph) const {
    return internal::CosineSimilaritiesImpl(graph);
  }
};

// TODO add comment
class ApproxCosineSimilaritiesFunctor {
 public:
  ApproxCosineSimilaritiesFunctor(size_t num_samples, size_t random_seed)
    : num_samples_(num_samples), random_seed_(random_seed) {}

  template <template <typename> class VertexTemplate>
  pbbs::sequence<internal::EdgeSimilarity>
  operator()(symmetric_graph<VertexTemplate, pbbslib::empty>* graph) const {
    return
      internal::ApproxCosineSimilaritiesImpl(graph, num_samples_, random_seed_);
  }

 private:
  const size_t num_samples_;
  const size_t random_seed_;
};

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
  //   similarities_func
  //     Determines what similarity function is used when determining the
  //     similarity between pairs of adjacent vertices. Can be
  //     `CosineSimilaritiesFunctor` or `ApproximateCosineSimilaritiesFunctor`.
  //     (More generally, this needs to be a functor taking a graph and
  //     outputting a `pbbs::sequence<internal::EdgeSimilarity>` of every
  //     directed edge in the graph along with a similarity score on each edge.)
  template <
    template <typename> class VertexTemplate,
    class SimilaritiesFunc = CosineSimilaritiesFunctor>
  explicit Index(
      symmetric_graph<VertexTemplate, pbbslib::empty>* graph,
      SimilaritiesFunc&& similarities_func = CosineSimilaritiesFunctor{})
    : num_vertices_{graph->n}
    , neighbor_order_{graph, std::forward<SimilaritiesFunc>(similarities_func)}
    , core_order_{neighbor_order_} {}

  // Compute a SCAN clustering of the indexed graph using SCAN parameters
  // mu and epsilon.
  //
  // Those who are familiar with SCAN may know that some "border" vertices of
  // clusters can belong to multiple clusters at once. This implementation
  // non-deterministically picks an arbitrary choice of a single cluster
  // assignment for those vertices.
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
  //
  // Returns:
  //   `graph->n`-length sequence S where S[i] is the cluster ID of vertex i or
  //   is `kUnclustered` if vertex i does not belong to any cluster. The cluster
  //   IDs will be in the range `[0, graph->n)` but will not necessarily be
  //   contiguous.
  Clustering Cluster(uint64_t mu, float epsilon) const;

 private:
  const size_t num_vertices_;
  const internal::NeighborOrder neighbor_order_;
  const internal::CoreOrder core_order_;
};

}  // namespace indexed_scan
