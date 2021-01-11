// Helper functions for SCAN logic.
//
// Mainly contains template functions that would otherwise clutter about the
// main SCAN header file.
#pragma once

#include <utility>

#include "benchmarks/SCAN/IndexBased/similarity_measure.h"
#include "gbbs/graph.h"
#include "pbbslib/get_time.h"
#include "pbbslib/sample_sort.h"
#include "pbbslib/seq.h"

namespace gbbs {
namespace indexed_scan {

namespace internal {

using EdgeSimilarity = scan::EdgeSimilarity;

// An adjacency list for the graph in which each vertex's neighbor list is
// sorted by descending similarity.
class NeighborOrder {
 public:
  // Constructor where the similarity between two adjacent vertices is
  // determined using a particular similarity measure (see
  // `similarity_measure.h` for options).
  //
  // The neighbor lists for each vertex in the graph must be sorted by ascending
  // neighbor ID.
  template <
    template <typename> class VertexTemplate,
    typename Weight,
    class SimilarityMeasure>
  NeighborOrder(
      symmetric_graph<VertexTemplate, Weight>* graph,
      const SimilarityMeasure& similarity_measure);

  NeighborOrder();

  // Get all similarity scores from vertex `source` to its neighbors (not
  // including `source` itself), sorted by descending similarity.
  const pbbs::range<EdgeSimilarity*>& operator[](size_t source) const;

  bool empty() const;
  // Returns the number of vertices.
  size_t size() const;

  pbbs::range<EdgeSimilarity*>* begin() const;
  pbbs::range<EdgeSimilarity*>* end() const;

 private:
  // Holds similarity scores for all edges, sorted by source and then by
  // similarity.
  pbbs::sequence<EdgeSimilarity> similarities_;
  pbbs::sequence<pbbs::range<EdgeSimilarity*>> similarities_by_source_;
};

struct CoreThreshold {
  uintE vertex_id;
  // Maximum value of the SCAN parameter epsilon for which `vertex_id` is a core
  // vertex (given some fixed reference value for SCAN parameter mu).
  float threshold;
};
std::ostream& operator<<(std::ostream& os, const CoreThreshold&);

class CoreOrder {
 public:
  explicit CoreOrder(const NeighborOrder& neighbor_order);

  // Return all vertices that are cores under SCAN parameters `mu` and
  // `epsilon`.
  pbbs::sequence<uintE> GetCores(uint64_t mu, float epsilon) const;

 private:
  size_t num_vertices_;
  pbbs::sequence<pbbs::sequence<CoreThreshold>> order_{};
};

// Prints the total time captured by `timer` to stderr if macro
// SCAN_DETAILED_TIMES is defined, otherwise does nothing.
void ReportTime(const timer&);

template <
  template <typename> class VertexTemplate,
  typename Weight,
  class SimilarityMeasure>
NeighborOrder::NeighborOrder(
    symmetric_graph<VertexTemplate, Weight>* graph,
    const SimilarityMeasure& similarity_measure) {
  timer function_timer{"Construct neighbor order"};
  similarities_ = similarity_measure.AllEdges(graph);
  pbbs::sample_sort_inplace(
      similarities_.slice(),
      [](const EdgeSimilarity& left, const EdgeSimilarity& right) {
        // Sort by ascending source, then descending similarity.
        return std::tie(left.source, right.similarity) <
          std::tie(right.source, left.similarity);
      });
  pbbs::sequence<uintT> vertex_offsets{
      graph->n,
      [&](const size_t i) { return graph->get_vertex(i).out_degree(); }};
  pbbslib::scan_add_inplace(vertex_offsets);
  similarities_by_source_ = pbbs::sequence<pbbs::range<EdgeSimilarity*>>(
      graph->n,
      [&](const size_t i) {
        return similarities_.slice(
          vertex_offsets[i],
          i + 1 == graph->n ? similarities_.size() : vertex_offsets[i + 1]);
      });
  internal::ReportTime(function_timer);
}

// Returns a sequence CO where CO[i] for i >= 2 is a list of vertices that
// can be a core when the SCAN parameter mu is set to i. The vertices in
// CO[i] are sorted by their core threshold values, the maximum value of SCAN
// parameter epsilon such that the vertex is a core when mu == i.
//
// CO[0] and CO[1] are left empty --- when mu is less than 2, all vertices are
// always cores and have a core threshold of 1.
pbbs::sequence<pbbs::sequence<CoreThreshold>> ComputeCoreOrder(
    const NeighborOrder& neighbor_order);

}  // namespace internal

}  // namespace indexed_scan
}  // namespace gbbs
