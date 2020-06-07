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

namespace indexed_scan {

namespace internal {

using EdgeSimilarity = scan::EdgeSimilarity;
using NoWeight = pbbslib::empty;

// An adjacency list for the graph in which each vertex's neighbor list is
// sorted by descending structural similarity.
//
// The structural similarity between two vertices u and v is
//   (size of intersection of closed neighborhoods of u and v) /
//   (geometric mean of size of closed neighborhoods of u and of v)
// where the closed neighborhood of a vertex x consists of all neighbors of x
// along with x itself.
//
// Unlike the presentation in "Efficient Structural Graph Clustering:  An
// Index-Based Approach", the neighbor list for a vertex `v` will not contain
// `v` itself.
class NeighborOrder {
 public:
  // Constructor.
  //
  // The neighbor lists for each vertex in the graph must be sorted by ascending
  // neighbor ID.
  //
  // TODO add comment
  template <template <typename> class VertexTemplate, class SimilarityMeasure>
  NeighborOrder(
      symmetric_graph<VertexTemplate, NoWeight>* graph,
      const SimilarityMeasure& similarity_measure);

  // Get all structural similarity scores from vertex `source` to its neighbors,
  // sorted by descending similarity.
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
// Beware that this equality operator compares the floating-point field
// approximately. This is convenient for unit tests but might not be appropriate
// for other uses (e.g., this operator is not transitive)
bool operator==(const CoreThreshold&, const CoreThreshold&);
std::ostream& operator<<(std::ostream& os, const CoreThreshold&);

class CoreOrder {
 public:
  explicit CoreOrder(const NeighborOrder& neighbor_order);

  // Return all vertices that are cores under SCAN parameters `mu` and
  // `epsilon`.
  pbbs::sequence<uintE> GetCores(uint64_t mu, float epsilon) const;

 private:
  const size_t num_vertices_;
  pbbs::sequence<pbbs::sequence<CoreThreshold>> order_{};
};

// Prints the total time captured by `timer` to stderr if macro
// SCAN_DETAILED_TIMES is defined, otherwise does nothing.
void ReportTime(const timer&);

// Finds the least `i` such that `predicate(sequence[i])` is false. If
// `predicate(sequence[i])` is true for all `i`, then this returns
// `sequence.size()`.
//
// `sequence` and `predicate` must be partitioned such that there is some `i`
// (which will be the return value) for which `predicate(sequence[i])` is true
// for all j < i and for which `predicate(sequence[i])` is false for all j >= i.
//
// `predicate` should take a `SeqElement` and return a boolean.
//
// Running time is O(log [return value]).
template <class Seq, class Func>
size_t BinarySearch(const Seq& sequence, Func&& predicate) {
  // Start off the binary search with an exponential search so that running time
  // while be O(log [return value]) rather than O(log sequence.size()).
  // In practice this probably doesn't matter much....
  size_t hi{0};
  while (hi < sequence.size() && predicate(sequence[hi])) {
    hi = 2 * hi + 1;
  }
  const size_t lo{hi / 2};
  if (hi > sequence.size()) {
    hi = sequence.size();
  }

  return lo +
    pbbs::binary_search(sequence.slice(lo, hi), std::forward<Func>(predicate));
}

template <template <typename> class VertexTemplate, class SimilarityMeasure>
NeighborOrder::NeighborOrder(
    symmetric_graph<VertexTemplate, NoWeight>* graph,
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
      [&](const size_t i) { return graph->get_vertex(i).getOutDegree(); }};
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
