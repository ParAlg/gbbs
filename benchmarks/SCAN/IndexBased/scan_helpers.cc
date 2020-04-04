#include "benchmarks/SCAN/IndexBased/scan_helpers.h"

#include <limits>

namespace indexed_scan {

namespace {

// Holds a vertex and its degree.
struct VertexDegree {
  uintE vertex_id;
  uintE degree;
};

}  // namespace

namespace internal {

bool operator==(const EdgeSimilarity& a, const EdgeSimilarity& b) {
  constexpr float kEpsilon{1e-6};
  return std::tie(a.source, a.neighbor) == std::tie(b.source, b.neighbor) &&
    std::abs(a.similarity - b.similarity) < kEpsilon;
}

std::ostream&
operator<<(std::ostream& os, const EdgeSimilarity& edge_similarity) {
  os << "{edge=(" << edge_similarity.source << ',' << edge_similarity.neighbor
     << "), similarity=" << edge_similarity.similarity << '}';
  return os;
}

bool operator==(const CoreThreshold& a, const CoreThreshold& b) {
  constexpr float kEpsilon{1e-6};
  return a.vertex_id == b.vertex_id &&
    std::abs(a.threshold - b.threshold) < kEpsilon;
}

std::ostream&
operator<<(std::ostream& os, const CoreThreshold& core_threshold) {
  os << "{vertex=" << core_threshold.vertex_id
     << ", threshold=" << core_threshold.threshold << '}';
  return os;
}

void ReportTime([[maybe_unused]] const timer& t) {
#ifdef SCAN_DETAILED_TIMES
  t.reportTotal("");
#endif
}

const pbbs::range<EdgeSimilarity*>&
NeighborOrder::operator[](size_t source) const {
  return similarities_by_source_[source];
}

bool NeighborOrder::empty() const {
  return size() == 0;
}

size_t NeighborOrder::size() const {
  return similarities_by_source_.size();
}

pbbs::range<EdgeSimilarity*>* NeighborOrder::begin() const {
  return similarities_by_source_.begin();
}

pbbs::range<EdgeSimilarity*>* NeighborOrder::end() const {
  return similarities_by_source_.end();
}

pbbs::sequence<pbbs::sequence<CoreThreshold>> ComputeCoreOrder(
    const NeighborOrder& neighbor_order) {
  if (neighbor_order.empty()) {
    return {};
  }

  timer function_timer{"Compute core order time"};

  pbbs::sequence<VertexDegree> vertex_degrees{
    pbbs::map_with_index<VertexDegree>(
        neighbor_order,
        [](const size_t v,
           const pbbs::range<EdgeSimilarity*>& neighbors) {
          return VertexDegree{
            .vertex_id = static_cast<uintE>(v),
            .degree = static_cast<uintE>(neighbors.size())};
        })
  };
  // Sort `vertex_degrees` by ascending degree.
  integer_sort_inplace(
      vertex_degrees.slice(),
      [](const VertexDegree& vertex_degree) { return vertex_degree.degree; });
  const size_t max_degree{vertex_degrees[vertex_degrees.size() - 1].degree};

  // `degree_offsets[j]` is the first index `i` at which
  // `vertex_degrees[i].degree >= j`.
  pbbs::sequence<uintE> degree_offsets{
    pbbs::sequence<uintE>::no_init(max_degree + 1)};
  const size_t min_degree{vertex_degrees[0].degree};
  par_for(0, min_degree + 1, [&](const size_t j) {
    degree_offsets[j] = 0;
  });
  par_for(1, vertex_degrees.size(), [&](const size_t i) {
    const size_t degree{vertex_degrees[i].degree};
    const size_t prev_degree{vertex_degrees[i - 1].degree};
    if (degree != prev_degree) {
      par_for(prev_degree + 1, degree + 1, [&](const size_t j) {
        degree_offsets[j] = i;
      });
    }
  });

  const auto get_core_order{[&](const size_t mu) {
    if (mu <= 1) {
      return pbbs::sequence<CoreThreshold>{};
    }

    // Only vertices with high enough degree can be cores.
    const pbbs::sequence<VertexDegree>& core_vertices{
      vertex_degrees.slice(degree_offsets[mu - 1], vertex_degrees.size())};

    pbbs::sequence<CoreThreshold> core_thresholds{
      pbbs::map<CoreThreshold>(
        core_vertices,
        [&](const VertexDegree& vertex_degree) {
          return CoreThreshold{
            .vertex_id = vertex_degree.vertex_id,
            .threshold =
                neighbor_order[vertex_degree.vertex_id][mu - 2].similarity};
        })};
    // Sort by descending threshold.
    const auto compare_threshold_descending{
      [](const CoreThreshold& a, const CoreThreshold& b) {
        return a.threshold > b.threshold;
      }};
    pbbs::sample_sort_inplace(
        core_thresholds.slice(), compare_threshold_descending);
    return core_thresholds;
  }};

  pbbs::sequence<pbbs::sequence<CoreThreshold>> core_order(
      max_degree + 2, get_core_order);
  internal::ReportTime(function_timer);
  return core_order;
}

CoreOrder::CoreOrder(const NeighborOrder& neighbor_order)
    : num_vertices_{neighbor_order.size()}
    , order_{ComputeCoreOrder(neighbor_order)} {}

pbbs::sequence<uintE>
CoreOrder::GetCores(const uint64_t mu, const float epsilon) const {
  if (mu <= 1) {  // All vertices are cores.
    return
      pbbs::sequence<uintE>{num_vertices_, [](const size_t i) { return i; }};
  }
  if (mu >= order_.size()) {  // No vertices are cores.
    return {};
  }

  const pbbs::sequence<CoreThreshold>& possible_cores(order_[mu]);
  const size_t cores_end{
    BinarySearch(
        possible_cores,
        [epsilon](const internal::CoreThreshold& core_threshold) {
          return core_threshold.threshold >= epsilon;
        })};
  return pbbs::map<uintE>(
      possible_cores.slice(0, cores_end),
      [](const internal::CoreThreshold& core_threshold) {
        return core_threshold.vertex_id;
      });
}

}  // namespace internal

}  // namespace indexed_scan
