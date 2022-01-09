#include "benchmarks/SCAN/IndexBased/scan_helpers.h"

#include <limits>

namespace gbbs {
namespace indexed_scan {

namespace {

// Holds a vertex and its degree.
struct VertexDegree {
  uintE vertex_id;
  uintE degree;
};

}  // namespace

namespace internal {

std::ostream& operator<<(std::ostream& os,
                         const CoreThreshold& core_threshold) {
  os << "{vertex=" << core_threshold.vertex_id
     << ", threshold=" << core_threshold.threshold << '}';
  return os;
}

void ReportTime([[maybe_unused]] const timer& t) {
#ifdef SCAN_DETAILED_TIMES
  t.next("");
#endif
}

NeighborOrder::NeighborOrder() : similarities_{}, similarities_by_source_{} {}

const gbbs::slice<EdgeSimilarity>& NeighborOrder::operator[](
    size_t source) const {
  return similarities_by_source_[source];
}

bool NeighborOrder::empty() const { return size() == 0; }

size_t NeighborOrder::size() const { return similarities_by_source_.size(); }

const gbbs::slice<EdgeSimilarity>* NeighborOrder::begin() const {
  return similarities_by_source_.begin();
}

const gbbs::slice<EdgeSimilarity>* NeighborOrder::end() const {
  return similarities_by_source_.end();
}

sequence<sequence<CoreThreshold>> ComputeCoreOrder(
    const NeighborOrder& neighbor_order) {
  if (neighbor_order.empty()) {
    return {};
  }

  timer function_timer{"Compute core order time"};

  sequence<VertexDegree> vertex_degrees{parlay::map_with_index<VertexDegree>(
      neighbor_order,
      [](const size_t v, const gbbs::slice<EdgeSimilarity>& neighbors) {
        return VertexDegree{.vertex_id = static_cast<uintE>(v),
                            .degree = static_cast<uintE>(neighbors.size())};
      })};
  // Sort `vertex_degrees` by ascending degree.
  integer_sort_inplace(
      make_slice(vertex_degrees),
      [](const VertexDegree& vertex_degree) { return vertex_degree.degree; });
  const size_t max_degree{vertex_degrees[vertex_degrees.size() - 1].degree};

  // `degree_offsets[j]` is the first index `i` at which
  // `vertex_degrees[i].degree >= j`.
  sequence<uintE> degree_offsets =
      sequence<uintE>::uninitialized(max_degree + 1);
  const size_t min_degree{vertex_degrees[0].degree};
  parallel_for(0, min_degree + 1,
               [&](const size_t j) { degree_offsets[j] = 0; });
  parallel_for(1, vertex_degrees.size(), [&](const size_t i) {
    const size_t degree{vertex_degrees[i].degree};
    const size_t prev_degree{vertex_degrees[i - 1].degree};
    if (degree != prev_degree) {
      parallel_for(prev_degree + 1, degree + 1,
                   [&](const size_t j) { degree_offsets[j] = i; });
    }
  });

  const auto get_core_order{[&](const size_t mu) {
    if (mu <= 1) {
      return sequence<CoreThreshold>{};
    }

    // Only vertices with high enough degree can be cores.
    auto core_vertices =
        vertex_degrees.cut(degree_offsets[mu - 1], vertex_degrees.size());

    sequence<CoreThreshold> core_thresholds{parlay::map<CoreThreshold>(
        core_vertices, [&](const VertexDegree& vertex_degree) {
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
    parlay::sample_sort_inplace(make_slice(core_thresholds),
                                compare_threshold_descending);
    return core_thresholds;
  }};

  auto core_order = sequence<sequence<CoreThreshold>>::from_function(
      max_degree + 2, get_core_order);

  internal::ReportTime(function_timer);
  return core_order;
}

CoreOrder::CoreOrder(const NeighborOrder& neighbor_order)
    : num_vertices_{neighbor_order.size()},
      order_{ComputeCoreOrder(neighbor_order)} {}

sequence<uintE> CoreOrder::GetCores(const uint64_t mu,
                                    const float epsilon) const {
  if (mu <= 1) {  // All vertices are cores.
    return sequence<uintE>::from_function(num_vertices_,
                                          [](const size_t i) { return i; });
  }
  if (mu >= order_.size()) {  // No vertices are cores.
    return {};
  }

  const sequence<CoreThreshold>& possible_cores(order_[mu]);
  const size_t cores_end{parlay::binary_search(
      possible_cores, [epsilon](const internal::CoreThreshold& core_threshold) {
        return core_threshold.threshold >= epsilon;
      })};
  return parlay::map<uintE>(possible_cores.cut(0, cores_end),
                            [](const internal::CoreThreshold& core_threshold) {
                              return core_threshold.vertex_id;
                            });
}

}  // namespace internal

}  // namespace indexed_scan
}  // namespace gbbs
