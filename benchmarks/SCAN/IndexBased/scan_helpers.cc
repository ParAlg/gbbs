#define NOTMAIN

#include "benchmarks/SCAN/scan_helpers.h"

namespace scan {

namespace {

// Holds a vertex and its degree.
struct VertexDegree {
  uintE vertex_id;
  uintE degree;
};

}  // namespace

namespace internal {

bool operator==(const NeighborSimilarity& a, const NeighborSimilarity& b) {
  return
    std::tie(a.neighbor, a.similarity) == std::tie(b.neighbor, b.similarity);
}

std::ostream&
operator<<(std::ostream& os, const NeighborSimilarity& neighbor_similarity) {
  os << "{neighbor=" << neighbor_similarity.neighbor
     << ", similarity=" << neighbor_similarity.similarity << '}';
  return os;
}

bool operator==(const CoreThreshold& a, const CoreThreshold& b) {
  return
    std::tie(a.vertex_id, a.threshold) == std::tie(b.vertex_id, b.threshold);
}

std::ostream&
operator<<(std::ostream& os, const CoreThreshold& core_threshold) {
  os << "{vertex=" << core_threshold.vertex_id
     << ", threshold=" << core_threshold.threshold << '}';
  return os;
}

pbbs::sequence<pbbs::sequence<CoreThreshold>> ComputeCoreOrder(
    const NeighborOrder& neighbor_order) {
  if (neighbor_order.empty()) {
    return {};
  }

  pbbs::sequence<VertexDegree> vertex_degrees{
    pbbs::map_with_index<VertexDegree>(
        neighbor_order,
        [](const size_t v,
           const pbbs::sequence<NeighborSimilarity>& neighbors) {
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

  // `bucket_offsets[i]` is the first index `j` at which
  // `vertex_degrees[j].degree >= i`.
  pbbs::sequence<uintE> degree_offsets{
    pbbs::sequence<uintE>::no_init(max_degree + 1)};
  const size_t min_degree{vertex_degrees[0].degree};
  par_for(0, min_degree + 1, [&](const size_t j) {
    degree_offsets[j] = 0;
  });
  par_for(1, vertex_degrees.size(), [&](const size_t i) {
    const size_t degree{vertex_degrees[i].degree};
    const size_t prev_degree{vertex_degrees[i-1].degree};
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

  return pbbs::sequence<pbbs::sequence<CoreThreshold>>{
      max_degree + 2,
      get_core_order};
}

CoreOrder::CoreOrder(const NeighborOrder& neighbor_order)
    : num_vertices_{neighbor_order.size()}
    , order_{ComputeCoreOrder(neighbor_order)} {}

pbbs::sequence<uintE>
CoreOrder::GetCores(const uint64_t mu, const float epsilon) const {
  if (mu <= 1) {  // All vertices are cores.
    return
      pbbs::sequence<uintE>(num_vertices_, [](const size_t i) { return i; });
  }
  if (mu >= order_.size()) {  // No vertices are cores.
    return pbbs::sequence<uintE>{};
  }

  const pbbs::sequence<CoreThreshold>& possible_cores{order_[mu]};
  const size_t cores_end{
    BinarySearch<internal::CoreThreshold>(
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

}  // namespace scan
