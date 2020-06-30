#include "benchmarks/SCAN/IndexBased/similarity_measure.h"

#include <cmath>
#include <type_traits>

namespace gbbs {

namespace scan {

namespace internal {

// Pseudorandomly generate `num_numbers` random normal numbers, each with zero
// mean and unit variance.
pbbs::sequence<float>
RandomNormalNumbers(const size_t num_numbers, const pbbs::random rng) {
  using RngInt =
    std::result_of<
      decltype(&pbbs::random::ith_rand)(pbbs::random*, uint64_t)>::type;
  constexpr float max_rng_val{
    static_cast<float>(std::numeric_limits<RngInt>::max())};

  pbbs::sequence<float> normals{pbbs::sequence<float>::no_init(num_numbers)};
  // Generate normal numbers with the Box–Muller transform.
  par_for(0, (num_numbers + 1) / 2, [&](const size_t i) {
    const float uniform_1{rng.ith_rand(2 * i) / max_rng_val};
    const float uniform_2{rng.ith_rand(2 * i + 1) / max_rng_val};
    const float radius{std::sqrt(-2.0f * std::log(uniform_1))};
    const float angle{static_cast<float>(2.0 * M_PI * uniform_2)};
    normals[2 * i] = radius * std::cos(angle);
    if (2 * i + 1 < num_numbers) {
      normals[2 * i + 1] = radius * std::sin(angle);
    }
  });
  return normals;
}

// Given a sequence of similarities for each undirected edge, construct a
// sequence of similarities for each directed edge, assuming that
// similarity(u, v) = similarity(v, u) = similarity({u, v}) for every pair of
// adjacent vertices u and v.
//
// Arguments:
//   num_directed_edges
//     The number of directed edges in the graph. (Should be `graph->m` for a
//     graph `graph`.)
//   undirectional_similarities
//     A sequence of similarities for each undirected edge {u, v}. Entries with
//     NaN as their `.similarity` field are skipped.
// Given `undirectional_similarities` filled with edge similarities
pbbs::sequence<EdgeSimilarity> BidirectionalSimilarities(
    const size_t num_directed_edges,
    const pbbs::sequence<EdgeSimilarity>& unidirectional_similarities) {
  // Copy similarities for edges (u, v) where u > v.
  const size_t half_num_edges{num_directed_edges / 2};
  pbbs::sequence<EdgeSimilarity> similarities{
    pbbs::sequence<EdgeSimilarity>::no_init(num_directed_edges)};
  constexpr auto is_valid_similarity{
    [](const EdgeSimilarity& edge) { return !std::isnan(edge.similarity); }};
  pbbs::filter_out(
      unidirectional_similarities, similarities.slice(), is_valid_similarity);
  par_for(0, half_num_edges, [&](const size_t i) {
      const EdgeSimilarity& edge{similarities[i]};
      similarities[i + half_num_edges] = {
        .source = edge.neighbor,
        .neighbor = edge.source,
        .similarity = edge.similarity};
  });
  return similarities;
}

}  // namespace internal

bool operator==(const EdgeSimilarity& a, const EdgeSimilarity& b) {
  return std::tie(a.source, a.neighbor, a.similarity) ==
    std::tie(b.source, b.neighbor, b.similarity);
}

std::ostream&
operator<<(std::ostream& os, const EdgeSimilarity& edge_similarity) {
  os << "{edge=(" << edge_similarity.source << ',' << edge_similarity.neighbor
     << "), similarity=" << edge_similarity.similarity << '}';
  return os;
}

ApproxCosineSimilarity::ApproxCosineSimilarity(
    const uint32_t num_samples,
    const size_t random_seed)
  : num_samples_{num_samples}
  , random_seed_{random_seed} {}

ApproxJaccardSimilarity::ApproxJaccardSimilarity(
    const uint32_t num_samples,
    const size_t random_seed)
  : num_samples_{num_samples}
  , random_seed_{random_seed} {}

}  // namespace scan

}  // namespace gbbs
