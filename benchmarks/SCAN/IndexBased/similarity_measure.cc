#include "benchmarks/SCAN/IndexBased/similarity_measure.h"

#include <cmath>
#include <type_traits>

namespace gbbs {

namespace scan {

namespace internal {

sequence<float> RandomNormalNumbers(const size_t num_numbers,
                                    const parlay::random rng) {
  using RngInt = std::result_of<decltype (&parlay::random::ith_rand)(
      parlay::random*, uint64_t)>::type;
  constexpr float max_rng_val{
      static_cast<float>(std::numeric_limits<RngInt>::max())};

  sequence<float> normals{sequence<float>::uninitialized(num_numbers)};
  // Generate normal numbers with the Box–Muller transform.
  parallel_for(0, (num_numbers + 1) / 2, [&](const size_t i) {
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

}  // namespace internal

bool operator==(const EdgeSimilarity& a, const EdgeSimilarity& b) {
  return std::tie(a.source, a.neighbor, a.similarity) ==
         std::tie(b.source, b.neighbor, b.similarity);
}

std::ostream& operator<<(std::ostream& os,
                         const EdgeSimilarity& edge_similarity) {
  os << "{edge=(" << edge_similarity.source << ',' << edge_similarity.neighbor
     << "), similarity=" << edge_similarity.similarity << '}';
  return os;
}

ApproxCosineSimilarity::ApproxCosineSimilarity(const uint32_t num_samples,
                                               const size_t random_seed)
    : num_samples_{num_samples}, random_seed_{random_seed} {}

ApproxJaccardSimilarity::ApproxJaccardSimilarity(const uint32_t num_samples,
                                                 const size_t random_seed)
    : num_samples_{num_samples}, random_seed_{random_seed} {}

}  // namespace scan

}  // namespace gbbs
