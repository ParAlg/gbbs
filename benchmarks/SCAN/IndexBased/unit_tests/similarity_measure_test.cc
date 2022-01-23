#include "benchmarks/SCAN/IndexBased/similarity_measure.h"

#include <limits>
#include <unordered_set>

#include "benchmarks/SCAN/IndexBased/unit_tests/similarity_measure_test_utils.h"
#include "gbbs/macros.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace gbbs {

namespace gt = graph_test;
namespace s = scan;
namespace si = scan::internal;

using ::testing::AllOf;
using ::testing::Field;
using ::testing::Eq;
using ::testing::FloatNear;
using ::testing::UnorderedElementsAre;
using ::testing::UnorderedElementsAreArray;

namespace {

// Googletest-like equality matcher for `scan::EdgeSimilarity` where the
// tolerance on the `similarity` field is a parameter.
auto EdgeSimilarityApproxEq(const uintE expected_source,
                            const uintE expected_neighbor,
                            const float expected_similarity,
                            const float similarity_max_absolute_error) {
  return AllOf(
      Field(&scan::EdgeSimilarity::source, Eq(expected_source)),
      Field(&scan::EdgeSimilarity::neighbor, Eq(expected_neighbor)),
      Field(&scan::EdgeSimilarity::similarity,
            FloatNear(expected_similarity, similarity_max_absolute_error)));
}

// Returns a small graph.
// Graph diagram with cosine similarity scores labeled:
//       .71    .67  .63
//     0 --- 1 ---- 2 -- 5
//           |    / |
//       .75 |  .89 | .77
//           | /    |
//           3 ---- 4
//             .87
// Graph diagram with Jaccard similarity scores labeled:
//       .5     .5    .4
//     0 --- 1 ---- 2 -- 5
//           |    / |
//       .6  |  .8  | .6
//           | /    |
//           3 ---- 4
//             .75
auto MakeBasicGraph() {
  constexpr size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {1, 2}, {1, 3}, {2, 3}, {2, 4}, {2, 5}, {3, 4},
  };
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  return graph;
}

}  // namespace

TEST(CosineSimilarity, AllEdges) {
  auto graph{MakeBasicGraph()};
  const s::CosineSimilarity similarity_measure{};
  // Check `CosineSimilarity::AllEdges` output.
  const sequence<s::EdgeSimilarity> similarities{
      similarity_measure.AllEdges(&graph)};
  EXPECT_THAT(similarities,
              UnorderedElementsAre(EdgeSimilarityEq(0, 1, 2.0 / sqrt(8)),
                                   EdgeSimilarityEq(1, 0, 2.0 / sqrt(8)),
                                   EdgeSimilarityEq(1, 2, 3.0 / sqrt(20)),
                                   EdgeSimilarityEq(2, 1, 3.0 / sqrt(20)),
                                   EdgeSimilarityEq(1, 3, 3.0 / sqrt(16)),
                                   EdgeSimilarityEq(3, 1, 3.0 / sqrt(16)),
                                   EdgeSimilarityEq(2, 3, 4.0 / sqrt(20)),
                                   EdgeSimilarityEq(3, 2, 4.0 / sqrt(20)),
                                   EdgeSimilarityEq(2, 4, 3.0 / sqrt(15)),
                                   EdgeSimilarityEq(4, 2, 3.0 / sqrt(15)),
                                   EdgeSimilarityEq(2, 5, 2.0 / sqrt(10)),
                                   EdgeSimilarityEq(5, 2, 2.0 / sqrt(10)),
                                   EdgeSimilarityEq(3, 4, 3.0 / sqrt(12)),
                                   EdgeSimilarityEq(4, 3, 3.0 / sqrt(12))));

  // Also check that `ApproxCosine::AllEdges` with its threshold tuned so that
  // it always outputs exact similarities also gives the same output.
  constexpr uint32_t kNumSamples{10};
  constexpr size_t kRandomSeed{0};
  constexpr size_t kDegreeThreshold{std::numeric_limits<size_t>::max()};
  const sequence<s::EdgeSimilarity> approx_similarities{
      si::ApproxCosineEdgeSimilarities(&graph, kNumSamples, kDegreeThreshold,
                                       kRandomSeed)};
  EXPECT_THAT(approx_similarities, UnorderedElementsAreArray(similarities));
}

TEST(ApproxCosineSimilarity, AllEdges) {
  auto graph{MakeBasicGraph()};
  // This tests `scan::ApproxCosineSimilarity::AllEdges`, which has a
  // pseudorandom output. Test failures might just be unlucky and fixable by
  // changing the random seed. At the time of writing this comment, the test
  // passes on all random seeds in the range [0, 99).
  constexpr uint32_t kNumSamples{400};
  constexpr size_t kRandomSeed{0};
  constexpr size_t kDegreeThreshold{0};  // Approximate all similarities
  const sequence<s::EdgeSimilarity> similarities{
      si::ApproxCosineEdgeSimilarities(&graph, kNumSamples, kDegreeThreshold,
                                       kRandomSeed)};
  constexpr float kTolerance{0.1};
  EXPECT_THAT(similarities,
              UnorderedElementsAre(
                  EdgeSimilarityApproxEq(0, 1, 2.0 / sqrt(8), kTolerance),
                  EdgeSimilarityApproxEq(1, 0, 2.0 / sqrt(8), kTolerance),
                  EdgeSimilarityApproxEq(1, 2, 3.0 / sqrt(20), kTolerance),
                  EdgeSimilarityApproxEq(2, 1, 3.0 / sqrt(20), kTolerance),
                  EdgeSimilarityApproxEq(1, 3, 3.0 / sqrt(16), kTolerance),
                  EdgeSimilarityApproxEq(3, 1, 3.0 / sqrt(16), kTolerance),
                  EdgeSimilarityApproxEq(2, 3, 4.0 / sqrt(20), kTolerance),
                  EdgeSimilarityApproxEq(3, 2, 4.0 / sqrt(20), kTolerance),
                  EdgeSimilarityApproxEq(2, 4, 3.0 / sqrt(15), kTolerance),
                  EdgeSimilarityApproxEq(4, 2, 3.0 / sqrt(15), kTolerance),
                  EdgeSimilarityApproxEq(2, 5, 2.0 / sqrt(10), kTolerance),
                  EdgeSimilarityApproxEq(5, 2, 2.0 / sqrt(10), kTolerance),
                  EdgeSimilarityApproxEq(3, 4, 3.0 / sqrt(12), kTolerance),
                  EdgeSimilarityApproxEq(4, 3, 3.0 / sqrt(12), kTolerance)));
}

TEST(JaccardSimilarity, AllEdges) {
  auto graph{MakeBasicGraph()};
  constexpr s::JaccardSimilarity similarity_measure{};
  const sequence<s::EdgeSimilarity> similarities{
      similarity_measure.AllEdges(&graph)};
  EXPECT_THAT(similarities,
              UnorderedElementsAre(
                  EdgeSimilarityEq(0, 1, 0.5), EdgeSimilarityEq(1, 0, 0.5),
                  EdgeSimilarityEq(1, 2, 0.5), EdgeSimilarityEq(2, 1, 0.5),
                  EdgeSimilarityEq(1, 3, 0.6), EdgeSimilarityEq(3, 1, 0.6),
                  EdgeSimilarityEq(2, 3, 0.8), EdgeSimilarityEq(3, 2, 0.8),
                  EdgeSimilarityEq(2, 4, 0.6), EdgeSimilarityEq(4, 2, 0.6),
                  EdgeSimilarityEq(2, 5, 0.4), EdgeSimilarityEq(5, 2, 0.4),
                  EdgeSimilarityEq(3, 4, 0.75), EdgeSimilarityEq(4, 3, 0.75)));

  // Also check that `ApproxJaccard::AllEdges` with its threshold tuned so that
  // it always outputs exact similarities also gives the same output.
  constexpr uint32_t kNumSamples{10};
  constexpr size_t kRandomSeed{0};
  constexpr size_t kDegreeThreshold{std::numeric_limits<size_t>::max()};
  const sequence<s::EdgeSimilarity> approx_similarities{
      si::ApproxJaccardEdgeSimilarities(&graph, kNumSamples, kDegreeThreshold,
                                        kRandomSeed)};
  EXPECT_THAT(approx_similarities, UnorderedElementsAreArray(similarities));
}

}  // namespace gbbs
