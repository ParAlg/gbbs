#include "benchmarks/SCAN/IndexBased/similarity_measure.h"

#include <unordered_set>

#include "benchmarks/SCAN/IndexBased/unit_tests/similarity_measure_test_utils.h"
#include "gbbs/graph_test_utils.h"
#include "gbbs/macros.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace gt = graph_test;
namespace s = scan;

using ::testing::AllOf;
using ::testing::Field;
using ::testing::Eq;
using ::testing::FloatNear;
using ::testing::UnorderedElementsAre;

namespace {

// Googletest-like equality matcher for `scan::EdgeSimilarity` where the
// tolerance on the `similarity` field is a parameter.
auto EdgeSimilarityApproxEq(
    const uintE expected_source,
    const uintE expected_neighbor,
    const float expected_similarity,
    const float similarity_max_absolute_error) {
  return AllOf(
    Field(&scan::EdgeSimilarity::source, Eq(expected_source)),
    Field(&scan::EdgeSimilarity::neighbor, Eq(expected_neighbor)),
    Field(
      &scan::EdgeSimilarity::similarity,
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
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{
    {0, 1},
    {1, 2},
    {1, 3},
    {2, 3},
    {2, 4},
    {2, 5},
    {3, 4},
  };
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  return graph;
}

}  // namespace

TEST(CosineSimilarity, AllEdges) {
  auto graph{MakeBasicGraph()};
  const s::CosineSimilarity similarity_measure{};
  const pbbs::sequence<s::EdgeSimilarity> similarities{
    similarity_measure.AllEdges(&graph)};
  EXPECT_THAT(similarities,
      UnorderedElementsAre(
        EdgeSimilarityEq(0, 1, 2.0 / sqrt(8)),
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
}

TEST(ApproxCosineSimilarity, AllEdges) {
  auto graph{MakeBasicGraph()};

  // TODO comment about this test being pseudo-random
  const uintE kNumSamples{400};
  const uintE kRandomSeed{0};
  const s::ApproxCosineSimilarity similarity_measure{kNumSamples, kRandomSeed};
  const pbbs::sequence<s::EdgeSimilarity> similarities{
    similarity_measure.AllEdges(&graph)};
  const float kTolerance{0.1};
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
  const s::JaccardSimilarity similarity_measure{};
  const pbbs::sequence<s::EdgeSimilarity> similarities{
    similarity_measure.AllEdges(&graph)};
  EXPECT_THAT(similarities,
      UnorderedElementsAre(
        EdgeSimilarityEq(0, 1, 0.5),
        EdgeSimilarityEq(1, 0, 0.5),
        EdgeSimilarityEq(1, 2, 0.5),
        EdgeSimilarityEq(2, 1, 0.5),
        EdgeSimilarityEq(1, 3, 0.6),
        EdgeSimilarityEq(3, 1, 0.6),
        EdgeSimilarityEq(2, 3, 0.8),
        EdgeSimilarityEq(3, 2, 0.8),
        EdgeSimilarityEq(2, 4, 0.6),
        EdgeSimilarityEq(4, 2, 0.6),
        EdgeSimilarityEq(2, 5, 0.4),
        EdgeSimilarityEq(5, 2, 0.4),
        EdgeSimilarityEq(3, 4, 0.75),
        EdgeSimilarityEq(4, 3, 0.75)));
}

TEST(ApproxJaccardSimilarity, AllEdges) {
  auto graph{MakeBasicGraph()};

  // TODO comment about this test being pseudo-random
  const uintE kNumSamples{300};
  const uintE kRandomSeed{0};
  const s::ApproxJaccardSimilarity similarity_measure{kNumSamples, kRandomSeed};
  const pbbs::sequence<s::EdgeSimilarity> similarities{
    similarity_measure.AllEdges(&graph)};
  const float kTolerance{0.1};
  EXPECT_THAT(similarities,
      UnorderedElementsAre(
        EdgeSimilarityApproxEq(0, 1, 0.5, kTolerance),
        EdgeSimilarityApproxEq(1, 0, 0.5, kTolerance),
        EdgeSimilarityApproxEq(1, 2, 0.5, kTolerance),
        EdgeSimilarityApproxEq(2, 1, 0.5, kTolerance),
        EdgeSimilarityApproxEq(1, 3, 0.6, kTolerance),
        EdgeSimilarityApproxEq(3, 1, 0.6, kTolerance),
        EdgeSimilarityApproxEq(2, 3, 0.8, kTolerance),
        EdgeSimilarityApproxEq(3, 2, 0.8, kTolerance),
        EdgeSimilarityApproxEq(2, 4, 0.6, kTolerance),
        EdgeSimilarityApproxEq(4, 2, 0.6, kTolerance),
        EdgeSimilarityApproxEq(2, 5, 0.4, kTolerance),
        EdgeSimilarityApproxEq(5, 2, 0.4, kTolerance),
        EdgeSimilarityApproxEq(3, 4, 0.75, kTolerance),
        EdgeSimilarityApproxEq(4, 3, 0.75, kTolerance)));
}
