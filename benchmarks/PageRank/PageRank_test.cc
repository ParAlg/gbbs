// TODO: Add tests for `PageRank_edgeMap` with `has_in_edges = false`.

#include "benchmarks/PageRank/PageRank.h"

#include <cmath>
#include <cstddef>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

#include "absl/status/statusor.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "benchmarks/PageRank/PageRank_edgeMapReduce.h"
#include "gbbs/bridge.h"
#include "gbbs/graph.h"
#include "gbbs/helpers/progress_reporting.h"
#include "gbbs/helpers/progress_reporting_mock.h"
#include "gbbs/helpers/status_macros.h"
#include "gbbs/helpers/undirected_edge.h"
#include "gbbs/macros.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gbbs/vertex.h"

using ::testing::DoubleNear;
using ::testing::ElementsAre;
using ::testing::ElementsAreArray;
using ::testing::IsEmpty;
using ::testing::Pointwise;

namespace gbbs {

struct PageRank_edgeMapFunctor {
  template <class Graph>
  static absl::StatusOr<sequence<double>> compute_pagerank(
      Graph& G, double eps = 0.000001, std::vector<uintE> sources = {},
      double damping_factor = 0.85, size_t max_iters = 100,
      std::optional<ReportProgressCallback> report_progress = std::nullopt) {
    return PageRank_edgeMap(G, eps, sources, damping_factor, max_iters,
                            /*has_in_edges=*/true, std::move(report_progress));
  }
};

struct PageRank_edgeMapReduceFunctor {
  template <class Graph>
  static absl::StatusOr<sequence<double>> compute_pagerank(
      Graph& G, double eps = 0.000001, std::vector<uintE> sources = {},
      double damping_factor = 0.85, size_t max_iters = 100,
      std::optional<ReportProgressCallback> report_progress = std::nullopt) {
    return PageRank_edgeMapReduce(G, eps, sources, damping_factor, max_iters,
                                  std::move(report_progress));
  }
};

template <typename T>
class PageRankFixture : public testing::Test, public ReportProgressMock {
 public:
  using Impl = T;
};

using Implementations =
    ::testing::Types<PageRank_edgeMapFunctor, PageRank_edgeMapReduceFunctor>;
TYPED_TEST_SUITE(PageRankFixture, Implementations);

TYPED_TEST(PageRankFixture, EmptyGraph) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_THAT(Impl::compute_pagerank(graph), IsOkAndHolds(IsEmpty()));
}

TYPED_TEST(PageRankFixture, EdgelessGraph) {
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_THAT(Impl::compute_pagerank(graph),
              IsOkAndHolds(ElementsAre(1.0 / 3, 1.0 / 3, 1.0 / 3)));
}

TYPED_TEST(PageRankFixture, ZeroIterations) {
  // Graph diagram:
  //     0 - 1 - 2
  //     3 (isolated)
  //
  constexpr uintE kNumVertices{4};
  const std::unordered_set<UndirectedEdge> kEdges{{0, 1}, {1, 2}};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_THAT(Impl::compute_pagerank(graph, /*eps=*/0.1, /*sources=*/{},
                                     /*damping_factor=*/0.85, /*max_iters=*/0,
                                     this->MockReportProgressCallback()),
              IsOkAndHolds(ElementsAre(0.25, 0.25, 0.25, 0.25)));
  EXPECT_THAT(this->GetProgressReports(), ElementsAre(1.0));
}

TYPED_TEST(PageRankFixture, Cycle) {
  // Graph diagram:
  //     0 - 1 - 2 - 3 - 4 (loops back to 0)
  //
  constexpr uintE kNumVertices{5};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  using Impl = typename TestFixture::Impl;

  {
    EXPECT_THAT(Impl::compute_pagerank(graph),
                IsOkAndHolds(ElementsAre(0.2, 0.2, 0.2, 0.2, 0.2)));
  }
}

TYPED_TEST(PageRankFixture, Path) {
  // Graph diagram:
  //     0 - 1 - 2
  //
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1},
      {1, 2},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  using Impl = typename TestFixture::Impl;

  EXPECT_THAT(
      Impl::compute_pagerank(graph),
      IsOkAndHolds(Pointwise(DoubleNear(1e-4),
                             {0.2567570878, 0.4864858243, 0.2567570878})));
}

TYPED_TEST(PageRankFixture, RespectsEps) {
  // Graph diagram:
  //     0 - 1 - 2
  //
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1},
      {1, 2},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  using Impl = typename TestFixture::Impl;

  {
    // Between iterations 4 and 5 and between iterations 5 and 6, the L1
    // distance between the PageRank values are approximately 0.295804 and
    // 0.251433 respectively. The next lines verify that assumption.
    ASSERT_OK_AND_ASSIGN(
        const sequence<double> result_with_4_iterations,
        Impl::compute_pagerank(graph, /*eps=*/1e-6,
                               /*sources=*/{},
                               /*damping_factor=*/0.85, /*max_iters=*/4,
                               this->MockReportProgressCallback()));
    EXPECT_THAT(this->GetProgressReports(),
                ElementsAre(0.2, 0.4, 0.6, 0.8, 1.0));
    ASSERT_OK_AND_ASSIGN(
        const sequence<double> result_with_5_iterations,
        Impl::compute_pagerank(graph, /*eps=*/1e-6,
                               /*sources=*/{},
                               /*damping_factor=*/0.85, /*max_iters=*/5,
                               this->MockReportProgressCallback()));
    EXPECT_THAT(
        this->GetProgressReports(),
        ElementsAre(1.0 / 6.0, 2.0 / 6.0, 0.5, 4.0 / 6.0, 5.0 / 6.0, 1.0));
    ASSERT_OK_AND_ASSIGN(
        const sequence<double> result_with_6_iterations,
        Impl::compute_pagerank(graph, /*eps=*/1e-6,
                               /*sources=*/{},
                               /*damping_factor=*/0.85, /*max_iters=*/6,
                               this->MockReportProgressCallback()));
    EXPECT_THAT(this->GetProgressReports(),
                ElementsAre(1.0 / 7.0, 2.0 / 7.0, 3.0 / 7.0, 4.0 / 7.0,
                            5.0 / 7.0, 6.0 / 7.0, 1.0));

    double l1_distance_iteration_4_to_5 = 0.0;
    for (int i = 0; i < 3; ++i) {
      l1_distance_iteration_4_to_5 +=
          fabs(result_with_4_iterations[i] - result_with_5_iterations[i]);
    }
    ASSERT_GT(l1_distance_iteration_4_to_5, 0.28);
    ASSERT_LT(l1_distance_iteration_4_to_5, 0.3);

    double l1_distance_iteration_5_to_6 = 0.0;
    for (int i = 0; i < 3; ++i) {
      l1_distance_iteration_5_to_6 +=
          fabs(result_with_5_iterations[i] - result_with_6_iterations[i]);
    }
    ASSERT_GT(l1_distance_iteration_5_to_6, 0.24);
    ASSERT_LT(l1_distance_iteration_5_to_6, 0.26);

    // Run the algorithm with `eps` small enough to reach iteration 5, but not
    // iteration 6.
    EXPECT_THAT(
        Impl::compute_pagerank(graph,
                               /*eps=*/0.31 / 3, /*sources=*/{},
                               /*damping_factor=*/0.85,
                               /*max_iters=*/9, this->MockReportProgressCallback()),
        IsOkAndHolds(ElementsAreArray(result_with_5_iterations)));
    EXPECT_THAT(this->GetProgressReports(),
                ElementsAre(0.1, 0.2, 0.3, 0.4, 0.5, 1.0));

    // Run the algorithm with `eps` small enough to reach iteration
    // 6, but not iteration 7.
    EXPECT_THAT(
        Impl::compute_pagerank(graph,
                               /*eps=*/0.27 / 3, /*sources=*/{},
                               /*damping_factor=*/0.85,
                               /*max_iters=*/9, this->MockReportProgressCallback()),
        IsOkAndHolds(ElementsAreArray(result_with_6_iterations)));
    EXPECT_THAT(this->GetProgressReports(),
                ElementsAre(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0));
  }
}

TYPED_TEST(PageRankFixture, BasicUndirected) {
  // Graph diagram:
  //     0 - 1    2 - 3 - 4
  //                    \ |
  //                      5 -- 6
  constexpr uintE kNumVertices{7};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {2, 3}, {3, 4}, {3, 5}, {4, 5}, {5, 6},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  using Impl = typename TestFixture::Impl;

  const sequence<double> expected{0.1428571429, 0.1428571429, 0.0802049539,
                                  0.2074472351, 0.1389813363, 0.2074472351,
                                  0.0802049539};

  EXPECT_THAT(Impl::compute_pagerank(graph),
              IsOkAndHolds(Pointwise(DoubleNear(1e-4), expected)));
}

TYPED_TEST(PageRankFixture, TwoSources) {
  // Graph diagram:
  //     0 - 1    2 - 3 - 4
  //                    \ |
  //                      5 -- 6
  constexpr uintE kNumVertices{7};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {2, 3}, {3, 4}, {3, 5}, {4, 5}, {5, 6},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  using Impl = typename TestFixture::Impl;

  // The results below were computed by using NetworkX's source-based PageRank
  // as a reference implementation.
  const sequence<double> expected{0.2702716450, 0.2297283550, 0.0384308511,
                                  0.1356382979, 0.1518617021, 0.1356382979,
                                  0.0384308511};
  EXPECT_THAT(Impl::compute_pagerank(graph,
                                     /*eps=*/0.000001, {0, 4}),
              IsOkAndHolds(Pointwise(DoubleNear(1e-4), expected)));
}

TYPED_TEST(PageRankFixture, AllButOneSource) {
  // Graph diagram:
  //     0 - 1 - 2 - 3 - 4 (loops back to 0)
  //
  constexpr uintE kNumVertices{5};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  using Impl = typename TestFixture::Impl;

  const sequence<double> expected{0.2137115266, 0.2137115266, 0.2009045292,
                                  0.1707678884, 0.2009045292};
  EXPECT_THAT(Impl::compute_pagerank(graph,
                                     /*eps=*/0.000001, {0, 1, 2, 4}),
              IsOkAndHolds(Pointwise(DoubleNear(1e-4), expected)));
}

TYPED_TEST(PageRankFixture, TwoNodesNoEdges) {
  constexpr uintE kNumVertices{2};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  using Impl = typename TestFixture::Impl;

  {
    // No mass is migrated.
    const sequence<double> expected{0.5, 0.5};
    EXPECT_THAT(Impl::compute_pagerank(graph,
                                       /*eps=*/0.000001, {0, 1}),
                IsOkAndHolds(Pointwise(DoubleNear(1e-4), expected)));
  }

  {
    // No mass is migrated.
    const sequence<double> expected{1.0, 0};
    EXPECT_THAT(Impl::compute_pagerank(graph,
                                       /*eps=*/0.000001, {0}),
                IsOkAndHolds(Pointwise(DoubleNear(1e-4), expected)));
  }
}

TYPED_TEST(PageRankFixture, TwoNodesOneEdge) {
  // Graph diagram:
  //     0 - 1
  constexpr uintE kNumVertices{2};
  const std::unordered_set<UndirectedEdge> kEdges{{0, 1}};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  using Impl = typename TestFixture::Impl;

  {
    // Both have equal mass; again no mass gets migrated.
    const sequence<double> expected{0.5, 0.5};
    EXPECT_THAT(Impl::compute_pagerank(graph,
                                       /*eps=*/0.000001, {0, 1}),
                IsOkAndHolds(Pointwise(DoubleNear(1e-4), expected)));
  }

  {
    // The steady state for this setting is when P[0] = 20/37 and P[1] = 17/37.
    // This can be verified by plugging in the values above into the PPR
    // equations
    // P[0] = 0.85*P[1] + 0.15
    // P[1] = 0.85*P[0]
    const sequence<double> expected{0.5405409317, 0.4594590684};
    EXPECT_THAT(Impl::compute_pagerank(graph,
                                       /*eps=*/0.000001, {0}),
                IsOkAndHolds(Pointwise(DoubleNear(1e-4), expected)));
  }
}

TYPED_TEST(PageRankFixture, FourNodePathOneSource) {
  // Graph diagram:
  //     0 - 1 - 2 - 3
  constexpr uintE kNumVertices{4};
  const std::unordered_set<UndirectedEdge> kEdges{{0, 1}, {1, 2}, {2, 3}};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  using Impl = typename TestFixture::Impl;

  {
    // The values below can be verified as being stationary for the following
    // equations (note that with only one source, the addedConstant is only
    // present for node 0)
    // P[0] = 0.85*P[1]/2 + 0.15/4
    // P[1] = 0.85*(P[0] + P[2]/2) + 0.15/4
    // P[2] = 0.85*(P[1]/2 + P[3]) + 0.15/4
    // P[3] = 0.85*P[2]/2 + 0.15/4
    const sequence<double> expected{0.17543856058862931, 0.32456143941137061,
                                    0.32456143941137061, 0.17543856058862931};
    EXPECT_THAT(Impl::compute_pagerank(graph,
                                       /*eps=*/0.000001, {0, 1, 2, 3}),
                IsOkAndHolds(Pointwise(DoubleNear(1e-4), expected)));
  }

  {
    // The values below can be verified as being stationary for the following
    // equations (note that with only one source, the addedConstant is only
    // present for node 0)
    // P[0] = 0.85*P[1]/2 + 0.15
    // P[1] = 0.85*(P[0] + P[2]/2)
    // P[2] = 0.85*(P[1]/2 + P[3])
    // P[3] = 0.85*P[2]/2
    const sequence<double> expected{0.3022241274, 0.3581756964, 0.2383155317,
                                    0.1012846445};
    EXPECT_THAT(Impl::compute_pagerank(graph,
                                       /*eps=*/0.000001, {0}),
                IsOkAndHolds(Pointwise(DoubleNear(1e-4), expected)));
  }

  {
    // The values below can be verified as being stationary for the following
    // equations (now that we have two sources, the added constant is split
    // between the two sources).
    // P[0] = 0.85*P[1]/2 + 0.15/2
    // P[1] = 0.85*(P[0] + P[2]/2)
    // P[2] = 0.85*(P[1]/2 + P[3])
    // P[3] = 0.85*P[2]/2 + 0.15/2
    const sequence<double> expected{0.2017542424, 0.2982457576, 0.2982457576,
                                    0.2017542424};
    EXPECT_THAT(Impl::compute_pagerank(graph,
                                       /*eps=*/0.000001, {0, 3}),
                IsOkAndHolds(Pointwise(DoubleNear(1e-4), expected)));
  }
}

class PageRankEdgeMapFunctorTest : public testing::Test,
                                   public ReportProgressMock {};

TEST_F(PageRankEdgeMapFunctorTest, WeightedDigraph) {
  // Graph diagram (each number wrapped in parentheses denotes the weight of its
  // nearest edge):
  //
  //            3
  //           ^ |
  //       (3) | | (4)
  //           | v
  //   0 ---->  1  <---- 2
  //      (1)   |    (2)
  //            | (5)
  //            v
  //            4
  //

  using WeightedDigraph = asymmetric_graph<asymmetric_vertex, float>;
  sequence<WeightedDigraph::edge> edges(
      {{0, 1, 1.0}, {1, 3, 3.0}, {1, 4, 5.0}, {2, 1, 2.0}, {3, 1, 4.0}});
  WeightedDigraph graph = WeightedDigraph::from_edges(edges, /*n=*/5);

  // The values below can be verified as being stationary for the following
  // equations:
  //
  // P[0] = (0.15 + 0.85 * P[4]) / 5
  // P[1] = (0.15 + 0.85 * P[4]) / 5
  //        + 0.85 * (P[0] * 1/1 + P[2] * 2/2 + P[3] * 4/4)
  // P[2] = (0.15 + 0.85 * P[4]) / 5
  // P[3] = (0.15 + 0.85 * P[4]) / 5 + 0.85 * (P[1] * 3/8)
  // P[4] = (0.15 + 0.85 * P[4]) / 5 + 0.85 * (P[1] * 5/8)
  EXPECT_THAT(PageRank_edgeMapFunctor::compute_pagerank(graph),
              IsOkAndHolds(Pointwise(
                  DoubleNear(1e-4),
                  {2333.0 / 30348.0, 11360.0 / 30348.0, 2333.0 / 30348.0,
                   5954.0 / 30348.0, 8368.0 / 30348.0})));
}

// Tests the case where a node (in this case, node 2) does have some outgoing
// edges, but they all have weight zero. An implementation that doesn't handle
// that case carefully could end up producing NaNs due to computing the ratio
// 0/0.
TEST_F(PageRankEdgeMapFunctorTest, NodeWithAllOutgoingEdgesHavingZeroWeight) {
  // Graph diagram (each number wrapped in parentheses denotes the weight of its
  // nearest edge):
  //
  //        (1)
  //     0 -----> 3
  //     |        ^
  // (2) |        | (0)
  //     v  (1)   |
  //     1 -----> 2
  //       <-----
  //         (0)
  //

  using WeightedDigraph = asymmetric_graph<asymmetric_vertex, float>;
  sequence<WeightedDigraph::edge> edges(
      {{0, 1, 2.0}, {0, 3, 1.0}, {1, 2, 1.0}, {2, 1, 0.0}, {2, 3, 0.0}});
  WeightedDigraph graph = WeightedDigraph::from_edges(edges, /*n=*/4);

  // The values below can be verified as being stationary for the following
  // equations:
  //
  // P[0] = (0.15 + 0.85 * P[2] + 0.85 * P[3]) / 4
  // P[1] = (0.15 + 0.85 * P[2] + 0.85 * P[3]) / 4 + 0.85 * (P[0] * 2/3)
  // P[2] = (0.15 + 0.85 * P[2] + 0.85 * P[3]) / 4 + 0.85 * (P[1] * 1/1)
  // P[3] = (0.15 + 0.85 * P[2] + 0.85 * P[3]) / 4 + 0.85 * (P[0] * 1/3)
  EXPECT_THAT(PageRank_edgeMapFunctor::compute_pagerank(graph),
              IsOkAndHolds(Pointwise(DoubleNear(1e-4),
                                     {600.0 / 3709.0, 940.0 / 3709.0,
                                      1399.0 / 3709.0, 770.0 / 3709.0})));
}

template <typename T>
using PageRankDeathTest = PageRankFixture<T>;

TYPED_TEST_SUITE(PageRankDeathTest, Implementations);

TYPED_TEST(PageRankDeathTest, EpsNegative) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_DEATH(auto result = Impl::compute_pagerank(graph, /*eps=*/-1.0),
               "Failed assertion `eps >= 0.0`");
}

TYPED_TEST(PageRankDeathTest, DampingFactorNegative) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_DEATH(
      auto result =
          Impl::compute_pagerank(graph,
                                 /*eps=*/0.0, {}, /*damping_factor=*/-1.0),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

TYPED_TEST(PageRankDeathTest, DampingFactorEqualToOne) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_DEATH(
      auto result =
          Impl::compute_pagerank(graph,
                                 /*eps=*/0.0, {}, /*damping_factor=*/1.0),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

TYPED_TEST(PageRankDeathTest, DampingFactorGreaterThanOne) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_DEATH(
      auto result =
          Impl::compute_pagerank(graph,
                                 /*eps=*/0.0, {}, /*damping_factor=*/1.1),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

}  // namespace gbbs