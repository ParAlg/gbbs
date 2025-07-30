#include "benchmarks/PageRank/PageRank.h"

#include <cstddef>
#include <unordered_set>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "benchmarks/PageRank/PageRank_edgeMapReduce.h"
#include "gbbs/bridge.h"
#include "gbbs/helpers/undirected_edge.h"
#include "gbbs/macros.h"
#include "gbbs/unit_tests/graph_test_utils.h"

using ::testing::ElementsAre;
using ::testing::IsEmpty;

namespace gbbs {

// TODO: add tests for directed graphs.

struct PageRank_ligra {
  template <class Graph>
  static sequence<double> compute_pagerank(Graph& G, double eps = 0.000001,
                                           std::vector<uintE> sources = {},
                                           double damping_factor = 0.85,
                                           size_t max_iters = 100) {
    return PageRank_edgeMap(G, eps, sources, damping_factor, max_iters);
  }
};

struct PageRank_opt {
  template <class Graph>
  static sequence<double> compute_pagerank(Graph& G, double eps = 0.000001,
                                           std::vector<uintE> sources = {},
                                           double damping_factor = 0.85,
                                           size_t max_iters = 100) {
    return PageRank_edgeMapReduce(G, eps, sources, damping_factor, max_iters);
  }
};

template <typename T>
class PageRankFixture : public testing::Test {
 public:
  using Impl = T;
};

using Implementations = ::testing::Types<PageRank_ligra, PageRank_opt>;
TYPED_TEST_SUITE(PageRankFixture, Implementations);

TYPED_TEST(PageRankFixture, EmptyGraph) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_THAT(Impl::compute_pagerank(graph), IsEmpty());
}

TYPED_TEST(PageRankFixture, EdgelessGraph) {
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  const sequence<double> result = Impl::compute_pagerank(graph);
  EXPECT_THAT(result, ElementsAre(1.0 / 3, 1.0 / 3, 1.0 / 3));
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
                                     /*damping_factor=*/0.85,
                                     /*max_iters=*/0),
              ElementsAre(0.25, 0.25, 0.25, 0.25));
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
    const sequence<double> result{Impl::compute_pagerank(graph)};
    EXPECT_THAT(result, ElementsAre(0.2, 0.2, 0.2, 0.2, 0.2));
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

  {
    const sequence<double> result{Impl::compute_pagerank(graph)};
    const sequence<double> expected{0.2567570878, 0.4864858243, 0.2567570878};
    EXPECT_THAT(result,
                testing::Pointwise(testing::DoubleNear(1e-4), expected));
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

  {
    const sequence<double> result{Impl::compute_pagerank(graph)};
    const sequence<double> expected{0.1428571429, 0.1428571429, 0.0802049539,
                                    0.2074472351, 0.1389813363, 0.2074472351,
                                    0.0802049539};

    EXPECT_THAT(result,
                testing::Pointwise(testing::DoubleNear(1e-4), expected));
  }
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

  {
    const sequence<double> result{
        Impl::compute_pagerank(graph, /*eps=*/0.000001, {0, 4})};
    // The results below were computed by using NetworkX's source-based PageRank
    // as a reference implementation.
    const sequence<double> expected{0.2702716450, 0.2297283550, 0.0384308511,
                                    0.1356382979, 0.1518617021, 0.1356382979,
                                    0.0384308511};

    EXPECT_THAT(result,
                testing::Pointwise(testing::DoubleNear(1e-4), expected));
  }
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

  {
    const sequence<double> result{
        Impl::compute_pagerank(graph, /*eps=*/0.000001, {0, 1, 2, 4})};
    const sequence<double> expected{0.2137115266, 0.2137115266, 0.2009045292,
                                    0.1707678884, 0.2009045292};

    EXPECT_THAT(result,
                testing::Pointwise(testing::DoubleNear(1e-4), expected));
  }
}

TYPED_TEST(PageRankFixture, TwoNodesNoEdges) {
  constexpr uintE kNumVertices{2};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  using Impl = typename TestFixture::Impl;

  {
    const sequence<double> result{
        Impl::compute_pagerank(graph, /*eps=*/0.000001, {0, 1})};
    // No mass is migrated.
    const sequence<double> expected{0.5, 0.5};

    EXPECT_THAT(result,
                testing::Pointwise(testing::DoubleNear(1e-4), expected));
  }

  {
    const sequence<double> result{
        Impl::compute_pagerank(graph, /*eps=*/0.000001, {0})};
    // No mass is migrated.
    const sequence<double> expected{1.0, 0};

    EXPECT_THAT(result,
                testing::Pointwise(testing::DoubleNear(1e-4), expected));
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
    const sequence<double> result{
        Impl::compute_pagerank(graph, /*eps=*/0.000001, {0, 1})};
    // Both have equal mass; again no mass gets migrated.
    const sequence<double> expected{0.5, 0.5};

    EXPECT_THAT(result,
                testing::Pointwise(testing::DoubleNear(1e-4), expected));
  }

  {
    const sequence<double> result{
        Impl::compute_pagerank(graph, /*eps=*/0.000001, {0})};

    // The steady state for this setting is when P[0] = 20/37 and P[1] = 17/37.
    // This can be verified by plugging in the values above into the PPR
    // equations
    // P[0] = 0.85*P[1] + 0.15
    // P[1] = 0.85*P[0]
    const sequence<double> expected{0.5405409317, 0.4594590684};

    EXPECT_THAT(result,
                testing::Pointwise(testing::DoubleNear(1e-4), expected));
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
    const sequence<double> result{
        Impl::compute_pagerank(graph, /*eps=*/0.000001, {0, 1, 2, 3})};

    // The values below can be verified as being stationary for the following
    // equations (note that with only one source, the addedConstant is only
    // present for node 0)
    // P[0] = 0.85*P[1]/2 + 0.15/4
    // P[1] = 0.85*(P[0] + P[2]/2) + 0.15/4
    // P[2] = 0.85*(P[1]/2 + P[3]) + 0.15/4
    // P[3] = 0.85*P[2]/2 + 0.15/4
    const sequence<double> expected{0.17543856058862931, 0.32456143941137061,
                                    0.32456143941137061, 0.17543856058862931};

    EXPECT_THAT(result,
                testing::Pointwise(testing::DoubleNear(1e-4), expected));
  }

  {
    const sequence<double> result{
        Impl::compute_pagerank(graph, /*eps=*/0.000001, {0})};

    // The values below can be verified as being stationary for the following
    // equations (note that with only one source, the addedConstant is only
    // present for node 0)
    // P[0] = 0.85*P[1]/2 + 0.15
    // P[1] = 0.85*(P[0] + P[2]/2)
    // P[2] = 0.85*(P[1]/2 + P[3])
    // P[3] = 0.85*P[2]/2
    const sequence<double> expected{0.3022241274, 0.3581756964, 0.2383155317,
                                    0.1012846445};

    EXPECT_THAT(result,
                testing::Pointwise(testing::DoubleNear(1e-4), expected));
  }

  {
    const sequence<double> result{
        Impl::compute_pagerank(graph, /*eps=*/0.000001, {0, 3})};

    // The values below can be verified as being stationary for the following
    // equations (now that we have two sources, the added constant is split
    // between the two sources).
    // P[0] = 0.85*P[1]/2 + 0.15/2
    // P[1] = 0.85*(P[0] + P[2]/2)
    // P[2] = 0.85*(P[1]/2 + P[3])
    // P[3] = 0.85*P[2]/2 + 0.15/2
    const sequence<double> expected{0.2017542424, 0.2982457576, 0.2982457576,
                                    0.2017542424};

    EXPECT_THAT(result,
                testing::Pointwise(testing::DoubleNear(1e-4), expected));
  }
}

template <typename T>
using PageRankDeathTest = PageRankFixture<T>;

TYPED_TEST_SUITE(PageRankDeathTest, Implementations);

TYPED_TEST(PageRankDeathTest, EpsNegative) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_DEATH(Impl::compute_pagerank(graph, /*eps=*/-1.0),
               "Failed assertion `eps >= 0.0`");
}

TYPED_TEST(PageRankDeathTest, DampingFactorNegative) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_DEATH(
      Impl::compute_pagerank(graph, /*eps=*/0.0, {}, /*damping_factor=*/-1.0),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

TYPED_TEST(PageRankDeathTest, DampingFactorEqualToOne) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_DEATH(
      Impl::compute_pagerank(graph, /*eps=*/0.0, {}, /*damping_factor=*/1.0),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

TYPED_TEST(PageRankDeathTest, DampingFactorGreaterThanOne) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_DEATH(
      Impl::compute_pagerank(graph, /*eps=*/0.0, {}, /*damping_factor=*/1.1),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

}  // namespace gbbs