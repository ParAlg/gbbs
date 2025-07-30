#include "benchmarks/PageRank/PageRank.h"

#include <cstddef>
#include <unordered_set>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
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
                                           double damping_factor = 0.85,
                                           size_t max_iters = 100) {
    return PageRank_edgeMap(G, eps, damping_factor, max_iters);
  }
};

struct PageRank_opt {
  template <class Graph>
  static sequence<double> compute_pagerank(Graph& G, double eps = 0.000001,
                                           double damping_factor = 0.85,
                                           size_t max_iters = 100) {
    return PageRank_edgeMapReduce(G, eps, damping_factor, max_iters);
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
  EXPECT_THAT(
      Impl::compute_pagerank(graph, /*eps=*/0.1, /*damping_factor=*/0.85,
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
      Impl::compute_pagerank(graph, /*eps=*/0.0, /*damping_factor=*/-1.0),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

TYPED_TEST(PageRankDeathTest, DampingFactorEqualToOne) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_DEATH(
      Impl::compute_pagerank(graph, /*eps=*/0.0, /*damping_factor=*/1.0),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

TYPED_TEST(PageRankDeathTest, DampingFactorGreaterThanOne) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  using Impl = typename TestFixture::Impl;
  EXPECT_DEATH(
      Impl::compute_pagerank(graph, /*eps=*/0.0, /*damping_factor=*/1.1),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

}  // namespace gbbs