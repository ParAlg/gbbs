#include "benchmarks/BFS/NonDeterministicBFS/BFS.h"

#include <unordered_set>

#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::AnyOf;
using ::testing::ElementsAre;

namespace gbbs {

TEST(NondeterministicBFS, EdgelessGraph) {
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  constexpr uintE source_vertex{1};
  const sequence<uintE> bfsResult{BFS(graph, source_vertex)};
  EXPECT_THAT(bfsResult, ElementsAre(UINT_E_MAX, source_vertex, UINT_E_MAX));
}

TEST(NondeterministicBFS, BasicUsage) {
  // Graph diagram:
  //     0 - 1    2 - 3 - 4
  //                    \ |
  //                      5 -- 6
  constexpr uintE kNumVertices{7};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {2, 3}, {3, 4}, {3, 5}, {4, 5}, {5, 6},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  {
    constexpr uintE source_vertex{1};
    const sequence<uintE> bfsResult{BFS(graph, source_vertex)};
    EXPECT_THAT(bfsResult,
                ElementsAre(source_vertex, source_vertex, UINT_E_MAX,
                            UINT_E_MAX, UINT_E_MAX, UINT_E_MAX, UINT_E_MAX));
  }
  {
    constexpr uintE source_vertex{2};
    const sequence<uintE> bfsResult{BFS(graph, source_vertex)};
    EXPECT_THAT(bfsResult,
                ElementsAre(UINT_E_MAX, UINT_E_MAX, source_vertex,
                            source_vertex, AnyOf(3, 5), AnyOf(3, 4), 5));
  }
}

}  // namespace gbbs
