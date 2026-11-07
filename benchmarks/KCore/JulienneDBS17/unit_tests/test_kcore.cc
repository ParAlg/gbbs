#include "benchmarks/KCore/JulienneDBS17/KCore.h"

#include <unordered_set>

#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::AnyOf;
using ::testing::ElementsAre;

namespace gbbs {

TEST(KCore, EdgelessGraph) {
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const sequence<uintE> kcoreResult{KCore(graph)};
  EXPECT_THAT(kcoreResult, ElementsAre(0, 0, 0));
}

TEST(KCore, BasicUsage) {
  // Graph diagram:
  //     0 - 1    2 - 3 - 4
  //                    \ |
  //                      5 -- 6    7
  constexpr uintE kNumVertices{8};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {2, 3}, {3, 4}, {3, 5}, {4, 5}, {5, 6},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const sequence<uintE> kcoreResult{KCore(graph)};
  EXPECT_THAT(kcoreResult, ElementsAre(1, 1, 1, 2, 2, 2, 1, 0));
}


TEST(KCore, FetchAddReportsFirstContributionOnly) {
  uintE pending[] = {0};
  uintE degree[] = {3};
  kcore_fetch_add<empty> update(pending, degree, 1);
  EXPECT_TRUE(update.updateAtomic(0,0,empty{}).has_value());
  EXPECT_FALSE(update.updateAtomic(0,0,empty{}).has_value());
  EXPECT_EQ(pending[0],2u);
}

TEST(KCore, FetchAddMatchesOrdinaryPeelingOnPath) {
  const std::unordered_set<UndirectedEdge> edges{{0,1},{1,2}};
  auto graph = graph_test::MakeUnweightedSymmetricGraph(3,edges);
  auto result = KCore_FA(graph);
  EXPECT_THAT(result,ElementsAre(1,1,1));
}

}  // namespace gbbs
