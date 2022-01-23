#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"

#include <unordered_set>

#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::AnyOf;
using ::testing::ElementsAre;

namespace gbbs {

TEST(WorkEfficientCC, EdgelessGraph) {
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const sequence<uintE> ccResult{workefficient_cc::CC(graph)};
  EXPECT_THAT(ccResult, ElementsAre(0, 1, 2));
}

TEST(WorkEfficientCC, BasicUsage) {
  // Graph diagram:
  //     0 - 1    2 - 3 - 4
  //                    \ |
  //                      5 -- 6
  constexpr uintE kNumVertices{7};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {2, 3}, {3, 4}, {3, 5}, {4, 5}, {5, 6},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const sequence<uintE> ccResult{workefficient_cc::CC(graph)};
  uintE class_one = ccResult[0];
  uintE class_two = ccResult[2];
  EXPECT_THAT(ccResult, ElementsAre(class_one, class_one, class_two, class_two,
                                    class_two, class_two, class_two));
  EXPECT_NE(class_one, class_two);
}

}  // namespace gbbs
