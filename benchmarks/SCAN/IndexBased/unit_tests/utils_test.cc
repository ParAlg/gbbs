#include "benchmarks/SCAN/IndexBased/utils.h"

#include "gbbs/unit_tests/graph_test_utils.h"
#include "gtest/gtest.h"

namespace gbbs {
namespace gt = graph_test;

TEST(Modularity, NullGraph) {
  constexpr size_t kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  const scan::Clustering clusters{};
  EXPECT_EQ(scan::Modularity(&graph, clusters), 0.0);
}

TEST(Modularity, EmptyGraph) {
  constexpr size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  const scan::Clustering clusters{0, 1, 2, 3, 4, 5};
  EXPECT_EQ(scan::Modularity(&graph, clusters), 0.0);
}

TEST(Modularity, BasicGraph) {
  // Graph diagram:
  //   0               5
  //   | \           / |
  //   |  2 -- 3 -- 4  |
  //   | /     |     \ |
  //   1       7       6
  //         /   \.
  //        8 --- 9
  constexpr size_t kNumVertices{10};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {0, 2}, {1, 2}, {2, 3}, {3, 4}, {3, 7},
      {4, 5}, {4, 6}, {5, 6}, {7, 8}, {7, 9}, {8, 9},
  };
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  {
    const scan::Clustering clusters{0, 0, 0, 0, 1, 1, 1, 2, 2, 2};
    EXPECT_FLOAT_EQ(scan::Modularity(&graph, clusters), 0.48958333);
  }
  {  // Test having unclustered vertices (vertices 3, 5, 6)
    const scan::Clustering clusters{
        0, 0, 0, scan::kUnclustered, 1, scan::kUnclustered, scan::kUnclustered,
        2, 2, 2};
    EXPECT_FLOAT_EQ(scan::Modularity(&graph, clusters), 0.28472222);
  }
}
}  // namespace gbbs
