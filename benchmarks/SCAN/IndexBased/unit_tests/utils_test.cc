#include "benchmarks/SCAN/IndexBased/utils.h"

#include <unordered_set>

#include "benchmarks/SCAN/IndexBased/utils.h"
#include "gtest/gtest.h"
#include "ligra/graph.h"
#include "ligra/graph_test_utils.h"
#include "pbbslib/seq.h"

namespace gt = graph_test;

TEST(Modularity, NullGraph) {
  constexpr size_t kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  const scan::Clustering clusters{};
  constexpr uintE kMaxClusterId{0};

  EXPECT_EQ(scan::Modularity(&graph, clusters, kMaxClusterId), 0.0);
}

TEST(Modularity, EmptyGraph) {
  constexpr size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  {
    const scan::Clustering clusters{0, 1, 2, 3, 4, 5};
    constexpr uintE kMaxClusterId{kNumVertices};
    EXPECT_EQ(scan::Modularity(&graph, clusters, kMaxClusterId), 0.0);
  }
  {
    const scan::Clustering clusters(kNumVertices, scan::kUnclustered);
    constexpr uintE kMaxClusterId{0};
    EXPECT_EQ(scan::Modularity(&graph, clusters, kMaxClusterId), 0.0);
  }
}
