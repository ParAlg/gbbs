#include "benchmarks/MaximalIndependentSet/RandomGreedy/MaximalIndependentSet.h"

#include <unordered_set>

#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::AnyOf;
using ::testing::UnorderedElementsAre;

namespace gbbs {

using namespace MaximalIndependentSet_rootset;

TEST(MaximalIndependentSet, EdgelessGraph) {
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const sequence<bool> misResult{MaximalIndependentSet(graph)};
  EXPECT_THAT(misResult, UnorderedElementsAre(true, true, true));
}

template <class Graph>
inline void CheckMIS(Graph& G, parlay::sequence<bool>& mis) {
  using W = typename Graph::weight_type;
  auto d = sequence<uintE>(G.n, (uintE)0);
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    if (!d[ngh]) {
      d[ngh] = 1;
    }
  };
  parallel_for(0, G.n, [&](size_t i) {
    if (mis[i]) {
      G.get_vertex(i).out_neighbors().map(map_f);
    }
  });
  parallel_for(0, G.n, [&](size_t i) {
    if (mis[i]) {
      EXPECT_FALSE(d[i]);
    }
  });

  auto mis_int = parlay::delayed_seq<size_t>(
      G.n, [&](size_t i) { return (size_t)mis[i]; });
  size_t mis_size = parlay::reduce(mis_int);
  EXPECT_EQ(parlay::reduce(d), (G.n - mis_size));
}

TEST(MaximalIndependentSet, BasicUsage) {
  // Graph diagram:
  //     0 - 1    2 - 3 - 4
  //                    \ |
  //                      5 -- 6    7
  constexpr uintE kNumVertices{8};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {2, 3}, {3, 4}, {3, 5}, {4, 5}, {5, 6},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  sequence<bool> misResult{MaximalIndependentSet(graph)};
  CheckMIS(graph, misResult);
}

}  // namespace gbbs
