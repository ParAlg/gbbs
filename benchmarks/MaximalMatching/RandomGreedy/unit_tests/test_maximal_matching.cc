#include "benchmarks/MaximalMatching/RandomGreedy/MaximalMatching.h"

#include <unordered_set>

#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::AnyOf;
using ::testing::ElementsAre;

namespace gbbs {

TEST(MaximalMatching, EdgelessGraph) {
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const sequence<std::tuple<uintE, uintE, gbbs::empty>> matchingResult{
      MaximalMatching(graph)};
  EXPECT_EQ(matchingResult.size(), 0);
}

template <template <class W> class vertex, class W, class Seq>
inline void CheckMatching(symmetric_graph<vertex, W>& G, Seq& matching) {
  size_t n = G.n;
  auto matched = sequence<uintE>::from_function(n, [](size_t i) { return 0; });

  // Check that this is a valid matching
  parallel_for(0, matching.size(), [&](size_t i) {
    const auto& edge = matching[i];
    gbbs::write_add(&matched[std::get<0>(edge)], 1);
    gbbs::write_add(&matched[std::get<1>(edge)], 1);
  });

  parallel_for(0, n, [&](size_t i) { EXPECT_LE(matched[i], 1); });

  // Check maximality of the matching
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    EXPECT_FALSE(!matched[src] && !matched[ngh]);
  };
  parallel_for(0, n,
               [&](size_t i) { G.get_vertex(i).out_neighbors().map(map_f); });
}

TEST(MaximalMatching, BasicUsage) {
  // Graph diagram:
  //     0 - 1    2 - 3 - 4
  //                    \ |
  //                      5 -- 6
  constexpr uintE kNumVertices{7};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {2, 3}, {3, 4}, {3, 5}, {4, 5}, {5, 6},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
  auto graph_copy = graph;

  const sequence<std::tuple<uintE, uintE, gbbs::empty>> matchingResult{
      MaximalMatching(graph)};
  for (size_t i = 0; i < matchingResult.size(); i++) {
    std::cout << std::get<0>(matchingResult[i]) << " "
              << std::get<1>(matchingResult[i]) << std::endl;
  }
  CheckMatching(graph_copy, matchingResult);
}

}  // namespace gbbs
