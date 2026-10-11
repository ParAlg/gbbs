#include "benchmarks/GraphColoring/Hasenplaugh14/GraphColoring.h"

#include <cstddef>
#include <set>
#include <tuple>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gbbs/bridge.h"
#include "gbbs/graph.h"
#include "gbbs/graph_generator.h"
#include "gbbs/macros.h"
#include "gbbs/vertex.h"
#include "parlay/monoid.h"
#include "parlay/primitives.h"

namespace gbbs {
namespace {

using ::testing::Bool;
using ::testing::Contains;
using ::testing::ElementsAre;
using ::testing::TestWithParam;

class GraphColoringTest : public TestWithParam<bool> {
 protected:
  bool use_lf_heuristic() const { return GetParam(); }
};

TEST_P(GraphColoringTest, SingleVertexColoring) {
  auto graph = graph_generator::clique_graph(1);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_EQ(colors[0], 0);
}

TEST_P(GraphColoringTest, TwoConnectedVertices) {
  auto graph = graph_generator::clique_graph(2);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_NE(colors[0], colors[1]);
  EXPECT_LE(parlay::reduce_max(colors), 1);
}

TEST_P(GraphColoringTest, CliqueColoring) {
  uintE n = 10;
  auto graph = graph_generator::clique_graph(n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_EQ(parlay::reduce_max(colors), n - 1);
}

TEST_P(GraphColoringTest, PathColoring) {
  uintE n = 10;
  auto graph = graph_generator::path_graph(n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_LE(parlay::reduce_max(colors), 2);
}

TEST_P(GraphColoringTest, EvenCycleColoring) {
  uintE n = 10;
  auto graph = graph_generator::cycle_graph(n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_LE(parlay::reduce_max(colors), 2);
}

TEST_P(GraphColoringTest, OddCycleColoring) {
  uintE n = 9;
  auto graph = graph_generator::cycle_graph(n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_LE(parlay::reduce_max(colors), 2);
}

TEST_P(GraphColoringTest, StarGraphColoring) {
  uintE n = 11;
  size_t m = 2 * (n - 1);
  auto edges = sequence<std::tuple<uintE, uintE, empty>>(m);
  for (uintE i = 1; i < n; ++i) {
    edges[2 * (i - 1)] = std::make_tuple(0, i, empty{});
    edges[2 * (i - 1) + 1] = std::make_tuple(i, 0, empty{});
  }
  auto graph = symmetric_graph<symmetric_vertex, empty>::from_edges(edges, n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_LE(parlay::reduce_max(colors), 1);
}

TEST_P(GraphColoringTest, ColorsStartFromZero) {
  uintE n = 20;
  auto graph = graph_generator::path_graph(n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  uintE min_color = parlay::reduce(colors, parlay::minimum<uintE>());
  EXPECT_EQ(min_color, 0);
}

TEST_P(GraphColoringTest, ColorsAreContiguous) {
  uintE n = 15;
  auto graph = graph_generator::clique_graph(n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);

  std::set<uintE> unique_colors(colors.begin(), colors.end());
  uintE max_color = parlay::reduce_max(colors);
  EXPECT_EQ(unique_colors.size(), max_color + 1);
  for (uintE c = 0; c <= max_color; c++) {
    EXPECT_THAT(unique_colors, Contains(c));
  }
}

TEST_P(GraphColoringTest, LargerPathGraph) {
  uintE n = 1000;
  auto graph = graph_generator::path_graph(n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_LE(parlay::reduce_max(colors), 2);
}

TEST_P(GraphColoringTest, LargerCliqueGraph) {
  uintE n = 50;
  auto graph = graph_generator::clique_graph(n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_EQ(parlay::reduce_max(colors), n - 1);
}

TEST_P(GraphColoringTest, SingleVertexWithSelfLoop) {
  uintE n = 1;
  auto edges = sequence<std::tuple<uintE, uintE, empty>>(1);
  edges[0] = std::make_tuple(0, 0, empty{});
  // Note: from_edges does not remove self-loops.
  auto graph = symmetric_graph<symmetric_vertex, empty>::from_edges(edges, n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_EQ(colors[0], 0);
}

TEST_P(GraphColoringTest, TriangleWithSelfLoops) {
  uintE n = 3;
  auto edges = sequence<std::tuple<uintE, uintE, empty>>(11);
  edges[0] = std::make_tuple(0, 1, empty{});
  edges[1] = std::make_tuple(1, 0, empty{});
  edges[2] = std::make_tuple(0, 2, empty{});
  edges[3] = std::make_tuple(2, 0, empty{});
  edges[4] = std::make_tuple(1, 2, empty{});
  edges[5] = std::make_tuple(2, 1, empty{});
  edges[6] = std::make_tuple(0, 0, empty{});
  edges[7] = std::make_tuple(0, 0, empty{});
  edges[8] = std::make_tuple(0, 0, empty{});
  edges[9] = std::make_tuple(1, 1, empty{});
  edges[10] = std::make_tuple(2, 2, empty{});
  // Note: from_edges does not remove self-loops.
  auto graph = symmetric_graph<symmetric_vertex, empty>::from_edges(edges, n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_EQ(parlay::reduce_max(colors), 2);
}

TEST_P(GraphColoringTest, DisconnectedVerticesWithSelfLoops) {
  uintE n = 5;
  auto edges = sequence<std::tuple<uintE, uintE, empty>>(2);
  edges[0] = std::make_tuple(1, 1, empty{});
  edges[1] = std::make_tuple(3, 3, empty{});
  // Note: from_edges does not remove self-loops.
  auto graph = symmetric_graph<symmetric_vertex, empty>::from_edges(edges, n);
  auto colors = Coloring(graph, use_lf_heuristic());
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  for (uintE i = 0; i < n; ++i) {
    EXPECT_EQ(colors[i], 0);
  }
}

INSTANTIATE_TEST_SUITE_P(GraphColoringTestSuite, GraphColoringTest, Bool(),
                         [](const ::testing::TestParamInfo<bool>& info) {
                           return info.param ? "LF" : "LLF";
                         });

using DeterministicColoringTest = ::testing::Test;

TEST_F(DeterministicColoringTest, CliqueGraphExactColorsLF) {
  uintE n = 4;
  auto graph = graph_generator::clique_graph(n);
  auto random_ranks = sequence<uintE>::from_function(
      n, [](size_t i) { return static_cast<uintE>(i); });
  auto colors =
      coloring::ColoringWithRanks(graph, /*use_lf=*/true, random_ranks);
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_THAT(colors, ElementsAre(0, 1, 2, 3));
}

TEST_F(DeterministicColoringTest, CliqueGraphReversedRanks) {
  uintE n = 4;
  auto graph = graph_generator::clique_graph(n);
  auto random_ranks = sequence<uintE>::from_function(
      n, [n](size_t i) { return static_cast<uintE>(n - 1 - i); });
  auto colors =
      coloring::ColoringWithRanks(graph, /*use_lf=*/true, random_ranks);
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_THAT(colors, ElementsAre(3, 2, 1, 0));
}

TEST_F(DeterministicColoringTest, StarGraphCenterFirst) {
  uintE n = 5;
  size_t m = 2 * (n - 1);
  auto edges = sequence<std::tuple<uintE, uintE, empty>>(m);
  for (uintE i = 1; i < n; ++i) {
    edges[2 * (i - 1)] = std::make_tuple(0, i, empty{});
    edges[2 * (i - 1) + 1] = std::make_tuple(i, 0, empty{});
  }
  auto graph = symmetric_graph<symmetric_vertex, empty>::from_edges(edges, n);
  auto random_ranks = sequence<uintE>::from_function(
      n, [](size_t i) { return static_cast<uintE>(i); });
  auto colors =
      coloring::ColoringWithRanks(graph, /*use_lf=*/true, random_ranks);
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_THAT(colors, ElementsAre(0, 1, 1, 1, 1));
}

TEST_F(DeterministicColoringTest, CycleGraphExactColors) {
  uintE n = 6;
  auto graph = graph_generator::cycle_graph(n);
  auto random_ranks = sequence<uintE>::from_function(
      n, [](size_t i) { return static_cast<uintE>(i); });
  auto colors =
      coloring::ColoringWithRanks(graph, /*use_lf=*/true, random_ranks);
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_THAT(colors, ElementsAre(0, 1, 0, 1, 0, 1));
}

TEST_F(DeterministicColoringTest, OddCycleRequiresThreeColors) {
  uintE n = 5;
  auto graph = graph_generator::cycle_graph(n);
  auto random_ranks = sequence<uintE>::from_function(
      n, [](size_t i) { return static_cast<uintE>(i); });
  auto colors =
      coloring::ColoringWithRanks(graph, /*use_lf=*/true, random_ranks);
  EXPECT_EQ(coloring::CountColoringViolations(graph, colors), 0);
  EXPECT_THAT(colors, ElementsAre(0, 1, 0, 1, 2));
}

TEST_F(DeterministicColoringTest, VerifyPriorityComputationCliqueReversed) {
  uintE n = 4;
  auto graph = graph_generator::clique_graph(n);
  auto random_ranks = sequence<uintE>::from_function(
      n, [n](size_t i) { return static_cast<uintE>(n - 1 - i); });
  auto priorities = coloring::ComputeVertexPrioritiesWithRanks(
      graph, /*use_lf=*/true, random_ranks);
  EXPECT_THAT(priorities, ElementsAre(3, 2, 1, 0));
}

TEST_F(DeterministicColoringTest, VerifyPriorityComputationClique) {
  uintE n = 4;
  auto graph = graph_generator::clique_graph(n);
  auto random_ranks = sequence<uintE>::from_function(
      n, [](size_t i) { return static_cast<uintE>(i); });
  auto priorities = coloring::ComputeVertexPrioritiesWithRanks(
      graph, /*use_lf=*/true, random_ranks);
  EXPECT_THAT(priorities, ElementsAre(0, 1, 2, 3));
}

TEST_F(DeterministicColoringTest, DeterministicResultsAreReproducible) {
  uintE n = 10;
  auto graph = graph_generator::clique_graph(n);
  auto random_ranks = sequence<uintE>::from_function(
      n, [](size_t i) { return static_cast<uintE>(i); });
  auto colors1 =
      coloring::ColoringWithRanks(graph, /*use_lf=*/true, random_ranks);
  auto colors2 =
      coloring::ColoringWithRanks(graph, /*use_lf=*/true, random_ranks);
  for (uintE i = 0; i < n; ++i) {
    EXPECT_EQ(colors1[i], colors2[i]);
  }
}

TEST_F(DeterministicColoringTest, DifferentRanksProduceDifferentColors) {
  uintE n = 4;
  auto graph = graph_generator::clique_graph(n);
  auto ranks_forward = sequence<uintE>::from_function(
      n, [](size_t i) { return static_cast<uintE>(i); });
  auto ranks_reversed = sequence<uintE>::from_function(
      n, [n](size_t i) { return static_cast<uintE>(n - 1 - i); });
  auto colors_forward =
      coloring::ColoringWithRanks(graph, /*use_lf=*/true, ranks_forward);
  auto colors_reversed =
      coloring::ColoringWithRanks(graph, /*use_lf=*/true, ranks_reversed);
  bool any_different = false;
  for (uintE i = 0; i < n; ++i) {
    if (colors_forward[i] != colors_reversed[i]) {
      any_different = true;
      break;
    }
  }
  EXPECT_TRUE(any_different);
}

}  // namespace
}  // namespace gbbs
