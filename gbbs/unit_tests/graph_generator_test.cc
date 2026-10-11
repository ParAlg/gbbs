#include "gbbs/graph_generator.h"

#include <set>
#include <utility>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gbbs/macros.h"

namespace gbbs {
namespace graph_generator {
namespace {

using ::testing::IsEmpty;
using ::testing::Pair;
using ::testing::SizeIs;
using ::testing::UnorderedElementsAre;

TEST(CliqueEdgesTest, EmptyClique) {
  auto edges = clique_edges(0);
  EXPECT_THAT(edges, IsEmpty());
}

TEST(CliqueEdgesTest, SingleVertex) {
  auto edges = clique_edges(1);
  EXPECT_THAT(edges, IsEmpty());
}

TEST(CliqueEdgesTest, TwoVertices) {
  auto edges = clique_edges(2);
  std::set<std::pair<uintE, uintE>> edge_set;
  for (const auto& [u, v, w] : edges) {
    edge_set.insert({u, v});
  }
  EXPECT_THAT(edge_set, UnorderedElementsAre(Pair(0, 1), Pair(1, 0)));
}

TEST(CliqueEdgesTest, FourVertices) {
  auto edges = clique_edges(4);
  std::set<std::pair<uintE, uintE>> edge_set;
  for (const auto& [u, v, w] : edges) {
    edge_set.insert({u, v});
  }
  EXPECT_THAT(edge_set, UnorderedElementsAre(
                            Pair(0, 1), Pair(0, 2), Pair(0, 3), Pair(1, 0),
                            Pair(1, 2), Pair(1, 3), Pair(2, 0), Pair(2, 1),
                            Pair(2, 3), Pair(3, 0), Pair(3, 1), Pair(3, 2)));
}

TEST(PathEdgesTest, EmptyPath) {
  auto edges = path_edges(0);
  EXPECT_THAT(edges, IsEmpty());
}

TEST(PathEdgesTest, SingleVertex) {
  auto edges = path_edges(1);
  EXPECT_THAT(edges, IsEmpty());
}

TEST(PathEdgesTest, TwoVertices) {
  auto edges = path_edges(2);
  std::set<std::pair<uintE, uintE>> edge_set;
  for (const auto& [u, v, w] : edges) {
    edge_set.insert({u, v});
  }
  EXPECT_THAT(edge_set, UnorderedElementsAre(Pair(0, 1), Pair(1, 0)));
}

TEST(PathEdgesTest, FiveVertices) {
  auto edges = path_edges(5);
  std::set<std::pair<uintE, uintE>> edge_set;
  for (const auto& [u, v, w] : edges) {
    edge_set.insert({u, v});
  }
  EXPECT_THAT(edge_set, UnorderedElementsAre(Pair(0, 1), Pair(1, 0), Pair(1, 2),
                                             Pair(2, 1), Pair(2, 3), Pair(3, 2),
                                             Pair(3, 4), Pair(4, 3)));
}

TEST(CycleEdgesTest, EmptyCycle) {
  auto edges = cycle_edges(0);
  EXPECT_THAT(edges, IsEmpty());
}

TEST(CycleEdgesTest, SingleVertex) {
  auto edges = cycle_edges(1);
  EXPECT_THAT(edges, IsEmpty());
}

TEST(CycleEdgesTest, TwoVertices) {
  auto edges = cycle_edges(2);
  EXPECT_THAT(edges, SizeIs(4));
}

TEST(CycleEdgesTest, FiveVertices) {
  auto edges = cycle_edges(5);
  std::set<std::pair<uintE, uintE>> edge_set;
  for (const auto& [u, v, w] : edges) {
    edge_set.insert({u, v});
  }
  EXPECT_THAT(edge_set, UnorderedElementsAre(Pair(0, 1), Pair(1, 0), Pair(1, 2),
                                             Pair(2, 1), Pair(2, 3), Pair(3, 2),
                                             Pair(3, 4), Pair(4, 3), Pair(4, 0),
                                             Pair(0, 4)));
}

TEST(CliqueGraphTest, SingleVertex) {
  auto graph = clique_graph(1);
  EXPECT_EQ(graph.n, 1);
  EXPECT_EQ(graph.get_vertex(0).out_degree(), 0);
}

TEST(CliqueGraphTest, FourVertices) {
  auto graph = clique_graph(4);
  EXPECT_EQ(graph.n, 4);
  for (uintE i = 0; i < 4; ++i) {
    EXPECT_EQ(graph.get_vertex(i).out_degree(), 3);
  }
}

TEST(PathGraphTest, SingleVertex) {
  auto graph = path_graph(1);
  EXPECT_EQ(graph.n, 1);
  EXPECT_EQ(graph.get_vertex(0).out_degree(), 0);
}

TEST(PathGraphTest, FiveVertices) {
  auto graph = path_graph(5);
  EXPECT_EQ(graph.n, 5);
  EXPECT_EQ(graph.get_vertex(0).out_degree(), 1);
  EXPECT_EQ(graph.get_vertex(1).out_degree(), 2);
  EXPECT_EQ(graph.get_vertex(2).out_degree(), 2);
  EXPECT_EQ(graph.get_vertex(3).out_degree(), 2);
  EXPECT_EQ(graph.get_vertex(4).out_degree(), 1);
}

TEST(CycleGraphTest, SingleVertex) {
  auto graph = cycle_graph(1);
  EXPECT_EQ(graph.n, 1);
  EXPECT_EQ(graph.get_vertex(0).out_degree(), 0);
}

TEST(CycleGraphTest, FiveVertices) {
  auto graph = cycle_graph(5);
  EXPECT_EQ(graph.n, 5);
  for (uintE i = 0; i < 5; ++i) {
    EXPECT_EQ(graph.get_vertex(i).out_degree(), 2);
  }
}

TEST(CycleGraphTest, ThreeVertices) {
  auto graph = cycle_graph(3);
  EXPECT_EQ(graph.n, 3);
  for (uintE i = 0; i < 3; ++i) {
    EXPECT_EQ(graph.get_vertex(i).out_degree(), 2);
  }
}

}  // namespace
}  // namespace graph_generator
}  // namespace gbbs
