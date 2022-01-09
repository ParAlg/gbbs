#include "gbbs/graph_io.h"

#include <vector>

#include "gbbs/unit_tests/graph_test_utils.h"
#include "gtest/gtest.h"

namespace gbbs {

namespace gi = gbbs_io;
namespace gt = graph_test;
using NoWeight = gbbs::empty;

TEST(EdgeListToAsymmetricGraph, NoEdges) {
  const std::vector<gi::Edge<NoWeight>> kEdges{};
  const auto graph{gi::edge_list_to_asymmetric_graph(kEdges)};
  EXPECT_EQ(graph.n, 0);
  EXPECT_EQ(graph.m, 0);
}

TEST(EdgeListToAsymmetricGraph, DuplicateEdges) {
  // Check that `EdgeListToAsymmetricGraph()` works even when there are
  // duplicate edges in the input.
  //
  // Graph diagram:
  // 0 --> 1
  const std::vector<gi::Edge<NoWeight>> kEdges{
      {0, 1}, {0, 1},
  };
  auto graph{gi::edge_list_to_asymmetric_graph(kEdges)};
  EXPECT_EQ(graph.n, 2);
  EXPECT_EQ(graph.m, 1);

  {
    auto vertex{graph.get_vertex(0)};
    const std::vector<uintE> kExpectedOutNeighbors{1};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
    EXPECT_EQ(vertex.in_degree(), 0);
  }
  {
    auto vertex{graph.get_vertex(1)};
    const std::vector<uintE> kExpectedInNeighbors{0};
    EXPECT_EQ(vertex.out_degree(), 0);
    gt::CheckUnweightedInNeighbors(vertex, kExpectedInNeighbors);
  }
}

TEST(EdgeListToAsymmetricGraph, SkipFirstVertex) {
  // Check that `EdgeListToAsymmetricGraph()` works even when the first vertex
  // has degree 0.
  //
  // Graph diagram:
  // 1 --> 2
  const std::vector<gi::Edge<NoWeight>> kEdges{
      {1, 2},
  };
  auto graph{gi::edge_list_to_asymmetric_graph(kEdges)};
  EXPECT_EQ(graph.n, 3);
  EXPECT_EQ(graph.m, kEdges.size());

  {
    auto vertex{graph.get_vertex(0)};
    EXPECT_EQ(vertex.out_degree(), 0);
    EXPECT_EQ(vertex.in_degree(), 0);
  }
  {
    auto vertex{graph.get_vertex(1)};
    const std::vector<uintE> kExpectedOutNeighbors{2};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
    EXPECT_EQ(vertex.in_degree(), 0);
  }
  {
    auto vertex{graph.get_vertex(2)};
    const std::vector<uintE> kExpectedInNeighbors{1};
    EXPECT_EQ(vertex.out_degree(), 0);
    gt::CheckUnweightedInNeighbors(vertex, kExpectedInNeighbors);
  }
}

TEST(EdgeListToAsymmetricGraph, OutOfOrderEdges) {
  // Check that `EdgeListToAsymmetricGraph()` works even when the input edge
  // list is in a scrambled order.
  //
  // Graph diagram:
  // 0 --> 1     3 --> 6
  // | ^   ^
  // |  \  |     4
  // v   \ |
  // 2 <-> 5

  const std::vector<gi::Edge<NoWeight>> kEdges{
      {3, 6}, {0, 2}, {5, 0}, {5, 1}, {2, 5}, {0, 1}, {5, 2},
  };
  auto graph{gi::edge_list_to_asymmetric_graph(kEdges)};
  EXPECT_EQ(graph.n, 7);
  EXPECT_EQ(graph.m, kEdges.size());

  {
    auto vertex{graph.get_vertex(0)};
    const std::vector<uintE> kExpectedOutNeighbors{1, 2};
    const std::vector<uintE> kExpectedInNeighbors{5};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(vertex, kExpectedInNeighbors);
  }
  {
    auto vertex{graph.get_vertex(1)};
    const std::vector<uintE> kExpectedOutNeighbors{};
    const std::vector<uintE> kExpectedInNeighbors{0, 5};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(vertex, kExpectedInNeighbors);
  }
  {
    auto vertex{graph.get_vertex(2)};
    const std::vector<uintE> kExpectedOutNeighbors{5};
    const std::vector<uintE> kExpectedInNeighbors{0, 5};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(vertex, kExpectedInNeighbors);
  }
  {
    auto vertex{graph.get_vertex(3)};
    const std::vector<uintE> kExpectedOutNeighbors{6};
    const std::vector<uintE> kExpectedInNeighbors{};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(vertex, kExpectedInNeighbors);
  }
  {
    auto vertex{graph.get_vertex(4)};
    const std::vector<uintE> kExpectedOutNeighbors{};
    const std::vector<uintE> kExpectedInNeighbors{};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(vertex, kExpectedInNeighbors);
  }
  {
    auto vertex{graph.get_vertex(5)};
    const std::vector<uintE> kExpectedOutNeighbors{0, 1, 2};
    const std::vector<uintE> kExpectedInNeighbors{2};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(vertex, kExpectedInNeighbors);
  }
  {
    auto vertex{graph.get_vertex(6)};
    const std::vector<uintE> kExpectedOutNeighbors{};
    const std::vector<uintE> kExpectedInNeighbors{3};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
    gt::CheckUnweightedInNeighbors(vertex, kExpectedInNeighbors);
  }
}

TEST(EdgeListToSymmetricGraph, NoEdges) {
  const std::vector<gi::Edge<NoWeight>> kEdges{};
  const auto graph{gi::edge_list_to_symmetric_graph(kEdges)};
  EXPECT_EQ(graph.n, 0);
  EXPECT_EQ(graph.m, 0);
}

TEST(EdgeListToSymmetricGraph, DuplicateEdges) {
  // Check that `EdgeListToSymmetricGraph()` works even when there are
  // duplicate edges in the input.
  //
  // Graph diagram:
  // 0 --- 1
  const std::vector<gi::Edge<NoWeight>> kEdges{
      {0, 1}, {1, 0}, {0, 1},
  };
  auto graph{gi::edge_list_to_symmetric_graph(kEdges)};
  EXPECT_EQ(graph.n, 2);
  EXPECT_EQ(graph.m, 2);

  {
    auto vertex{graph.get_vertex(0)};
    const std::vector<uintE> kExpectedOutNeighbors{1};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
  }
  {
    auto vertex{graph.get_vertex(1)};
    const std::vector<uintE> kExpectedOutNeighbors{0};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
  }
}

TEST(EdgeListToSymmetricGraph, SkipFirstVertex) {
  // Check that `EdgeListToSymmetricGraph()` works even when the first vertex
  // has degree 0.
  //
  // Graph diagram:
  // 1 --- 2
  const std::vector<gi::Edge<NoWeight>> kEdges{
      {1, 2},
  };
  auto graph{gi::edge_list_to_symmetric_graph(kEdges)};
  EXPECT_EQ(graph.n, 3);
  EXPECT_EQ(graph.m, 2);

  {
    auto vertex{graph.get_vertex(0)};
    EXPECT_EQ(vertex.out_degree(), 0);
  }
  {
    auto vertex{graph.get_vertex(1)};
    const std::vector<uintE> kExpectedOutNeighbors{2};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
  }
  {
    auto vertex{graph.get_vertex(2)};
    const std::vector<uintE> kExpectedOutNeighbors{1};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
  }
}

TEST(EdgeListToSymmetricGraph, OutOfOrderEdges) {
  // Check that `EdgeListToSymmetricGraph()` works even when the input edge
  // list is in a scrambled order.
  //
  // Graph diagram:
  // 0 --- 1     3 --- 6
  // | \   |
  // |  \  |     4
  // |   \ |
  // 2     5

  const std::vector<gi::Edge<NoWeight>> kEdges{
      {1, 5}, {0, 5}, {6, 3}, {2, 0}, {1, 0},
  };
  auto graph{gi::edge_list_to_symmetric_graph(kEdges)};
  EXPECT_EQ(graph.n, 7);
  EXPECT_EQ(graph.m, 10);

  {
    auto vertex{graph.get_vertex(0)};
    const std::vector<uintE> kExpectedOutNeighbors{1, 2, 5};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
  }
  {
    auto vertex{graph.get_vertex(1)};
    const std::vector<uintE> kExpectedOutNeighbors{0, 5};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
  }
  {
    auto vertex{graph.get_vertex(2)};
    const std::vector<uintE> kExpectedOutNeighbors{0};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
  }
  {
    auto vertex{graph.get_vertex(3)};
    const std::vector<uintE> kExpectedOutNeighbors{6};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
  }
  {
    auto vertex{graph.get_vertex(4)};
    const std::vector<uintE> kExpectedOutNeighbors{};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
  }
  {
    auto vertex{graph.get_vertex(5)};
    const std::vector<uintE> kExpectedOutNeighbors{0, 1};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
  }
  {
    auto vertex{graph.get_vertex(6)};
    const std::vector<uintE> kExpectedOutNeighbors{3};
    gt::CheckUnweightedOutNeighbors(vertex, kExpectedOutNeighbors);
  }
}

}  // namespace gbbs
