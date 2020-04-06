#include "ligra/graph_io.h"

#include <vector>

#include <gtest/gtest.h>

namespace g = gbbs_io;

TEST(EdgeListToAsymmetricGraph, NoEdges) {
  const std::vector<g::Edge<pbbslib::empty>> kEdges;
  const auto graph{g::edge_list_to_asymmetric_graph(kEdges)};
  EXPECT_EQ(graph.n, 0);
  EXPECT_EQ(graph.m, 0);
}

TEST(EdgeListToSymmetricGraph, NoEdges) {
  const std::vector<g::Edge<pbbslib::empty>> kEdges;
  const auto graph{g::edge_list_to_symmetric_graph(kEdges)};
  EXPECT_EQ(graph.n, 0);
  EXPECT_EQ(graph.m, 0);
}
