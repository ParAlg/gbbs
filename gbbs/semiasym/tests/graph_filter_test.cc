#include <gtest/gtest.h>

#include "gbbs/gbbs.h"
#include "gbbs/graph.h"
#include "gbbs/semiasym/graph_filter.h"
#include "gbbs/unit_tests/graph_test_utils.h"

namespace gbbs {

symmetric_graph<symmetric_vertex, gbbs::empty> CreateStar(size_t n) {
  using edge = std::tuple<uintE, uintE, gbbs::empty>;
  auto edges = sequence<edge>(2 * (n - 1));
  for (size_t i = 0; i < (n - 1); i++) {
    edges[2 * i] = {0, i + 1, gbbs::empty()};
    edges[2 * i + 1] = {i + 1, 0, gbbs::empty()};
    std::cout << "(" << 0 << "," << (i + 1) << ")" << std::endl;
  }
  auto G = sym_graph_from_edges(edges, n, /* is_sorted = */ false);
  return G;
}

TEST(TestGraphFilter, TestCreation) {
  uintE n = 4;
  // Build a star graph
  auto G = CreateStar(n);
  auto PG = sage::build_symmetric_packed_graph(G);
  EXPECT_EQ(PG.n, G.n);
  EXPECT_EQ(PG.m, G.m);
  for (size_t i = 1; i < n; i++) {
    EXPECT_EQ(G.get_vertex(i).out_degree(), 1);
    EXPECT_EQ(PG.get_vertex(i).out_degree(), 1);
  }
  EXPECT_EQ(G.get_vertex(0).out_degree(), n - 1);
  EXPECT_EQ(PG.get_vertex(0).out_degree(), n - 1);
  std::cout << "m = " << PG.m << std::endl;
}

TEST(TestGraphFilter, TestFilter) {
  using W = gbbs::empty;
  uintE n = 20;
  auto G = CreateStar(n);
  auto predicate = [&](const uintE& u, const uintE& v, const W& wgh) -> bool {
    return (u % 2 == v % 2) && (u % 2 == 0);
  };
  auto PG = sage::filter_graph(G, predicate);
  for (size_t i = 1; i < n; i++) {
    if (i % 2 == 1) {
      EXPECT_EQ(PG.get_vertex(i).out_degree(), 0);
    } else {
      EXPECT_EQ(PG.get_vertex(i).out_degree(), 1);
    }
  }
  EXPECT_EQ(PG.get_vertex(0).out_degree(), n / 2 - 1);
}

asymmetric_graph<asymmetric_vertex, gbbs::empty> CreateDirectedStar(size_t n) {
  using edge = std::tuple<uintE, uintE, gbbs::empty>;
  auto edges = sequence<edge>(n - 1);
  for (size_t i = 0; i < (n - 1); i++) {
    edges[i] = {0, i + 1,
                gbbs::empty()};  // directed arc from center -> satellite
    std::cout << "(" << 0 << "," << (i + 1) << ")" << std::endl;
  }
  auto G = asym_graph_from_edges(edges, n, /* is_sorted = */ false);
  return G;
}

TEST(TestGraphFilter, TestAsymmetricCreation) {
  uintE n = 10;
  // Build a star graph
  auto G = CreateDirectedStar(n);
  auto PG = sage::build_asymmetric_packed_graph(G);
  EXPECT_EQ(PG.n, G.n);
  EXPECT_EQ(PG.m, G.m);
  for (size_t i = 1; i < n; i++) {
    EXPECT_EQ(G.get_vertex(i).out_degree(), 0);
    EXPECT_EQ(PG.get_vertex(i).out_degree(), 0);

    EXPECT_EQ(G.get_vertex(i).in_degree(), 1);
    EXPECT_EQ(PG.get_vertex(i).in_degree(), 1);
  }
  EXPECT_EQ(G.get_vertex(0).out_degree(), n - 1);
  EXPECT_EQ(PG.get_vertex(0).out_degree(), n - 1);

  EXPECT_EQ(G.get_vertex(0).in_degree(), 0);
  EXPECT_EQ(PG.get_vertex(0).in_degree(), 0);
  std::cout << "m = " << PG.m << std::endl;
}

TEST(TestGraphFilter, TestAsymmetricFilter) {
  using W = gbbs::empty;
  uintE n = 20;
  auto G = CreateDirectedStar(n);
  auto predicate = [&](const uintE& u, const uintE& v, const W& wgh) -> bool {
    return (u % 2 == v % 2) && (u % 2 == 0);
  };
  auto PG = sage::filter_graph(G, predicate);
  for (size_t i = 1; i < n; i++) {
    EXPECT_EQ(G.get_vertex(i).out_degree(), 0);
    EXPECT_EQ(PG.get_vertex(i).out_degree(), 0);
    if (i % 2 == 1) {
      EXPECT_EQ(G.get_vertex(i).in_degree(), 1);
      EXPECT_EQ(PG.get_vertex(i).in_degree(), 0);
    } else {
      EXPECT_EQ(G.get_vertex(i).in_degree(), 1);
      EXPECT_EQ(PG.get_vertex(i).in_degree(), 1);
    }
  }
  EXPECT_EQ(PG.get_vertex(0).out_degree(), n / 2 - 1);
}

}  // namespace gbbs
