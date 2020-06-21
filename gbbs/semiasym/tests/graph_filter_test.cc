#include <gtest/gtest.h>

#include "gbbs/graph.h"
#include "gbbs/gbbs.h"
#include "gbbs/graph_test_utils.h"
#include "gbbs/semiasym/graph_filter.h"
#include "pbbslib/seq.h"

namespace gbbs {

symmetric_graph<symmetric_vertex, pbbs::empty> CreateStar(size_t n) {
  using edge = std::tuple<uintE, uintE, pbbs::empty>;
  auto edges = pbbs::sequence<edge>(2*(n-1));
  for (size_t i=0; i<(n-1); i++) {
    edges[2*i] = {0, i+1, pbbs::empty()};
    edges[2*i+1] = {i+1, 0, pbbs::empty()};
    std::cout << "(" << 0 << "," << (i+1) << ")" << std::endl;
  }
  auto G = sym_graph_from_edges(edges, n, /* is_sorted = */false);
  return G;
}

TEST(TestGraphFilter, TestCreation) {
  uintE n = 4;
  // Build a star graph
  auto G = CreateStar(n);
  auto PG = sage::build_packed_graph(G);
  EXPECT_EQ(PG.n, G.n);
  EXPECT_EQ(PG.m, G.m);
  for (size_t i=1; i < n; i++) {
    EXPECT_EQ(G.get_vertex(i).getOutDegree(), 1);
    EXPECT_EQ(PG.get_vertex(i).getOutDegree(), 1);
  }
  EXPECT_EQ(G.get_vertex(0).getOutDegree(), n-1);
  EXPECT_EQ(PG.get_vertex(0).getOutDegree(), n-1);
  std::cout << "m = " << PG.m << std::endl;
}

TEST(TestGraphFilter, TestFilter) {
  using W = pbbs::empty;
  uintE n = 4;
  auto G = CreateStar(n);
  auto predicate = [&] (const uintE& u, const uintE& v, const W& wgh) -> bool {
    return (u % 2 == v % 2) && (u % 2 == 0);
  };
  auto PG = sage::filter_graph(G, predicate);
  for (size_t i=1; i < n; i++) {
    if (i % 2 == 1) {
      EXPECT_EQ(PG.get_vertex(i).getOutDegree(), 0);
    }
  }
  EXPECT_EQ(PG.get_vertex(0).getOutDegree(), n/2 - 1);
}

//TEST(TestGraphFilter, TestEdgeMapFilter) {
//  using W = pbbs::empty;
//  using edge = std::tuple<uintE, uintE, W>;
//  uintE n = 4;
//  auto G = CreateStar(n);
//  auto PG = sage::build_packed_graph(G);
//  auto predicate = [&] (const uintE& u, const uintE& v, const W& wgh) -> bool {
//    return (u % 2 == v % 2) && (u % 2 == 0);
//  };
//  auto vset = vertexSubset(0);
//  auto packed_vset = edgeMapFilter(PG, vset, predicate, pack_edges);
//  EXPECT_EQ(PG.get_vertex(0).getOutDegree(), n/2 - 1);
//}

}  // namespace gbbs
