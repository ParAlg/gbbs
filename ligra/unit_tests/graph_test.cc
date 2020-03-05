#include <gtest/gtest.h>
#include "ligra/graph.h"
#include "ligra/graph_test_utils.h"
#include "pbbslib/seq.h"

namespace gt = graph_test;

TEST(TestSymGraphFromEdges, TestBrokenPath) {
  using edge = std::tuple<uintE, uintE, int>;
  uintE n = 11;
  uintE last_vtx_id = n-1;
  auto edges = pbbs::sequence<edge>((n-1)*2);
  /* Builds a path 0--1--2--....--9--10 */
  for (size_t i=0; i<last_vtx_id; i++) {
    edges[2*i] = std::make_tuple(i, i+1, 1);
    edges[2*i + 1] = std::make_tuple(i+1, i, 1);
  }
  /* destroy 1--2 and 2--3 edge. replace with 0--10 and 1--10 */
  edges[2] = std::make_tuple(0, last_vtx_id, 1);
  edges[3] = std::make_tuple(last_vtx_id, 0, 1);
  edges[4] = std::make_tuple(1, last_vtx_id, 1);
  edges[5] = std::make_tuple(last_vtx_id, 1, 1);
  auto G = sym_graph_from_edges(edges, n, /* is_sorted = */false);

  ASSERT_EQ(G.n, 11);
  ASSERT_EQ(G.get_vertex(0).getOutDegree(), 2);
  ASSERT_EQ(G.get_vertex(1).getOutDegree(), 2);
  ASSERT_EQ(G.get_vertex(2).getOutDegree(), 0);
  ASSERT_EQ(G.get_vertex(3).getOutDegree(), 1);
  ASSERT_EQ(G.get_vertex(4).getOutDegree(), 2);
  ASSERT_EQ(G.get_vertex(5).getOutDegree(), 2);
  ASSERT_EQ(G.get_vertex(6).getOutDegree(), 2);
  ASSERT_EQ(G.get_vertex(7).getOutDegree(), 2);
  ASSERT_EQ(G.get_vertex(8).getOutDegree(), 2);
  ASSERT_EQ(G.get_vertex(9).getOutDegree(), 2);
  ASSERT_EQ(G.get_vertex(10).getOutDegree(), 3);
}

TEST(symmetric_vertex, Intersect) {
  using Vertex = symmetric_vertex<pbbs::empty>;

  // Graph diagram:
  //   0 - 1 - 2
  //    \ / \ /
  //     3 - 4 -- 5
  constexpr uintE kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{
    {0, 1},
    {0, 3},
    {1, 2},
    {1, 3},
    {1, 4},
    {2, 4},
    {3, 4},
    {4, 5},
  };
  auto graph{gt::MakeUnweightedSymmetricGraph(
      kNumVertices, kEdges, gt::ShouldSortNeighbors::kYes)};

  {
    const uintE u_id{0};
    const uintE v_id{5};
    Vertex u{graph.get_vertex(u_id)};
    Vertex v{graph.get_vertex(v_id)};
    EXPECT_EQ((u.intersect(&v, u_id, v_id)), 0);
    EXPECT_EQ((v.intersect(&u, v_id, u_id)), 0);
  }
  {
    const uintE u_id{0};
    const uintE v_id{1};
    Vertex u{graph.get_vertex(u_id)};
    Vertex v{graph.get_vertex(v_id)};
    EXPECT_EQ((u.intersect(&v, u_id, v_id)), 1);
    EXPECT_EQ((v.intersect(&u, v_id, u_id)), 1);
  }
  {
    const uintE u_id{1};
    const uintE v_id{3};
    Vertex u{graph.get_vertex(u_id)};
    Vertex v{graph.get_vertex(v_id)};
    EXPECT_EQ((u.intersect(&v, u_id, v_id)), 2);
    EXPECT_EQ((v.intersect(&u, v_id, u_id)), 2);
  }
}
