#include <gtest/gtest.h>
#include "gbbs/graph.h"
#include "pbbslib/seq.h"

namespace gbbs {

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
  ASSERT_EQ(G.get_vertex(0).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(1).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(G.get_vertex(3).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(4).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(5).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(6).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(7).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(8).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(9).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(10).out_degree(), 3);
}

TEST(TestSymGraphFromEdges, TestGraphWithSingletons) {
  // Graph diagram:
  // 0 -- 1    2    3
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 4;
  pbbs::sequence<edge> edges(2);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  auto graph = sym_graph_from_edges(edges, n);

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
}

TEST(TestSymGraphCopy, TestCopyGraphWithSingletons) {
  // Graph diagram:
  // 0 -- 1    2    3
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 4;
  pbbs::sequence<edge> edges(2);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  auto graph = sym_graph_from_edges(edges, n);

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
  auto G = graph.copy();

  ASSERT_EQ(G.n, n);
  ASSERT_EQ(G.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(G.get_vertex(3).out_degree(), 0);
}

TEST(TestSymPtrGraphFromEdges, TestGraphWithSingletons) {
  // Graph diagram:
  // 0 -- 1    2    3
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 4;
  pbbs::sequence<edge> edges(2);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  auto graph = sym_ptr_graph_from_edges(edges, n);

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
}

TEST(TestSymPtrGraphFromEdges, TestBrokenPath) {
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
  auto G = sym_ptr_graph_from_edges(edges, n, /* is_sorted = */false);

  ASSERT_EQ(G.n, 11);
  ASSERT_EQ(G.get_vertex(0).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(1).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(G.get_vertex(3).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(4).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(5).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(6).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(7).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(8).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(9).out_degree(), 2);
  ASSERT_EQ(G.get_vertex(10).out_degree(), 3);
}

TEST(TestSymPtrGraphCopy, TestGraphWithSingletons) {
  // Graph diagram:
  // 0 -- 1    2    3
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 4;
  pbbs::sequence<edge> edges(2);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  auto graph = sym_ptr_graph_from_edges(edges, n);

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);

  auto G = graph.copy();

  ASSERT_EQ(G.n, n);
  ASSERT_EQ(G.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(G.get_vertex(3).out_degree(), 0);
}

}  // namespace gbbs
