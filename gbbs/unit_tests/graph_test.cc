#include <gtest/gtest.h>
#include "gbbs/graph.h"
#include "pbbslib/seq.h"

namespace gbbs {

namespace {

  // Purely for the sake of quickly testing the _ptr_ graph objects.
  template <class Wgh>
  static inline symmetric_ptr_graph<symmetric_vertex, Wgh> sym_ptr_graph_from_edges(
      pbbs::sequence<std::tuple<uintE, uintE, Wgh>>& A, size_t n,
      bool is_sorted = false) {
    using edge = std::tuple<uintE, uintE, Wgh>;
    auto get_u = [&](const edge& e) { return std::get<0>(e); };
    auto get_v = [&](const edge& e) { return std::get<1>(e); };
    auto get_w = [&](const edge& e) { return std::get<2>(e); };
    auto G = sym_graph_from_edges<Wgh>(A, n, get_u, get_v, get_w, is_sorted);
    using vertex = symmetric_vertex<Wgh>;
    vertex* vertices = pbbs::new_array_no_init<vertex>(G.n);
    for (size_t i=0; i<G.n; i++) {
      vertices[i] = G.get_vertex(i);
    }
    auto GP = symmetric_ptr_graph<symmetric_vertex, Wgh>(G.n, G.m, vertices, G.deletion_fn);
    return GP;
  }

}

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
  G.del();
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
  ASSERT_EQ(graph.get_vertex(0).getOutDegree(), 1);
  ASSERT_EQ(graph.get_vertex(1).getOutDegree(), 1);
  ASSERT_EQ(graph.get_vertex(2).getOutDegree(), 0);
  ASSERT_EQ(graph.get_vertex(3).getOutDegree(), 0);
  graph.del();
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
  ASSERT_EQ(graph.get_vertex(0).getOutDegree(), 1);
  ASSERT_EQ(graph.get_vertex(1).getOutDegree(), 1);
  ASSERT_EQ(graph.get_vertex(2).getOutDegree(), 0);
  ASSERT_EQ(graph.get_vertex(3).getOutDegree(), 0);
  graph.del();
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
  G.del();
}

}  // namespace gbbs
