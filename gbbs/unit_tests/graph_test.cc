#include "gbbs/graph.h"
#include "parlay/primitives.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace gbbs {

namespace {

// Purely for the sake of quickly testing the _ptr_ graph objects.
template <class Wgh>
static inline std::pair<symmetric_ptr_graph<symmetric_vertex, Wgh>,
                        symmetric_vertex<Wgh>*>
sym_ptr_graph_from_edges(sequence<std::tuple<uintE, uintE, Wgh>>& A, size_t n,
                         bool is_sorted = false) {
  using edge = std::tuple<uintE, uintE, Wgh>;
  auto get_u = [&](const edge& e) { return std::get<0>(e); };
  auto get_v = [&](const edge& e) { return std::get<1>(e); };
  auto get_w = [&](const edge& e) { return std::get<2>(e); };
  auto G = sym_graph_from_edges<Wgh>(A, n, get_u, get_v, get_w, is_sorted);
  using vertex = symmetric_vertex<Wgh>;
  vertex* vertices = gbbs::new_array_no_init<vertex>(G.n);
  for (size_t i = 0; i < G.n; i++) {
    vertices[i] = G.get_vertex(i);
  }
  auto GP = symmetric_ptr_graph<symmetric_vertex, Wgh>(
      G.n, G.m, vertices, std::move(G.deletion_fn));
  G.deletion_fn = [=]() {};
  return {std::move(GP), vertices};
}

// Constructs and returns the following test graph with edge weights:
//     2     3
//  0 --- 1 --- 2
static inline std::pair<gbbs::symmetric_ptr_graph<gbbs::symmetric_vertex, int>, symmetric_vertex<int>*>
ThreeNodeSymPtrGraphFromEdges() {
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 3;
  sequence<edge> edges(4);
  edges[0] = std::make_tuple(0, 1, 2);
  edges[1] = std::make_tuple(1, 0, 2);
  edges[2] = std::make_tuple(1, 2, 3);
  edges[3] = std::make_tuple(2, 1, 3);
  return sym_ptr_graph_from_edges(edges, n);
}

}  // namespace

TEST(TestSymGraphFromEdges, TestBrokenPath) {
  using edge = std::tuple<uintE, uintE, int>;
  uintE n = 11;
  uintE last_vtx_id = n - 1;
  auto edges = sequence<edge>((n - 1) * 2);
  /* Builds a path 0--1--2--....--9--10 */
  for (size_t i = 0; i < last_vtx_id; i++) {
    edges[2 * i] = std::make_tuple(i, i + 1, 1);
    edges[2 * i + 1] = std::make_tuple(i + 1, i, 1);
  }
  /* destroy 1--2 and 2--3 edge. replace with 0--10 and 1--10 */
  edges[2] = std::make_tuple(0, last_vtx_id, 1);
  edges[3] = std::make_tuple(last_vtx_id, 0, 1);
  edges[4] = std::make_tuple(1, last_vtx_id, 1);
  edges[5] = std::make_tuple(last_vtx_id, 1, 1);
  auto G = sym_graph_from_edges(edges, n, /* is_sorted = */ false);

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
  sequence<edge> edges(2);
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
  sequence<edge> edges(2);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  auto graph = sym_graph_from_edges(edges, n);

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
  auto G = graph;

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
  sequence<edge> edges(2);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  auto P = sym_ptr_graph_from_edges(edges, n);
  auto& graph = P.first;

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
  gbbs::free_array(P.second, graph.n);
}

TEST(TestSymPtrGraphFromEdges, EdgeSequenceConstructionWorks) {
  // Graph diagram:
  //     2     3
  //  0 --- 1 --- 2
  auto graph_and_vtxs = ThreeNodeSymPtrGraphFromEdges();
  auto& graph = graph_and_vtxs.first;
  auto edges = graph.edges();
  EXPECT_EQ(edges[0], std::make_tuple(0, 1, 2));
  EXPECT_EQ(edges[1], std::make_tuple(1, 0, 2));
  EXPECT_EQ(edges[2], std::make_tuple(1, 2, 3));
  EXPECT_EQ(edges[3], std::make_tuple(2, 1, 3));
  gbbs::free_array(graph_and_vtxs.second, graph.n);
}

TEST(TestSymPtrGraphFromEdges, MapEdgesWorks) {
  // Graph diagram:
  //     2     3
  //  0 --- 1 --- 2
  auto graph_and_vtxs = ThreeNodeSymPtrGraphFromEdges();
  auto& graph = graph_and_vtxs.first;
  // NOTE: graph.mapEdges() does not return anything so we are testing it by
  // disabling parallelism and sum the edge weights into a single variable for
  // testing purposes.
  int sum = 0;
  auto map_fn = [&](uintE x, uintE y, int w) { sum += w; };
  graph.mapEdges(map_fn, /*parallel_inner_map=*/false);
  EXPECT_EQ(sum, /* 2 + 2 + 3 + 3 = */ 10);
  gbbs::free_array(graph_and_vtxs.second, graph.n);
}

TEST(TestSymPtrGraphFromEdges, ReduceEdgesWorks) {
  // Graph diagram:
  //     2     3
  //  0 --- 1 --- 2
  auto graph_and_vtxs = ThreeNodeSymPtrGraphFromEdges();
  auto& graph = graph_and_vtxs.first;
  auto map_fn = [](uintE x, uintE y, int w) -> int { return w; };
  auto reduce_fn = parlay::make_monoid([](int a, int b) { return a + b; }, 0);
  auto result = graph.reduceEdges(map_fn, reduce_fn);
  EXPECT_EQ(result, /* 2 + 2 + 3 + 3 = */ 10);
  gbbs::free_array(graph_and_vtxs.second, graph.n);
}

TEST(TestSymPtrGraphFromEdges, TestBrokenPath) {
  using edge = std::tuple<uintE, uintE, int>;
  uintE n = 11;
  uintE last_vtx_id = n - 1;
  auto edges = sequence<edge>((n - 1) * 2);
  /* Builds a path 0--1--2--....--9--10 */
  for (size_t i = 0; i < last_vtx_id; i++) {
    edges[2 * i] = std::make_tuple(i, i + 1, 1);
    edges[2 * i + 1] = std::make_tuple(i + 1, i, 1);
  }
  /* destroy 1--2 and 2--3 edge. replace with 0--10 and 1--10 */
  edges[2] = std::make_tuple(0, last_vtx_id, 1);
  edges[3] = std::make_tuple(last_vtx_id, 0, 1);
  edges[4] = std::make_tuple(1, last_vtx_id, 1);
  edges[5] = std::make_tuple(last_vtx_id, 1, 1);
  auto P = sym_ptr_graph_from_edges(edges, n, /* is_sorted = */ false);
  auto& G = P.first;

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
  gbbs::free_array(P.second, G.n);
}

TEST(TestSymPtrGraphCopy, TestGraphWithSingletons) {
  // Graph diagram:
  // 0 -- 1    2    3
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 4;
  sequence<edge> edges(2);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  auto P = sym_ptr_graph_from_edges(edges, n);
  auto& graph = P.first;

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);

  auto G = graph;

  ASSERT_EQ(G.n, n);
  ASSERT_EQ(G.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(G.get_vertex(3).out_degree(), 0);
  gbbs::free_array(P.second, graph.n);
}

}  // namespace gbbs
