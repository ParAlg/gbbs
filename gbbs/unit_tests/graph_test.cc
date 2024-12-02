#include "gbbs/graph.h"

#include <atomic>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "parlay/primitives.h"

namespace gbbs {

namespace {

// Constructs and returns the following test graph with edge weights:
//     2     3
//  0 --- 1 --- 2
static inline gbbs::symmetric_ptr_graph<gbbs::symmetric_vertex, int>
ThreeNodeSymPtrGraphFromEdges() {
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 3;
  sequence<edge> edges(4);
  edges[0] = std::make_tuple(0, 1, 2);
  edges[1] = std::make_tuple(1, 0, 2);
  edges[2] = std::make_tuple(1, 2, 3);
  edges[3] = std::make_tuple(2, 1, 3);
  return gbbs::symmetric_ptr_graph<gbbs::symmetric_vertex, int>::from_edges(
      edges, n);
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
  auto G = symmetric_graph<symmetric_vertex, int>::from_edges(edges, n);

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
  auto graph = symmetric_graph<symmetric_vertex, int>::from_edges(edges, n);

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
  auto graph = symmetric_graph<symmetric_vertex, int>::from_edges(edges, n);

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

TEST(TestSymGraphMove, TestMoveGraphWithSingletons) {
  // Graph diagram:
  // 0 -- 1    2    3
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 4;
  sequence<edge> edges(2);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  symmetric_graph<symmetric_vertex, int> G;
  {
    auto graph = symmetric_graph<symmetric_vertex, int>::from_edges(edges, n);

    ASSERT_EQ(graph.n, n);
    ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
    ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
    ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
    ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
    G = std::move(graph);
  }

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
  auto graph = symmetric_ptr_graph<symmetric_vertex, int>::from_edges(edges, n);

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
}

TEST(TestSymPtrGraphFromEdges, EdgeSequenceConstructionWorks) {
  // Graph diagram:
  //     2     3
  //  0 --- 1 --- 2
  auto graph = ThreeNodeSymPtrGraphFromEdges();
  auto edges = graph.edges();
  EXPECT_EQ(edges[0], std::make_tuple(0, 1, 2));
  EXPECT_EQ(edges[1], std::make_tuple(1, 0, 2));
  EXPECT_EQ(edges[2], std::make_tuple(1, 2, 3));
  EXPECT_EQ(edges[3], std::make_tuple(2, 1, 3));
}

TEST(TestSymPtrGraphFromEdges, MapEdgesWorks) {
  // Graph diagram:
  //     2     3
  //  0 --- 1 --- 2
  auto graph = ThreeNodeSymPtrGraphFromEdges();
  // NOTE: graph.mapEdges() does not return anything so we are testing it by
  // disabling parallelism and sum the edge weights into a single variable for
  // testing purposes.
  std::atomic<int> sum = 0;
  auto map_fn = [&](uintE x, uintE y, int w) { sum += w; };
  graph.mapEdges(map_fn);
  EXPECT_EQ(sum, /* 2 + 2 + 3 + 3 = */ 10);
}

TEST(TestSymPtrGraphFromEdges, ReduceEdgesWorks) {
  // Graph diagram:
  //     2     3
  //  0 --- 1 --- 2
  auto graph = ThreeNodeSymPtrGraphFromEdges();
  auto map_fn = [](uintE x, uintE y, int w) -> int { return w; };
  auto reduce_fn = parlay::make_monoid([](int a, int b) { return a + b; }, 0);
  auto result = graph.reduceEdges(map_fn, reduce_fn);
  EXPECT_EQ(result, /* 2 + 2 + 3 + 3 = */ 10);
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
  auto G = symmetric_ptr_graph<symmetric_vertex, int>::from_edges(edges, n);

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
  sequence<edge> edges(2);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  auto graph = symmetric_ptr_graph<symmetric_vertex, int>::from_edges(edges, n);

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

TEST(TestSymPtrGraphMove, TestGraphWithSingletons) {
  // Graph diagram:
  // 0 -- 1    2    3
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 4;
  sequence<edge> edges(2);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  symmetric_ptr_graph<symmetric_vertex, int> G;
  auto graph = symmetric_ptr_graph<symmetric_vertex, int>::from_edges(edges, n);

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
  G = std::move(graph);

  ASSERT_EQ(G.n, n);
  ASSERT_EQ(G.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(G.get_vertex(3).out_degree(), 0);
}

TEST(TestAsymGraphFromEdges, TestPath) {
  using edge = std::tuple<uintE, uintE, int>;
  uintE n = 11;
  uintE last_vtx_id = n - 1;
  auto edges = sequence<edge>((n - 1));
  /* Builds a path 0->1->2->....->9->10 */
  for (size_t i = 0; i < last_vtx_id; i++) {
    edges[i] = std::make_tuple(i, i + 1, 1);
  }
  auto G = asymmetric_graph<asymmetric_vertex, int>::from_edges(edges, n);

  ASSERT_EQ(G.n, 11);
  ASSERT_EQ(G.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(3).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(4).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(5).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(6).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(7).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(8).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(9).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(10).out_degree(), 0);

  ASSERT_EQ(G.get_vertex(0).in_degree(), 0);
  ASSERT_EQ(G.get_vertex(1).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(3).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(4).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(5).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(6).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(7).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(8).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(9).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(10).in_degree(), 1);
}

TEST(TestAsymGraphFromEdges, TestGraphWithSingletons) {
  // Graph diagram:
  // 0 <--> 1    2    3 <-- 4
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 5;
  sequence<edge> edges(3);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  edges[2] = std::make_tuple(4, 3, 1);
  auto graph = asymmetric_graph<asymmetric_vertex, int>::from_edges(edges, n);

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(4).out_degree(), 1);

  ASSERT_EQ(graph.get_vertex(0).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).in_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(4).in_degree(), 0);
}

TEST(TestAsymGraphCopy, TestCopyGraphWithSingletons) {
  // Graph diagram:
  // 0 <--> 1    2    3 <-- 4
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 5;
  sequence<edge> edges(3);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  edges[2] = std::make_tuple(4, 3, 1);
  auto graph = asymmetric_graph<asymmetric_vertex, int>::from_edges(edges, n);

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(4).out_degree(), 1);

  ASSERT_EQ(graph.get_vertex(0).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).in_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(4).in_degree(), 0);

  // Copy
  auto G = graph;
  ASSERT_EQ(G.n, n);
  ASSERT_EQ(G.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(G.get_vertex(3).out_degree(), 0);
  ASSERT_EQ(G.get_vertex(4).out_degree(), 1);

  ASSERT_EQ(G.get_vertex(0).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(1).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).in_degree(), 0);
  ASSERT_EQ(G.get_vertex(3).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(4).in_degree(), 0);
}

TEST(TestAsymPtrGraphFromEdges, TestPath) {
  using edge = std::tuple<uintE, uintE, int>;
  uintE n = 11;
  uintE last_vtx_id = n - 1;
  auto edges = sequence<edge>((n - 1));
  /* Builds a path 0->1->2->....->9->10 */
  for (size_t i = 0; i < last_vtx_id; i++) {
    edges[i] = std::make_tuple(i, i + 1, 1);
  }
  auto G = asymmetric_ptr_graph<asymmetric_vertex, int>::from_edges(edges, n);

  ASSERT_EQ(G.n, 11);
  ASSERT_EQ(G.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(3).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(4).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(5).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(6).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(7).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(8).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(9).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(10).out_degree(), 0);

  ASSERT_EQ(G.get_vertex(0).in_degree(), 0);
  ASSERT_EQ(G.get_vertex(1).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(3).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(4).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(5).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(6).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(7).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(8).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(9).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(10).in_degree(), 1);
}

TEST(TestAsymPtrGraphFromEdges, TestGraphWithSingletons) {
  // Graph diagram:
  // 0 <--> 1    2    3 <-- 4
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 5;
  sequence<edge> edges(3);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  edges[2] = std::make_tuple(4, 3, 1);
  auto graph =
      asymmetric_ptr_graph<asymmetric_vertex, int>::from_edges(edges, n);

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(4).out_degree(), 1);

  ASSERT_EQ(graph.get_vertex(0).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).in_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(4).in_degree(), 0);
}

TEST(TestAsymPtrGraphCopy, TestCopyGraphWithSingletons) {
  // Graph diagram:
  // 0 <--> 1    2    3 <-- 4
  using edge = std::tuple<uintE, uintE, int>;
  const uintE n = 5;
  sequence<edge> edges(3);
  edges[0] = std::make_tuple(0, 1, 1);
  edges[1] = std::make_tuple(1, 0, 1);
  edges[2] = std::make_tuple(4, 3, 1);
  auto graph =
      asymmetric_ptr_graph<asymmetric_vertex, int>::from_edges(edges, n);

  ASSERT_EQ(graph.n, n);
  ASSERT_EQ(graph.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).out_degree(), 0);
  ASSERT_EQ(graph.get_vertex(4).out_degree(), 1);

  ASSERT_EQ(graph.get_vertex(0).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(1).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(2).in_degree(), 0);
  ASSERT_EQ(graph.get_vertex(3).in_degree(), 1);
  ASSERT_EQ(graph.get_vertex(4).in_degree(), 0);

  // Copy
  auto G = graph;
  ASSERT_EQ(G.n, n);
  ASSERT_EQ(G.get_vertex(0).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(1).out_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).out_degree(), 0);
  ASSERT_EQ(G.get_vertex(3).out_degree(), 0);
  ASSERT_EQ(G.get_vertex(4).out_degree(), 1);

  ASSERT_EQ(G.get_vertex(0).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(1).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(2).in_degree(), 0);
  ASSERT_EQ(G.get_vertex(3).in_degree(), 1);
  ASSERT_EQ(G.get_vertex(4).in_degree(), 0);
}

}  // namespace gbbs
