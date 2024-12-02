#include "gbbs/unit_tests/graph_test_utils.h"

#include <tuple>

namespace gbbs {
namespace graph_test {

symmetric_graph<symmetric_vertex, gbbs::empty> MakeUnweightedSymmetricGraph(
    const uintE num_vertices, const std::unordered_set<UndirectedEdge>& edges) {
  using Edge = std::tuple<uintE, uintE, gbbs::empty>;

  constexpr gbbs::empty weight{};
  sequence<Edge> edge_sequence(edges.size());
  auto edges_it{edges.cbegin()};
  for (size_t i = 0; i < edges.size(); i++) {
    edge_sequence[i] = std::make_tuple(
        edges_it->endpoints().first, edges_it->endpoints().second, weight);
    ++edges_it;
  }
  return symmetric_graph<symmetric_vertex, gbbs::empty>::from_edges(edge_sequence, num_vertices);
}

symmetric_ptr_graph<symmetric_vertex, gbbs::empty> MakeUnweightedSymmetricPtrGraph(
    const uintE num_vertices, const std::unordered_set<UndirectedEdge>& edges) {
  using Edge = std::tuple<uintE, uintE, gbbs::empty>;
  using Weight = gbbs::empty;

  constexpr gbbs::empty weight{};
  sequence<Edge> edge_sequence(edges.size());
  auto edges_it{edges.cbegin()};
  for (size_t i = 0; i < edges.size(); i++) {
    edge_sequence[i] = std::make_tuple(
        edges_it->endpoints().first, edges_it->endpoints().second, weight);
    ++edges_it;
  }

  auto G = symmetric_graph<symmetric_vertex, gbbs::empty>::from_edges(edge_sequence, num_vertices);
  size_t n = G.n;

  using vertex = symmetric_vertex<Weight>;
  vertex* vertices = gbbs::new_array_no_init<vertex>(n);
  parlay::parallel_for(0, n, [&] (size_t i) {
    vertices[i] = G.get_vertex(i);
  });

  auto deletion_fn = G.deletion_fn;

  auto GP = symmetric_ptr_graph<symmetric_vertex, Weight>(
      G.n, G.m, vertices, [deletion_fn, vertices, n] () {
        deletion_fn();
        gbbs::free_array(vertices, n);
      });
  G.deletion_fn = [=]() {};
  return GP;
}

}  // namespace graph_test
}  // namespace gbbs
