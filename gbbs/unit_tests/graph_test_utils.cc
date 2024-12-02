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

  return symmetric_ptr_graph<symmetric_vertex, gbbs::empty>::from_edges(edge_sequence, num_vertices);
}

asymmetric_graph<asymmetric_vertex, gbbs::empty> MakeUnweightedAsymmetricGraph(
    const uintE num_vertices, const std::unordered_set<DirectedEdge>& edges) {
  using Edge = std::tuple<uintE, uintE, gbbs::empty>;

  constexpr gbbs::empty weight{};
  sequence<Edge> edge_sequence(edges.size());
  auto edges_it{edges.cbegin()};
  for (size_t i = 0; i < edges.size(); i++) {
    edge_sequence[i] = std::make_tuple(
        edges_it->endpoints().first, edges_it->endpoints().second, weight);
    ++edges_it;
  }
  return asymmetric_graph<asymmetric_vertex, gbbs::empty>::from_edges(edge_sequence, num_vertices);
}

}  // namespace graph_test
}  // namespace gbbs
