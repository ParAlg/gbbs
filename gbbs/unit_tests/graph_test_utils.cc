#include "gbbs/unit_tests/graph_test_utils.h"

#include <tuple>

namespace gbbs {
namespace graph_test {

symmetric_graph<symmetric_vertex, gbbs::empty> MakeUnweightedSymmetricGraph(
    const uintE num_vertices, const std::unordered_set<UndirectedEdge>& edges) {
  using Edge = std::tuple<uintE, uintE, gbbs::empty>;

  constexpr gbbs::empty weight{};
  sequence<Edge> edge_sequence(edges.size() * 2);
  auto edges_it{edges.cbegin()};
  for (size_t i = 0; i < edges.size(); i++) {
    edge_sequence[2 * i] = std::make_tuple(
        edges_it->endpoints().first, edges_it->endpoints().second, weight);
    edge_sequence[2 * i + 1] = std::make_tuple(
        edges_it->endpoints().second, edges_it->endpoints().first, weight);
    ++edges_it;
  }
  // TODO(tomtseng): Some graph operations assume that the neighbor lists are
  // sorted, so we sort here. But maybe this sorting belongs in
  // `sym_graph_from_edges` or as an option to `sym_graph_from_edges`.
  // See https://github.com/ldhulipala/gbbs/pull/21.
  parlay::sample_sort_inplace(
      make_slice(edge_sequence), [](const Edge& left, const Edge& right) {
        return std::tie(std::get<0>(left), std::get<1>(left)) <
               std::tie(std::get<0>(right), std::get<1>(right));
      });
  constexpr bool kEdgesAreSorted{true};
  return sym_graph_from_edges(edge_sequence, num_vertices, kEdgesAreSorted);
}

}  // namespace graph_test
}  // namespace gbbs
