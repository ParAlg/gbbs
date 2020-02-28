#include "ligra/graph_test_utils.h"

#include <tuple>

namespace graph_test {

symmetric_graph<symmetric_vertex, pbbslib::empty> MakeUnweightedSymmetricGraph(
    const uintE num_vertices,
    const std::unordered_set<UndirectedEdge>& edges) {
  constexpr pbbs::empty weight{};
  pbbs::sequence<std::tuple<uintE, uintE, pbbslib::empty>> edge_sequence(
      edges.size() * 2);
  auto edges_it{edges.cbegin()};
  for (size_t i = 0; i < edges.size(); i++) {
    edge_sequence[2 * i] =
      std::make_tuple(
          edges_it->endpoints().first,
          edges_it->endpoints().second,
          weight);
    edge_sequence[2 * i + 1] =
      std::make_tuple(
          edges_it->endpoints().second,
          edges_it->endpoints().first,
          weight);
    ++edges_it;
  }
  return sym_graph_from_edges(edge_sequence, num_vertices);
}

}  // namespace graph_test
