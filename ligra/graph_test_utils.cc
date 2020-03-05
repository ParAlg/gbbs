#include "ligra/graph_test_utils.h"

#include <tuple>

#include "pbbslib/sample_sort.h"

namespace graph_test {

symmetric_graph<symmetric_vertex, pbbslib::empty> MakeUnweightedSymmetricGraph(
    const uintE num_vertices,
    const std::unordered_set<UndirectedEdge>& edges,
    const ShouldSortNeighbors should_sort_neighbors) {
  using Edge = std::tuple<uintE, uintE, pbbslib::empty>;
  constexpr pbbs::empty weight{};

  pbbs::sequence<Edge> edge_sequence(edges.size() * 2);
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

  switch (should_sort_neighbors) {
    case ShouldSortNeighbors::kYes: {
      pbbs::sample_sort_inplace(
          edge_sequence.slice(),
          [](const Edge& left, const Edge& right) {
            return std::tie(std::get<0>(left), std::get<1>(left))
              < std::tie(std::get<0>(right), std::get<1>(right));
          });
      const bool edges_are_sorted{true};
      return sym_graph_from_edges(edge_sequence, num_vertices, edges_are_sorted);
    }
    case ShouldSortNeighbors::kNo: {
      return sym_graph_from_edges(edge_sequence, num_vertices);
    }
  }
}

}  // namespace graph_test
