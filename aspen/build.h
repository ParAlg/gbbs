#pragma once

#include "utils.h"

namespace aspen {
namespace build {

  // Builds an aspen graph from the given static graph.
  template <class Graph>
  auto graph_to_edges(Graph& G, size_t num_batches = 1) {
    using weight = typename Graph::weight_type;
    size_t n = G.n;

    auto degs = parlay::sequence<size_t>::from_function(
        n, [&](size_t i) { return G.get_vertex(i).out_degree(); });
    size_t sum_degs = parlay::scan_inplace(parlay::make_slice(degs));
    assert(sum_degs == G.m);

    using ngh_and_weight = std::tuple<vertex_id, weight>;
    using edge = std::pair<vertex_id, ngh_and_weight>;

    auto edges = parlay::sequence<edge>::uninitialized(sum_degs);

    parlay::parallel_for(0, n, [&](size_t i) {
      size_t k = degs[i];
      auto map_f = [&](const vertex_id& u, const vertex_id& v, const weight& wgh) {
       edges[k++] = std::make_pair(u, std::make_tuple(v, wgh));
      };
      G.get_vertex(i).out_neighbors().map(map_f, false);
    }, 1);
    return edges;
  }

}  // namespace build
}  // namespace aspen
