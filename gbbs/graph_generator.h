#pragma once

#include <cstddef>
#include <tuple>

#include "bridge.h"
#include "graph.h"
#include "macros.h"
#include "gbbs/vertex.h"

namespace gbbs {
namespace graph_generator {

// Build edges for clique with nodes 0, ..., n - 1.
inline sequence<std::tuple<uintE, uintE, empty>> clique_edges(uintE n) {
  if (n == 0) {
    return sequence<std::tuple<uintE, uintE, empty>>();
  }
  size_t m = n * (n - 1);
  auto edges = sequence<std::tuple<uintE, uintE, empty>>(m);
  parallel_for(0, n, [&](size_t i) {
    size_t offset = i * (n - 1);
    size_t k = 0;
    for (size_t j = 0; j < n; ++j) {
      if (i != j) {
        edges[offset + k++] = std::make_tuple(i, j, empty{});
      }
    }
  });
  return edges;
}

// Build edges for path graph 0-1-2-...-(n-1)
inline sequence<std::tuple<uintE, uintE, empty>> path_edges(uintE n) {
  if (n <= 1) {
    return sequence<std::tuple<uintE, uintE, empty>>();
  }
  size_t m = 2 * (n - 1);
  auto edges = sequence<std::tuple<uintE, uintE, empty>>(m);
  parallel_for(0, n - 1, [&](size_t i) {
    edges[2 * i] = std::make_tuple(i, i + 1, empty{});
    edges[2 * i + 1] = std::make_tuple(i + 1, i, empty{});
  });
  return edges;
}

// Build edges for cycle graph 0-1-2-...-(n-1)-0
inline sequence<std::tuple<uintE, uintE, empty>> cycle_edges(uintE n) {
  if (n <= 1) {
    return sequence<std::tuple<uintE, uintE, empty>>();
  }
  size_t m = 2 * n;
  auto edges = sequence<std::tuple<uintE, uintE, empty>>(m);
  parallel_for(0, n, [&](size_t i) {
    edges[2 * i] = std::make_tuple(i, (i + 1) % n, empty{});
    edges[2 * i + 1] = std::make_tuple((i + 1) % n, i, empty{});
  });
  return edges;
}

inline symmetric_graph<symmetric_vertex, empty> clique_graph(uintE n) {
  return symmetric_graph<symmetric_vertex, empty>::from_edges(clique_edges(n),
                                                              n);
}

inline symmetric_graph<symmetric_vertex, empty> path_graph(uintE n) {
  return symmetric_graph<symmetric_vertex, empty>::from_edges(path_edges(n), n);
}

inline symmetric_graph<symmetric_vertex, empty> cycle_graph(uintE n) {
  return symmetric_graph<symmetric_vertex, empty>::from_edges(cycle_edges(n),
                                                              n);
}

}  // namespace graph_generator
}  // namespace gbbs
