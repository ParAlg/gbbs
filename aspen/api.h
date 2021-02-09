#pragma once

#include "immutable_graph.h"
#include "traversable_graph.h"

namespace aspen {

template <class Graph>
auto symmetric_graph_from_static_graph(Graph& GA) {
  using W = typename Graph::weight_type;
  using inner_graph = symmetric_graph<W>;
  using outer_graph = traversable_graph<inner_graph>;
  return outer_graph(GA);
}

}  // namespace aspen
