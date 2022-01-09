#pragma once

#include "gbbs/graph.h"

namespace gbbs {
using edge = std::pair<uintE, uintE>;
using parent = uintE;

constexpr edge empty_edge = std::make_pair(UINT_E_MAX, UINT_E_MAX);

namespace spanning_forest {
constexpr uintE largest_comp = UINT_E_MAX;

inline sequence<edge> parents_to_edges(sequence<parent>& parents) {
  auto all_edges = parlay::delayed_seq<edge>(
      parents.size(), [&](uintE i) { return std::make_pair(i, parents[i]); });
  return parlay::filter(all_edges, [&](const edge& e) {
    return (e.first != e.second) && (e.second != UINT_E_MAX);
  });
}
}
}  // namespace gbbs
