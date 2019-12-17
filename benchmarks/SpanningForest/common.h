#pragma once

#include "ligra/graph.h"

using edge = std::pair<uintE, uintE>;
using parent = uintE;

namespace spanning_forest {

  auto parents_to_edges(pbbs::sequence<parent>& parents) -> pbbs::sequence<edge> {
    auto all_edges = pbbs::delayed_seq<edge>(parents.size(), [&] (uintE i) {
      return std::make_pair(i,parents[i]);
    });
    return pbbs::filter(all_edges, [&] (const edge& e) {
      return (e.first != e.second) && (e.second != UINT_E_MAX);
    });
  }
}
