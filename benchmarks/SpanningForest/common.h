#pragma once

#define parent uintE
using edge = std::pair<uintE, uintE>;


auto parents_to_edges(pbbs::sequence<parent>& parents) -> pbbs::sequence<edge> {
  auto all_edges = pbbs::delayed_seq<edge>(parents.size(), [&] (uintE i) {
    return std::make_pair(i,parents[i]);
  });
  return pbbs::filter(all_edges, [&] (const edge& e) {
    return e.first != e.second && e.second != UINT_E_MAX;
  });
}

void root_forest(size_t n, pbbs::sequence<edge>& F) {
  pbbs::sequence<uintE> parents = pbbs::sequence<uintE>(n, UINT_E_MAX);
  parallel_for(0, F.size(), [&] (size_t i) {
    auto [u, v] = F[i];
    uintE m_u = std::min(u, v);
    uintE m_v = std::max(u, v);
    if (m_u == m_v) {
      std::cout << "self-loop in forest" << std::endl;
      exit(-1);
    }
  });
}

void check_spanning_forest(size_t n, pbbs::sequence<edge>& correct, pbbs::sequence<edge>& check) {
  // check sizes
  if (correct.size() != check.size()) {
    std::cout << "Correct forest has: " << correct.size() << " many edges, but supplied forest has: " << check.size() << " edges." << std::endl;
    exit(-1);
  }
  // root forests
}
