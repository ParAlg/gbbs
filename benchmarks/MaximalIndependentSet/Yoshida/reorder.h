#pragma once

#include "gbbs/gbbs.h"

namespace gbbs {

template <class Graph, class F>
auto reorder_graph(Graph& G, F& vertex_pri) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  auto offs = pbbs::sequence<size_t>(n+1);
  parallel_for(0, n, [&] (size_t i) {
    offs[i] = G.get_vertex(i).getOutDegree();
  });
  offs[n] = 0;
  size_t m = pbbslib::scan_add_inplace(offs.slice());
  assert(G.m == m);


  using edge_w = std::tuple<uintE, W>;
  auto edges = pbbs::new_array_no_init<edge_w>(m);
  parallel_for(0, n, [&] (size_t i) {
    size_t ctr = 0;
    size_t off = offs[i];
    auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      edges[off + ctr] = std::make_tuple(v, wgh);
      ctr++;
    };
    G.get_vertex(i).mapOutNgh(i, map_f, false);

    auto comp_f = [&] (const edge_w& l, const edge_w& r) {
      auto pri_l = vertex_pri(std::get<0>(l));
      auto pri_r = vertex_pri(std::get<0>(r));
      if (pri_l < pri_r) {
        return true;
      } else if (pri_l > pri_r) {
        return false;
      } else { // lex
        return std::get<0>(l) < std::get<0>(r);
      }
    };
    auto ngh_seq = pbbslib::make_sequence(edges + off, ctr);
    pbbs::sample_sort_inplace(ngh_seq, comp_f);
  }, 1);

  auto v_data = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    size_t o = offs[i];
    v_data[i].offset = o;
    v_data[i].degree = offs[i+1] - o;
  });

  return symmetric_graph<symmetric_vertex, W>(v_data, n, m, [=]() {pbbslib::free_arrays(v_data, edges);}, edges);
}

}  // namespace gbbs
