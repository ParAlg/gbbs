#pragma once

#include "gbbs/gbbs.h"

namespace gbbs {

template <class Graph, class F>
auto reorder_graph(Graph& G, F& edge_pri) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  auto offs = sequence<size_t>(n + 1);
  parallel_for(0, n, [&](size_t i) { offs[i] = G.get_vertex(i).out_degree(); });
  offs[n] = 0;
  size_t m = parlay::scan_inplace(make_slice(offs));
  assert(G.m == m);

  using edge_w = std::tuple<uintE, W>;
  auto edges = gbbs::new_array_no_init<edge_w>(m);
  parallel_for(0, n,
               [&](size_t i) {
                 size_t ctr = 0;
                 size_t off = offs[i];
                 auto map_f = [&](const uintE& u, const uintE& v,
                                  const W& wgh) {
                   edges[off + ctr] = std::make_tuple(v, wgh);
                   ctr++;
                 };
                 G.get_vertex(i).out_neighbors().map(map_f, false);

                 auto comp_f = [&](const edge_w& l, const edge_w& r) {
                   auto ep_l = edge_pri(i, l);
                   auto ep_r = edge_pri(i, r);
                   if (ep_l < ep_r) {
                     return true;
                   } else if (ep_l > ep_r) {
                     return false;
                   } else {  // lex
                     return std::get<0>(l) < std::get<0>(r);
                   }
                 };
                 auto ngh_seq = gbbs::make_slice(edges + off, ctr);
                 parlay::sample_sort_inplace(ngh_seq, comp_f);
               },
               1);

  auto v_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    size_t o = offs[i];
    v_data[i].offset = o;
    v_data[i].degree = offs[i + 1] - o;
  });

  return symmetric_graph<symmetric_vertex, W>(v_data, n, m,
                                              [=]() {
                                                gbbs::free_array(v_data, n);
                                                gbbs::free_array(edges, m);
                                              },
                                              edges);
}

}  // namespace gbbs
