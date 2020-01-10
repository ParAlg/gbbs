#pragma once

#include "ligra/ligra.h"
#include "ligra/graph.h"
#include "pbbslib/seq.h"

template <
    template <class W> class vertex, class W, typename P,
    typename std::enable_if<std::is_same<vertex<W>, asymmetric_vertex<W>>::value,
                            int>::type = 0>
inline auto relabel_graph(asymmetric_graph<vertex, W>& G, uintE* rank, P& pred) -> decltype(G) {
  std::cout << "Filter graph not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return G;
}

template <
    template <class W> class vertex, class W, typename P,
    typename std::enable_if<
        std::is_same<vertex<W>, cav_bytepd_amortized<W>>::value, int>::type = 0>
inline auto relabel_graph(asymmetric_graph<vertex, W>& G,uintE* rank, P& pred) -> decltype(G) {
  std::cout << "Filter graph not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return G;
}




template <template <class W> class vertex, class W, typename P,
          typename std::enable_if<
              std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value,
              int>::type = 0>
inline auto relabel_graph(symmetric_graph<vertex, W>& G, uintE* rank, P& pred) -> decltype(G) {
  assert(false);  // Not implemented for directed graphs
  return G;
  //TODO: basically same as filter, but we have to write out neighbors to do the sort per vert, then use this to compute
  //the byte offsets since it is dependent on rank
}


template <template <class W> class vertex, class W, typename P,
    typename std::enable_if<std::is_same<vertex<W>, symmetric_vertex<W>>::value,
                            int>::type = 0>
inline symmetric_graph<symmetric_vertex, W> relabel_graph(symmetric_graph<vertex, W>& GA, uintE* rank, P& pred) {
  using w_vertex = vertex<W>;
  auto G = filter_graph(GA, pred);
  size_t n = G.n;
  auto outOffsets = sequence<uintT>(n + 1);

  parallel_for(0, n, [&] (size_t i) {
    w_vertex u = G.get_vertex(i);
    outOffsets[rank[i]] = u.getOutDegree();
  }, 1);

  outOffsets[n] = 0;
  uintT outEdgeCount = pbbslib::scan_add_inplace(outOffsets);

  using edge = std::tuple<uintE, W>;

  auto out_edges = sequence<edge>(outEdgeCount);

  parallel_for(0, n, [&] (size_t i) {
    w_vertex u = G.get_vertex(i);
    size_t out_offset = outOffsets[rank[i]];
    uintE d = u.getOutDegree();
    edge* nghs = u.getOutNeighbors();
    edge* dir_nghs = out_edges.begin() + out_offset;
    for (size_t j=0; j < d; j++) {
      dir_nghs[j] = std::make_tuple(rank[std::get<0>(nghs[j])], std::get<1>(nghs[j]));
    }
    pbbslib::sample_sort (dir_nghs, d, [&](const edge u, const edge v) {
      return std::get<0>(u) < std::get<0>(v);
    }, true);
  }, 1);

  auto out_vdata = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    out_vdata[i].offset = outOffsets[i];
    out_vdata[i].degree = outOffsets[i+1]-outOffsets[i];
  });
  outOffsets.clear();

  auto out_edge_arr = out_edges.to_array();
  G.del();
  return symmetric_graph<symmetric_vertex, W>(
      out_vdata, n, outEdgeCount,
      get_deletion_fn(out_vdata, out_edge_arr),
      out_edge_arr);
}


