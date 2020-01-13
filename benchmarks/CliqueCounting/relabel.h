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
inline auto relabel_graph(symmetric_graph<vertex, W>& GA, uintE* rank, P& pred) -> decltype(G) {
  //assert(false);  // Not implemented for directed graphs
  //return G;
  //TODO: basically same as filter, but we have to write out neighbors to do the sort per vert, then use this to compute
  //the byte offsets since it is dependent on rank
  using edge = std::tuple<uintE, W>;
  auto G = filter_graph(GA, pred);
  size_t n = G.n;

  // 1. Calculate total size
  auto byte_offsets = sequence<uintT>(n + 1);
  auto degrees = sequence<uintE>(n);
  parallel_for(0, n, [&] (size_t i) {
    size_t total_bytes = 0;
    uintE last_ngh = 0;
    size_t deg = G.get_vertex(i).getOutDegree();
    uchar tmp[16];
    // here write out all of G's outneighbors to an uncompressed array, and then relabel and sort
    auto tmp_edges = sequence<edge>(deg);
    size_t j = 0;
    auto f = [&](uintE u, uintE v, W w) {
      tmp_edges[j] = std::make_tuple(rank[v], w);
      j++;
      return false;
    };
    G.get_vertex(i).mapOutNgh(i, f, false);
    // need to sort tmp_edges
    pbbslib::sample_sort (tmp_edges, deg, [&](const edge u, const edge v) {
      return std::get<0>(u) < std::get<0>(v);
    }, true);
    // now need to compute total_bytes (TODO: this could be parallelized, b/c we know what last_ngh is)
    for (j=0; j < deg; j++) {
      size_t bytes = 0;
      auto v = std::get<0>(tmp_edges[j]);
      auto wgh = std::get<1>(tmp_edges[j]);
      if (j == 0) {
        bytes = byte::compressFirstEdge(tmp, bytes, rank[i], v);
        bytes = byte::compressWeight<W>(tmp, bytes, wgh);
      } else {
        bytes = byte::compressEdge(tmp, bytes, v - last_ngh);
        bytes = byte::compressWeight<W>(tmp, bytes, wgh);
      }
      last_ngh = v;
      total_bytes += bytes;
    }
    tmp_edges.clear();

    byte_offsets[rank[i]] = total_bytes;
    degrees[rank[i]] = deg;
  }, 1);
  byte_offsets[n] = 0;
  size_t last_offset = pbbslib::scan_add_inplace(byte_offsets);
  //std::cout << "# size is: " << last_offset << "\n";

  auto edges = sequence<uchar>(last_offset);
  // redo the sort from above, but now actually store your edges
  // WIP
  parallel_for(0, n, [&] (size_t i) {
    uintE deg = G.get_vertex(i).getOutDegree();
    if (deg > 0) {
      // here write out all of G's outneighbors to an uncompressed array, and then relabel and sort
      auto tmp_edges = sequence<edge>(deg);
      size_t j = 0;
      auto f = [&](uintE u, uintE v, W w) {
        tmp_edges[j] = std::make_tuple(rank[v], w);
        j++;
        return false;
      };
      G.get_vertex(i).mapOutNgh(i, f, false);
      // need to sort tmp_edges
      pbbslib::sample_sort (tmp_edges, deg, [&](const edge u, const edge v) {
        return std::get<0>(u) < std::get<0>(v);
      }, true);

      auto iter = vertex_ops::get_iter(tmp_edges, deg);
      size_t nbytes = byte::sequentialCompressEdgeSet<W>(
          edges.begin() + byte_offsets[rank[i]], 0, deg, rank[i], iter);
    }
  }, 1);

  auto out_vdata = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    out_vdata[i].offset = byte_offsets[i];
    out_vdata[i].degree = degrees[i];
  });
  byte_offsets.clear();

  auto deg_f = [&](size_t i) { return degrees[i]; };
  auto deg_map = pbbslib::make_sequence<size_t>(n, deg_f);
  uintT total_deg = pbbslib::reduce_add(deg_map);
  auto edge_arr = edges.to_array();
  //std::cout << "# Filtered, total_deg = " << total_deg << "\n";
  return symmetric_graph<csv_byte, W>(out_vdata, G.n, total_deg,
                            get_deletion_fn(out_vdata, edge_arr),
                            edge_arr, edge_arr);
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


