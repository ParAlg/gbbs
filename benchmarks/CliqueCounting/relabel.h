#pragma once

#include "gbbs/gbbs.h"
#include "gbbs/graph.h"

namespace gbbs {

template <
    template <class W> class vertex, class W, typename P,
    typename std::enable_if<
        std::is_same<vertex<W>, asymmetric_vertex<W>>::value, int>::type = 0>
inline auto relabel_graph(asymmetric_graph<vertex, W>& G, uintE* rank, P& pred)
    -> decltype(G) {
  std::cout << "Filter graph not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return G;
}

template <
    template <class W> class vertex, class W, typename P,
    typename std::enable_if<
        std::is_same<vertex<W>, cav_bytepd_amortized<W>>::value, int>::type = 0>
inline auto relabel_graph(asymmetric_graph<vertex, W>& G, uintE* rank, P& pred)
    -> decltype(G) {
  std::cout << "Filter graph not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return G;
}

template <
    template <class W> class vertex, class W, typename P,
    typename std::enable_if<
        std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value, int>::type = 0>
inline symmetric_graph<csv_byte, W> relabel_graph(
    symmetric_graph<vertex, W>& GA, uintE* rank, P& pred) {  // -> decltype(GA)
  size_t n = GA.n;
  using edge = std::tuple<uintE, W>;

  // 1. Calculate total size
  auto degrees = sequence<uintE>::uninitialized(n);
  auto byte_offsets = sequence<uintT>::uninitialized(n + 1);
  parallel_for(0, n, 1, [&](size_t i) {
    size_t total_bytes = 0;
    uintE last_ngh = 0;
    size_t deg = 0;
    uchar tmp[16];

    size_t prev_deg = GA.get_vertex(i).out_degree();
    // here write out all of G's outneighbors to an uncompressed array, and then
    // relabel and sort
    auto tmp_edges = sequence<edge>::uninitialized(prev_deg);
    auto f = [&](uintE u, uintE v, W w) {
      if (pred(u, v, w)) {
        tmp_edges[deg] = std::make_tuple(rank[v], w);
        deg++;
      }
      return false;
    };
    GA.get_vertex(i).out_neighbors().map(f, false);
    // need to sort tmp_edges
    tmp_edges.resize(deg);
    parlay::stable_sort_inplace(make_slice(tmp_edges),
                                [&](const edge u, const edge v) {
                                  return std::get<0>(u) < std::get<0>(v);
                                });

    // now need to compute total_bytes (TODO: this could be parallelized, b/c we
    // know what last_ngh is)
    for (size_t j = 0; j < deg; j++) {
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
    degrees[i] = deg;
  });
  byte_offsets[n] = 0;
  size_t last_offset = parlay::scan_inplace(byte_offsets);
  std::cout << "# size is: " << last_offset << "\n";

  auto edges = gbbs::new_array_no_init<uchar>(last_offset);

  // redo the sort from above, but now actually store your edges
  parallel_for(0, n, 1, [&](size_t i) {
    uintE deg = degrees[i];
    if (deg > 0) {
      // here write out all of G's outneighbors to an uncompressed array, and
      // then relabel and sort
      auto tmp_edges = sequence<edge>::uninitialized(deg);
      size_t j = 0;
      auto f = [&](uintE u, uintE v, W w) {
        if (pred(u, v, w)) {
          tmp_edges[j] = std::make_tuple(rank[v], w);
          j++;
        }
        return false;
      };
      GA.get_vertex(i).out_neighbors().map(f, false);
      // need to sort tmp_edges
      parlay::stable_sort_inplace(make_slice(tmp_edges),
                                  [&](const edge u, const edge v) {
                                    return std::get<0>(u) < std::get<0>(v);
                                  });

      auto iter = vertex_ops::get_iter(tmp_edges.begin(), deg);
      byte::sequentialCompressEdgeSet<W>(edges + byte_offsets[rank[i]], 0, deg,
                                         rank[i], iter);
    }
  });

  auto out_vdata = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    out_vdata[i].offset = byte_offsets[i];
    out_vdata[rank[i]].degree = degrees[i];
  });
  byte_offsets.clear();

  auto deg_f = [&](size_t i) { return degrees[i]; };
  auto deg_map = parlay::delayed_seq<size_t>(n, deg_f);
  uintT total_deg = parlay::reduce(deg_map);
  std::cout << "# Filtered, total_deg = " << total_deg << "\n";
  return symmetric_graph<csv_byte, W>(out_vdata, GA.n, total_deg,
                                      [=]() {
                                        gbbs::free_array(out_vdata, n);
                                        free_array(edges, last_offset);
                                      },
                                      edges);
}

template <
    template <class W> class vertex, class W, typename P,
    typename std::enable_if<std::is_same<vertex<W>, symmetric_vertex<W>>::value,
                            int>::type = 0>
inline symmetric_graph<symmetric_vertex, W> relabel_graph(
    symmetric_graph<vertex, W>& GA, uintE* rank, P& pred) {
  using w_vertex = vertex<W>;
  size_t n = GA.n;
  auto outOffsets = sequence<uintT>::uninitialized(n + 1);

  parallel_for(0, n, 1, [&](size_t i) {
    w_vertex u = GA.get_vertex(i);
    auto out_nghs = u.out_neighbors();
    auto out_f = [&](uintE j) {
      return static_cast<int>(
          pred(i, out_nghs.get_neighbor(j), out_nghs.get_weight(j)));
    };
    auto out_im = parlay::delayed_seq<int>(u.out_degree(), out_f);

    if (out_im.size() > 0)
      outOffsets[rank[i]] = parlay::reduce(out_im);
    else
      outOffsets[rank[i]] = 0;
  });

  outOffsets[n] = 0;
  uintT outEdgeCount = parlay::scan_inplace(outOffsets);

  using edge = std::tuple<uintE, W>;

  auto out_edges = gbbs::new_array_no_init<edge>(outEdgeCount);

  parallel_for(0, n, 1, [&](size_t i) {
    w_vertex u = GA.get_vertex(i);
    size_t out_offset = outOffsets[rank[i]];
    uintE d = u.out_degree();
    uintE true_deg = outOffsets[rank[i] + 1] - out_offset;
    if (d > 0) {
      edge* nghs = u.neighbors;
      edge* dir_nghs = out_edges + out_offset;
      auto pred_c = [&](const edge& e) {
        return pred(i, std::get<0>(e), std::get<1>(e));
      };
      auto n_im_f = [&](size_t j) { return nghs[j]; };
      auto n_im = parlay::delayed_seq<edge>(d, n_im_f);
      parlay::filter_out(n_im, gbbs::make_slice(dir_nghs, d), pred_c,
                         parlay::no_flag);
      parallel_for(0, true_deg, [&](size_t j) {
        dir_nghs[j] = std::make_tuple(rank[std::get<0>(dir_nghs[j])],
                                      std::get<1>(dir_nghs[j]));
      });
      parlay::stable_sort_inplace(gbbs::make_slice(dir_nghs, true_deg),
                                  [&](const edge left, const edge right) {
                                    return std::get<0>(left) <
                                           std::get<0>(right);
                                  });
    }
  });

  auto out_vdata = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    out_vdata[i].offset = outOffsets[i];
    out_vdata[i].degree = outOffsets[i + 1] - outOffsets[i];
  });
  outOffsets.clear();

  return symmetric_graph<symmetric_vertex, W>(out_vdata, GA.n, outEdgeCount,
                                              [=]() {
                                                gbbs::free_array(out_vdata, n);
                                                gbbs::free_array(out_edges,
                                                                 outEdgeCount);
                                              },
                                              out_edges);
}

template <class Graph>
void clr_sparsify_graph(Graph& GA, size_t denom, long seed) {
  using W = typename Graph::weight_type;
  size_t n = GA.n;
  // Color vertices with denom colors
  uintE numColors = std::max((size_t)1, denom);
  sequence<uintE> colors = sequence<uintE>::from_function(n, [&](size_t i) {
    return parlay::hash64_2((uintE)seed + i) % numColors;
  });
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    if (colors[u] == colors[v]) return 0;
    return 1;
  };
  // Mutates in-place.
  filter_edges(GA, pack_predicate);
}

}  // namespace gbbs
