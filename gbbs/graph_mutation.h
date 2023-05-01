#pragma once

#include "bridge.h"
#include "compressed_vertex.h"
#include "edge_array.h"
#include "vertex.h"
#include "vertex_subset.h"

namespace gbbs {

/* Filters a symmetric graph, G, with a predicate function pred.  Note
 * that the predicate does not have to be symmetric, i.e. f(u,v) is
 * not necesssarily equal to f(v,u), but we only represent the out-edges of this
 * (possibly) directed graph. */
template <
    template <class W> class vertex, class W, class Graph, typename P,
    typename std::enable_if<std::is_same<vertex<W>, symmetric_vertex<W>>::value,
                            int>::type = 0>
inline std::tuple<size_t, size_t, vertex_data*,
                  typename symmetric_vertex<W>::edge_type*>
filter_graph(Graph& G, P& pred) {
  using w_vertex = vertex<W>;
  size_t n = G.num_vertices();
  auto outOffsets = sequence<uintT>(n + 1);

  parallel_for(0, n, 1, [&](size_t i) {
    w_vertex u = G.get_vertex(i);
    auto u_out_nghs = u.out_neighbors();
    auto out_f = [&](uintE j) {
      return static_cast<int>(
          pred(i, u_out_nghs.get_neighbor(j), u_out_nghs.get_weight(j)));
    };
    auto out_im = parlay::delayed_seq<int>(u.out_degree(), out_f);

    if (out_im.size() > 0)
      outOffsets[i] = parlay::reduce(out_im);
    else
      outOffsets[i] = 0;
  });

  outOffsets[n] = 0;
  uintT outEdgeCount = parlay::scan_inplace(outOffsets);

  // assert(G.m / 2 == outEdgeCount);

  using edge = std::tuple<uintE, W>;

  auto out_edges = gbbs::new_array_no_init<edge>(outEdgeCount);

  parallel_for(0, n, 1, [&](size_t i) {
    w_vertex u = G.get_vertex(i);
    size_t out_offset = outOffsets[i];
    uintE d = u.out_degree();
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
    }
  });

  auto out_vdata = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    out_vdata[i].offset = outOffsets[i];
    out_vdata[i].degree = outOffsets[i + 1] - outOffsets[i];
  });
  outOffsets.clear();

  return std::make_tuple(G.num_vertices(), outEdgeCount, out_vdata, out_edges);
}

// byte version
template <
    template <class W> class vertex, class W, class Graph, typename P,
    typename std::enable_if<
        std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value, int>::type = 0>
inline auto filter_graph(Graph& G, P& pred) {
  size_t n = G.num_vertices();

  gbbs_debug(std::cout << "# Filtering"
                  << "\n");

  // 1. Calculate total size
  auto degrees = sequence<uintE>(n);
  auto byte_offsets = sequence<uintT>(n + 1);
  parallel_for(0, n, 1, [&](size_t i) {
    size_t total_bytes = 0;
    uintE last_ngh = 0;
    size_t deg = 0;
    uchar tmp[16];
    auto f = [&](uintE u, uintE v, W w) {
      if (pred(u, v, w)) {
        size_t bytes = 0;
        if (deg == 0) {
          bytes = byte::compressFirstEdge(tmp, bytes, u, v);
          bytes = byte::compressWeight<W>(tmp, bytes, w);
        } else {
          bytes = byte::compressEdge(tmp, bytes, v - last_ngh);
          bytes = byte::compressWeight<W>(tmp, bytes, w);
        }
        last_ngh = v;
        total_bytes += bytes;
        deg++;
      }
      return false;
    };
    G.get_vertex(i).out_neighbors().map(f, false);

    degrees[i] = deg;
    byte_offsets[i] = total_bytes;
  });
  byte_offsets[n] = 0;
  size_t last_offset = parlay::scan_inplace(byte_offsets);
  std::cout << "# size is: " << last_offset << "\n";

  size_t edges_size = last_offset;
  auto edges = gbbs::new_array_no_init<uchar>(edges_size);

  parallel_for(0, n, 1, [&](size_t i) {
    uintE new_deg = degrees[i];
    if (new_deg > 0) {
      auto app_pred = [&](std::tuple<uintE, W> val) {
        return pred(i, std::get<0>(val), std::get<1>(val));
      };

      auto iter = G.get_vertex(i).out_neighbors().get_iter();
      auto f_it = gbbs::make_filter_iter<std::tuple<uintE, W>>(iter, app_pred);
      size_t nbytes = byte::sequentialCompressEdgeSet<W>(
          edges + byte_offsets[i], 0, new_deg, i, f_it);
      if (nbytes != (byte_offsets[i + 1] - byte_offsets[i])) {
        std::cout << "# degree is: " << new_deg << " nbytes should be: "
                  << (byte_offsets[i + 1] - byte_offsets[i])
                  << " but is: " << nbytes << "\n";
        assert(nbytes == (byte_offsets[i + 1] - byte_offsets[i]));
      }
    }
  });

  auto out_vdata = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    out_vdata[i].offset = byte_offsets[i];
    out_vdata[i].degree = degrees[i];
  });
  byte_offsets.clear();

  auto deg_f = [&](size_t i) { return degrees[i]; };
  auto deg_map = parlay::delayed_seq<size_t>(n, deg_f);
  uintT total_deg = parlay::reduce(deg_map);
  std::cout << "# Filtered, total_deg = " << total_deg << "\n";
  return std::make_tuple(G.num_vertices(), edges_size, out_vdata, edges);
}

template <
    template <class W> class vertex, class W, class Graph, typename P,
    typename std::enable_if<
        std::is_same<vertex<W>, asymmetric_vertex<W>>::value, int>::type = 0>
inline auto filter_graph(Graph& G, P& pred) -> decltype(G) {
  std::cout << "# Filter graph not implemented for directed graphs"
            << std::endl;
  assert(false);  // Not implemented for directed graphs
  return G;
}

template <
    template <class W> class vertex, class W, class Graph, typename P,
    typename std::enable_if<
        std::is_same<vertex<W>, cav_bytepd_amortized<W>>::value, int>::type = 0>
inline auto filter_graph(Graph& G, P& pred) -> decltype(G) {
  std::cout << "# Filter graph not implemented for directed graphs"
            << std::endl;
  assert(false);  // Not implemented for directed graphs
  return G;
}

// Used by MST and MaximalMatching
// Predicate returns three values:
// 0 : keep in graph, do not return in edge array
// 1 : remove from graph, do not return in edge array
// 2 : remove from graph, return in edge array
// Cost: O(n+m) work
template <class Graph, class P>
edge_array<typename Graph::weight_type> filter_edges(Graph& G, P& pred,
                                                     const flags fl = 0) {
  using W = typename Graph::weight_type;
  using edge = std::tuple<uintE, uintE, W>;
  size_t n = G.num_vertices();
  auto vtx_offs = sequence<std::tuple<size_t, size_t, size_t>>(n + 1);

  // 1. map and write the # filtered edges for each vtx into vtx_offs
  std::tuple<uintT, uintT> id =
      std::make_tuple(0, 0);  // #vals == 1, #vals == 2, space to allocate
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    uintE pr = pred(src, ngh, wgh);
    return std::make_tuple<uintT, uintT>(pr == 1, pr == 2);
  };
  auto red_f = [](const std::tuple<uintT, uintT>& l,
                  const std::tuple<uintT, uintT>& r)
      __attribute__((always_inline)) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r));
  };
  auto red_monoid = parlay::make_monoid(red_f, id);
  timer reduce_t;
  reduce_t.start();
  parallel_for(0, n, 1, [&](size_t i) {
    auto res = G.get_vertex(i).out_neighbors().reduce(map_f, red_monoid);
    if (std::get<0>(res) > 0 || std::get<1>(res) > 0) {
      vtx_offs[i] = std::make_tuple(
          std::get<0>(res), std::get<1>(res),
          G.get_vertex(i).out_neighbors().calculateTemporarySpace());
    } else {
      vtx_offs[i] = std::make_tuple(std::get<0>(res), std::get<1>(res), 0);
    }
  });
  reduce_t.stop();
  reduce_t.next("reduce time");
  vtx_offs[n] = std::make_tuple(0, 0, 0);
  auto scan_f = [](const std::tuple<uintT, uintT, uintT>& l,
                   const std::tuple<uintT, uintT, uintT>& r) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r),
                           std::get<2>(l) + std::get<2>(r));
  };
  parlay::scan_inplace(make_slice(vtx_offs),
                       parlay::make_monoid(scan_f, std::make_tuple(0, 0, 0)));

  size_t total_space =
      std::get<2>(vtx_offs[n]);  // total space needed for all vertices
  size_t output_size = std::get<1>(vtx_offs[n]);
  std::cout << "# tmp space to allocate = " << total_space
            << " output size = " << output_size << "\n";
  auto arr = sequence<edge>::uninitialized(output_size);
  auto tmp = sequence<std::tuple<uintE, W>>::uninitialized(total_space);

  auto pred_two = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return pred(src, ngh, wgh) == 2;
  };

  // Only keeps zero degree vertices (pred == 1 and pred == 2 are packed out)
  auto pred_zero = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return pred(src, ngh, wgh) == 0;
  };

  // 2. pack and write out
  parallel_for(0, n, 1, [&](size_t i) {
    size_t deg = G.get_vertex(i).out_degree();
    size_t off = std::get<1>(vtx_offs[i]);
    size_t n_one = std::get<0>(vtx_offs[i + 1]) - std::get<0>(vtx_offs[i]);
    size_t n_two = std::get<1>(vtx_offs[i + 1]) - off;
    size_t n_to_pack = n_one + n_two;
    if (n_to_pack > 0) {
      std::tuple<uintE, W>* tmp_v = tmp.begin() + std::get<2>(vtx_offs[i]);
      auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
        arr[off + j] = std::make_tuple(i, std::get<0>(nw), std::get<1>(nw));
      };
      // Filter out edges where pred == 2.
      if (n_two > 0 && !(fl & no_output)) {
        G.get_vertex(i).out_neighbors().filter(pred_two, out_f, tmp_v);
      }
      // Pack out non-zero edges. This method updates the degree in G.
      if (n_to_pack < deg) {
        G.packNeighbors(i, pred_zero, (uint8_t*)tmp_v);
      } else {
        G.zeroVertexDegree(i);
      }
    }
  });
  auto degree_imap = parlay::delayed_seq<size_t>(
      n, [&](size_t i) { return G.get_vertex(i).out_degree(); });

  G.m = parlay::reduce(degree_imap);
  std::cout << "# G.m is now = " << G.m << "\n";

  return edge_array<W>(std::move(arr), n);
}

// Used by MaximalMatching.
template <class Graph, class P>
edge_array<typename Graph::weight_type> filter_all_edges(Graph& G, P& p,
                                                         flags fl = 0) {
  using W = typename Graph::weight_type;
  using edge = std::tuple<uintE, uintE, W>;
  size_t n = G.n;
  auto offs = sequence<std::tuple<uintT, uintT>>(n + 1);
  parallel_for(0, n, [&](size_t i) {
    offs[i] = std::make_tuple(
        G.get_vertex(i).out_neighbors().count(p),
        G.get_vertex(i).out_neighbors().calculateTemporarySpace());
  });
  //  std::cout << "fall e" << "\n";
  offs[n] = std::make_tuple(0, 0);
  auto scan_f = [](const std::tuple<uintT, uintT>& l,
                   const std::tuple<uintT, uintT>& r) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r));
  };
  parlay::scan_inplace(make_slice(offs),
                       parlay::make_monoid(scan_f, std::make_tuple(0, 0)));
  size_t total_space = std::get<1>(offs[n]);
  auto tmp = sequence<std::tuple<uintE, W>>(total_space);
  std::cout << "# tmp space allocated = " << total_space << "\n";

  size_t total_edges = std::get<0>(offs[n]);
  auto arr = sequence<edge>(total_edges);

  parallel_for(0, n, 1, [&](size_t i) {
    size_t off = std::get<0>(offs[i]);
    if (G.get_vertex(i).out_degree() > 0) {
      std::tuple<uintE, W>* tmp_v = tmp.begin() + std::get<1>(offs[i]);
      auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
        arr[off + j] = std::make_tuple(i, std::get<0>(nw), std::get<1>(nw));
      };
      G.get_vertex(i).out_neighbors().filter(p, out_f, tmp_v);
      G.zeroVertexDegree(i);
    }
  });
  //  std::cout << "G.m = " << G.m << "arr.size = " << arr.size() << "\n";
  G.m = 0;
  return edge_array<W>(std::move(arr), n);
}

// Similar to filter_edges, except we only filter (no packing). Any edge s.t.
// pred(src, ngh, wgh) == 1 is returned in the output edge array.
template <class Graph, class P>
edge_array<typename Graph::weight_type> sample_edges(Graph& G, P& pred) {
  using W = typename Graph::weight_type;
  using edge = std::tuple<uintE, uintE, W>;
  size_t n = G.num_vertices();
  auto vtx_offs = sequence<std::tuple<size_t, size_t>>(n + 1);

  // 1. Compute the # filtered edges and tmp-space required for each vtx.
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return pred(src, ngh, wgh);
  };
  auto red_monoid = parlay::plus<size_t>();
  parallel_for(0, n, [&](size_t i) {
    uintE ct = G.get_vertex(i).out_neighbors().reduce(map_f, red_monoid);
    if (ct > 0) {
      vtx_offs[i] = std::make_tuple(
          ct, G.get_vertex(i).out_neighbors().calculateTemporarySpace());
    } else {
      vtx_offs[i] = std::make_tuple(0, 0);
    }
  });
  vtx_offs[n] = std::make_tuple(0, 0);
  auto scan_f = [](const std::tuple<size_t, size_t>& l,
                   const std::tuple<size_t, size_t>& r) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r));
  };
  parlay::scan_inplace(make_slice(vtx_offs),
                       parlay::make_monoid(scan_f, std::make_tuple(0, 0)));

  size_t output_size = std::get<0>(vtx_offs[n]);
  auto output_arr = sequence<edge>(output_size);

  size_t tmp_space = std::get<1>(vtx_offs[n]);
  auto tmp = sequence<std::tuple<uintE, W>>(tmp_space);

  // 2. Filter edges into output arr.
  {
    parallel_for(0, n, 1, [&](size_t i) {
      size_t off = std::get<0>(vtx_offs[i]);
      size_t n_to_pack = std::get<0>(vtx_offs[i + 1]) - off;
      if (n_to_pack > 0) {
        std::tuple<uintE, W>* tmp_v = tmp.begin() + std::get<1>(vtx_offs[i]);
        auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
          output_arr[off + j] =
              std::make_tuple(i, std::get<0>(nw), std::get<1>(nw));
        };
        G.get_vertex(i).out_neighbors().filter(pred, out_f, tmp_v);
      }
    });
  }
  return edge_array<W>(std::move(output_arr), n);
}

// Packs out the adjacency lists of all vertex in vs. A neighbor, ngh, is kept
// in the new adjacency list if p(ngh) is true.
template <class Graph, class P>
inline vertexSubsetData<uintE> packEdges(Graph& G, vertexSubset& vs, P& p,
                                         const flags& fl = 0) {
  using S = std::tuple<uintE, uintE>;

  vs.toSparse();
  size_t m = vs.numNonzeros();
  size_t n = vs.numRows();
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  auto space = sequence<uintT>(m + 1);
  parallel_for(0, m, [&](size_t i) {
    uintE v = vs.vtx(i);
    space[i] = G.get_vertex(v).out_neighbors().calculateTemporarySpaceBytes();
  });
  space[m] = 0;
  size_t total_space = parlay::scan_inplace(space);
  uint8_t* tmp = nullptr;
  if (total_space > 0) {
    tmp = gbbs::new_array_no_init<uint8_t>(total_space);
    gbbs_debug(std::cout << "Allocated " << total_space << " temporary space."
                    << std::endl;);
  }
  if (should_output(fl)) {
    auto outV = sequence<S>::uninitialized(vs.size());
    parallel_for(0, m, 1, [&](size_t i) {
      uintE v = vs.vtx(i);
      uint8_t* tmp_v = nullptr;
      if (space[i + 1] > space[i]) {
        tmp_v = tmp + space[i];
      }
      uintE new_degree = G.packNeighbors(v, p, tmp_v);
      outV[i] = std::make_tuple(v, new_degree);
    });
    if (tmp) {
      gbbs::free_array(tmp, total_space);
    }
    return vertexSubsetData<uintE>(n, std::move(outV));
  } else {
    parallel_for(0, m, 1, [&](size_t i) {
      uintE v = vs.vtx(i);
      uint8_t* tmp_v = nullptr;
      if (space[i + 1] > space[i]) {
        tmp_v = tmp + space[i];
      }
      G.packNeighbors(v, p, tmp_v);
    });
    if (tmp) {
      gbbs::free_array(tmp, total_space);
    }
    return vertexSubsetData<uintE>(n);
  }
}

template <class Graph, class P>
inline vertexSubsetData<uintE> edgeMapFilter(Graph& G, vertexSubset& vs, P& p,
                                             const flags& fl = 0) {
  vs.toSparse();
  if (fl & pack_edges) {
    return packEdges(G, vs, p, fl);
  }
  size_t m = vs.numNonzeros();
  size_t n = vs.numRows();
  using S = std::tuple<uintE, uintE>;
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  if (should_output(fl)) {
    auto outV = sequence<S>(vs.size());
    parallel_for(0, m, 1, [&](size_t i) {
      uintE v = vs.vtx(i);
      size_t ct = G.get_vertex(v).out_neighbors().count(p);
      outV[i] = std::make_tuple(v, ct);
    });
    return vertexSubsetData<uintE>(n, std::move(outV));
  } else {
    parallel_for(0, m, 1, [&](size_t i) {
      uintE v = vs.vtx(i);
      G.get_vertex(v).out_neighbors().count(p);
    });
    return vertexSubsetData<uintE>(n);
  }
}

}  // namespace gbbs
