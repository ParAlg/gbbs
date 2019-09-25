#pragma once

#include "bridge.h"
#include "compressed_vertex.h"
#include "graph.h"
#include "vertex.h"

/* Filters a symmetric graph, G, with a predicate function pred.  Note
 * that the predicate does not have to be symmetric, i.e. f(u,v) is
 * not necesssarily equal to f(v,u). The output graph is therefore
 * directed. */
template <template <class W> class vertex, class W, typename P,
    typename std::enable_if<std::is_same<vertex<W>, symmetric_vertex<W>>::value,
                            int>::type = 0>
inline asymmetric_graph<asymmetric_vertex, W> filter_graph(symmetric_graph<vertex, W>& G, P& pred) {
  using w_vertex = vertex<W>;
  size_t n = G.n;
  auto outOffsets = sequence<uintT>(n + 1);
  auto inOffsets = sequence<uintT>(n + 1);

  par_for(0, n, 1, [&] (size_t i) {
    w_vertex u = G.get_vertex(i);
    auto out_f = [&](uintE j) {
      return static_cast<int>(pred(i, u.getOutNeighbor(j), u.getOutWeight(j)));
    };
    auto out_im = pbbslib::make_sequence<int>(u.getOutDegree(), out_f);
    auto in_f = [&](uintE j) {
      return static_cast<int>(pred(u.getInNeighbor(j), i, u.getInWeight(j)));
    };
    auto in_im = pbbslib::make_sequence<int>(u.getInDegree(), in_f);

    if (out_im.size() > 0)
      outOffsets[i] = pbbslib::reduce_add(out_im);
    else
      outOffsets[i] = 0;
    if (in_im.size() > 0)
      inOffsets[i] = pbbslib::reduce_add(in_im);
    else
      inOffsets[i] = 0;
  });

  outOffsets[n] = 0;
  inOffsets[n] = 0;
  uintT outEdgeCount = pbbslib::scan_add_inplace(outOffsets);
  uintT inEdgeCount = pbbslib::scan_add_inplace(inOffsets);

  assert(G.m / 2 == outEdgeCount);
  assert(G.m / 2 == inEdgeCount);

  using edge = std::tuple<uintE, W>;

  auto out_edges = sequence<edge>(outEdgeCount);
  auto in_edges = sequence<edge>(inEdgeCount);

  par_for(0, n, 1, [&] (size_t i) {
    w_vertex u = G.get_vertex(i);
    size_t out_offset = outOffsets[i];
    uintE d = u.getOutDegree();
    if (d > 0) {
      edge* nghs = u.getOutNeighbors();
      edge* dir_nghs = out_edges.begin() + out_offset;
      auto pred_c = [&](const edge& e) {
        return pred(i, std::get<0>(e), std::get<1>(e));
      };
      auto n_im_f = [&](size_t i) { return nghs[i]; };
      auto n_im = pbbslib::make_sequence<edge>(d, n_im_f);
      pbbslib::filter_out(n_im, pbbslib::make_sequence(dir_nghs, d), pred_c, pbbslib::no_flag);
    }
  });

  par_for(0, n, 1, [&] (size_t i) {
    w_vertex u = G.get_vertex(i);
    size_t in_offset = inOffsets[i];
    uintE d = u.getInDegree();
    if (d > 0) {
      edge* nghs = u.getInNeighbors();
      edge* dir_nghs = in_edges.begin() + in_offset;

      auto pred_c = [&](const edge& e) {
        return pred(std::get<0>(e), i, std::get<1>(e));
      };
      auto n_im_f = [&](size_t i) { return nghs[i]; };
      auto n_im = pbbslib::make_sequence<edge>(d, n_im_f);
      pbbslib::filter_out(n_im, pbbslib::make_sequence(dir_nghs, d), pred_c, pbbslib::no_flag);
    }
  });

  auto out_edge_arr = out_edges.to_array();
  auto in_edge_arr = in_edges.to_array();
  return graph<asymmetric_vertex<W>>(
      outOffsets, inOffsets, G.n, outEdgeCount,
      get_deletion_fn(outOffsets, inOffsets, out_edge_arr, in_edge_arr),
      out_edge_arr, in_edge_arr);
}

// byte version
template <template <class W> class vertex, class W, typename P,
          typename std::enable_if<
              std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value,
              int>::type = 0>
inline graph<cav_byte<W>> filter_graph(graph<vertex<W>>& G, P& pred) {
  using w_vertex = vertex<W>;
  size_t n = G.n;
  w_vertex* V = G.V;

  debug(std::cout << "Filtering" << "\n");

  // 1. Calculate total size
  auto degrees = sequence<uintE>(n);
  auto byte_offsets = sequence<uintT>(n + 1);
  par_for(0, n, 1, [&] (size_t i) {
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
    G.get_vertex(i).mapOutNgh(i, f, false);

    degrees[i] = deg;
    byte_offsets[i] = total_bytes;
  });
  byte_offsets[n] = 0;
  size_t last_offset = pbbslib::scan_add_inplace(byte_offsets);
  std::cout << " size is: " << last_offset << "\n";

  auto edges = sequence<uchar>(last_offset);

  par_for(0, n, 1, [&] (size_t i) {
    uintE new_deg = degrees[i];
    if (new_deg > 0) {
      auto app_pred = [&](std::tuple<uintE, W> val) {
        return pred(i, std::get<0>(val), std::get<1>(val));
      };

      auto iter = G.get_vertex(i).getOutIter(i);
      auto f_it =
          pbbslib::make_filter_iter<std::tuple<uintE, W>>(iter, app_pred);
      size_t nbytes = byte::sequentialCompressEdgeSet<W>(
          edges.begin() + byte_offsets[i], 0, new_deg, i, f_it);
      if (nbytes != (byte_offsets[i + 1] - byte_offsets[i])) {
        std::cout << "degree is: " << new_deg << " nbytes should be: "
                  << (byte_offsets[i + 1] - byte_offsets[i])
                  << " but is: " << nbytes << "\n";
        assert(nbytes == (byte_offsets[i + 1] - byte_offsets[i]));
      }
    }
  });

  auto deg_f = [&](size_t i) { return degrees[i]; };
  auto deg_map = pbbslib::make_sequence<size_t>(n, deg_f);
  uintT total_deg = pbbslib::reduce_add(deg_map);
  auto edge_arr = edges.to_array();
  std::cout << "Filtered, total_deg = " << total_deg << "\n";
  return graph<cav_byte<W>>(byte_offsets, byte_offsets, G.n, total_deg,
                            get_deletion_fn(byte_offsets, edge_arr),
                            edge_arr, edge_arr);
}

template <
    template <class W> class vertex, class W, typename P,
    typename std::enable_if<std::is_same<vertex<W>, asymmetric_vertex<W>>::value,
                            int>::type = 0>
inline auto filter_graph(asymmetric_graph<vertex<W>>& G, P& pred) -> decltype(G) {
  std::cout << "Filter graph not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return G;
}

template <
    template <class W> class vertex, class W, typename P,
    typename std::enable_if<
        std::is_same<vertex<W>, cav_bytepd_amortized<W>>::value, int>::type = 0>
inline auto filter_graph(asymmetric_graph<vertex<W>>& G, P& pred) -> decltype(G) {
  std::cout << "Filter graph not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return G;
}

// Used by MST and MaximalMatching
// Predicate returns three values:
// 0 : keep in graph, do not return in edge array
// 1 : remove from graph, do not return in edge array
// 2 : remove from graph, return in edge array
// Cost: O(n+m) work
template <template <class W> class vertex, class W, class P>
inline edge_array<W> filter_edges(symmetric_graph<vertex, W>& G, P& pred, const flags fl = 0) {
  using edge = std::tuple<uintE, uintE, W>;
  using T = std::tuple<uintT, uintT>;
  size_t n = G.n;
  auto vtx_offs = sequence<std::tuple<size_t, size_t, size_t>>(n + 1);

  // 1. map and write the # filtered edges for each vtx into vtx_offs
  std::tuple<uintT, uintT> id =
      std::make_tuple(0, 0);  // #vals == 1, #vals == 2, space to allocate
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    uintE pr = pred(src, ngh, wgh);
    return std::make_tuple<uintT, uintT>(pr == 1, pr == 2);
  };
  auto red_f = [](const std::tuple<uintT, uintT>& l,
                  const std::tuple<uintT, uintT>& r) __attribute__((always_inline)) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r));
  };
  auto red_monoid = pbbslib::make_monoid(red_f, id);
  timer reduce_t; reduce_t.start();
  par_for(0, n, 1, [&] (size_t i) {
    auto res = G.get_vertex(i).template reduceOutNgh<T>(i, map_f, red_monoid);
    if (std::get<0>(res) > 0 || std::get<1>(res) > 0) {
      vtx_offs[i] = std::make_tuple(std::get<0>(res), std::get<1>(res),
                                    G.get_vertex(i).calculateOutTemporarySpace());
    } else {
      vtx_offs[i] = std::make_tuple(std::get<0>(res), std::get<1>(res), 0);
    }
  });
  reduce_t.stop(); reduce_t.reportTotal("reduce time");
  vtx_offs[n] = std::make_tuple(0, 0, 0);
  auto scan_f = [](const std::tuple<uintT, uintT, uintT>& l,
                   const std::tuple<uintT, uintT, uintT>& r) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r),
                           std::get<2>(l) + std::get<2>(r));
  };
  pbbslib::scan_inplace(vtx_offs.slice(), pbbslib::make_monoid(scan_f, std::make_tuple(0, 0, 0)));

  size_t total_space =
      std::get<2>(vtx_offs[n]);  // total space needed for all vertices
  size_t output_size = std::get<1>(vtx_offs[n]);
  std::cout << "tmp space to allocate = " << total_space
            << " output size = " << output_size << "\n";
  auto arr = sequence<edge>(output_size);
  auto tmp = sequence<std::tuple<uintE, W>>(total_space);

  auto pred_two = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return pred(src, ngh, wgh) == 2;
  };

  // Only keeps zero degree vertices (pred == 1 and pred == 2 are packed out)
  auto pred_zero = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return pred(src, ngh, wgh) == 0;
  };

  // 2. pack and write out
  {
    par_for(0, n, 1, [&] (size_t i) {
      size_t deg = G.get_vertex(i).getOutDegree();
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
          G.get_vertex(i).filterOutNgh(i, pred_two, out_f, tmp_v);
        }
        // Pack out non-zero edges. This method updates the degree in G.
        if (n_to_pack < deg) {
          G.get_vertex(i).packOutNgh(i, pred_zero, tmp_v);
        } else {
          G.get_vertex(i).setOutDegree(0);
        }
      }
    });
  }
  auto degree_imap = pbbslib::make_sequence<size_t>(n,
      [&](size_t i) { return G.get_vertex(i).getOutDegree(); });

  G.m = pbbslib::reduce_add(degree_imap);
  std::cout << "G.m is now = " << G.m << "\n";

  return edge_array<W>(arr.to_array(), n, n, arr.size());
}

// Used by MaximalMatching.
template <template <class W> class vertex, class W, class P>
inline edge_array<W> filter_all_edges(graph<vertex<W>>& G, P& p) {
  using edge = std::tuple<uintE, uintE, W>;
  size_t n = G.n;
  auto offs = sequence<std::tuple<uintT, uintT>>(n + 1);
  par_for(0, n, [&] (size_t i) {
    offs[i] = std::make_tuple(G.get_vertex(i).countOutNgh(i, p),
                              G.get_vertex(i).calculateOutTemporarySpace());
  });
  //  std::cout << "fall e" << "\n";
  offs[n] = std::make_tuple(0, 0);
  auto scan_f = [](const std::tuple<uintT, uintT>& l,
                   const std::tuple<uintT, uintT>& r) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r));
  };
  pbbslib::scan_inplace(offs.slice(), pbbslib::make_monoid(scan_f, std::make_tuple(0, 0)));
  size_t total_space = std::get<1>(offs[n]);
  auto tmp = sequence<std::tuple<uintE, W>>(total_space);
  std::cout << "tmp space allocated = " << total_space << "\n";

  size_t total_edges = std::get<0>(offs[n]);
  auto arr = sequence<edge>(total_edges);

  {
    auto for_inner = [&](size_t i) {
      size_t off = std::get<0>(offs[i]);
      if (G.get_vertex(i).getOutDegree() > 0) {
        std::tuple<uintE, W>* tmp_v = tmp.begin() + std::get<1>(offs[i]);
        auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
          arr[off + j] = std::make_tuple(i, std::get<0>(nw), std::get<1>(nw));
        };
        G.get_vertex(i).filterOutNgh(i, p, out_f, tmp_v);
        G.get_vertex(i).setOutDegree(0); /* TODO */
        G.get_vertex(i).setInDegree(0);
      }
    };
    par_for(0, n, [&] (size_t i) { for_inner(i); });
  }
  //  std::cout << "G.m = " << G.m << "arr.size = " << arr.size() << "\n";
  G.m = 0;
  return edge_array<W>(arr.to_array(), n, n, arr.size());
}

// Similar to filter_edges, except we only filter (no packing). Any edge s.t.
// pred(src, ngh, wgh) == 1 is returned in the output edge array.
template <template <class W> class vertex, class W, class P>
inline edge_array<W> sample_edges(graph<vertex<W>>& G, P& pred) {
  using edge = std::tuple<uintE, uintE, W>;
  size_t n = G.n;
  auto vtx_offs = sequence<std::tuple<size_t, size_t>>(n + 1);

  // 1. Compute the # filtered edges and tmp-space required for each vtx.
  uintE id = 0;
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return pred(src, ngh, wgh);
  };
  auto red_f = [](size_t l, size_t r) { return l + r; };
  auto red_monoid = pbbslib::make_monoid(red_f, id);
  par_for(0, n, [&] (size_t i) {
    uintE ct = G.get_vertex(i).template reduceOutNgh<uintE>(i, map_f, red_monoid);
    if (ct > 0) {
      vtx_offs[i] = std::make_tuple(ct, G.get_vertex(i).calculateOutTemporarySpace());
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
  pbbslib::scan_inplace(vtx_offs.slice(), pbbslib::make_monoid(scan_f, std::make_tuple(0, 0)));

  size_t output_size = std::get<0>(vtx_offs[n]);
  auto output_arr = sequence<edge>(output_size);

  size_t tmp_space = std::get<1>(vtx_offs[n]);
  auto tmp = sequence<std::tuple<uintE, W>>(tmp_space);

  // 2. Filter edges into output arr.
  {
    auto for_inner = [&](size_t i) {
      size_t off = std::get<0>(vtx_offs[i]);
      size_t n_to_pack = std::get<0>(vtx_offs[i + 1]) - off;
      if (n_to_pack > 0) {
        std::tuple<uintE, W>* tmp_v = tmp.begin() + std::get<1>(vtx_offs[i]);
        auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
          output_arr[off + j] =
              std::make_tuple(i, std::get<0>(nw), std::get<1>(nw));
        };
        G.get_vertex(i).filterOutNgh(i, pred, out_f, tmp_v);
      }
    };
    par_for(0, n, [&] (size_t i) { for_inner(i); });
  }
  return edge_array<W>(output_arr.to_array(), n, n, output_arr.size());
}



// Packs out the adjacency lists of all vertex in vs. A neighbor, ngh, is kept
// in the new adjacency list if p(ngh) is true.
template <template <class W> class wvertex, class W, class P>
inline void packAllEdges(graph<wvertex<W>>& G, P& p, const flags& fl = 0) {
  using vertex = wvertex<W>;
  vertex* G = G.V;
  size_t n = G.n;
  auto space = sequence<uintT>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    space[i] = G[i].calculateOutTemporarySpace();
  });
  size_t total_space = pbbslib::scan_add_inplace(space);
  auto tmp = sequence<std::tuple<uintE, W>>(total_space);

  auto for_inner = [&](size_t i) {
    std::tuple<uintE, W>* tmp_v = tmp.begin() + space[i];
    G[i].packOutNgh(i, p, tmp_v);
  };
  par_for(0, n, 1, [&] (size_t i) { for_inner(i); });
}


// Packs out the adjacency lists of all vertex in vs. A neighbor, ngh, is kept
// in the new adjacency list if p(ngh) is true.
template <class Graph, class P>
inline vertexSubsetData<uintE> packEdges(Graph& G,
                                         vertexSubset& vs, P& p,
                                         const flags& fl = 0) {
  using S = std::tuple<uintE, uintE>;
  vs.toSparse();
  size_t m = vs.numNonzeros();
  size_t n = vs.numRows();
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  auto space = sequence<uintT>(m);
  par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintE v = vs.vtx(i);
    space[i] = G.get_vertex(v).calculateOutTemporarySpace();
  });
  size_t total_space = pbbslib::scan_add_inplace(space);
  //std::cout << "packNghs: total space allocated = " << total_space << "\n";
  auto tmp = sequence<std::tuple<uintE, W>>(
      total_space);  // careful when total_space == 0
  S* outV;
  if (should_output(fl)) {
    outV = pbbslib::new_array_no_init<S>(vs.size());
    {
      auto for_inner = [&](size_t i) {
        uintE v = vs.vtx(i);
        std::tuple<uintE, W>* tmp_v = tmp.begin() + space[i];

        uintE new_degree = G.get_vertex(v).packOutNgh(v, p, tmp_v);
        G.decrease_degree(new_degree);

        outV[i] = std::make_tuple(v, new_degree);
      };
      par_for(0, m, 1, [&] (size_t i) { for_inner(i); });
    }
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    {
      auto for_inner = [&](size_t i) {
        uintE v = vs.vtx(i);
        std::tuple<uintE, W>* tmp_v = tmp.begin() + space[i];

        uintE new_degree = G[v].packOutNgh(v, p, tmp_v);
        G.decrease_degree(new_degree);

      };
      par_for(0, m, 1, [&] (size_t i) { for_inner(i); });
    }
    return vertexSubsetData<uintE>(n);
  }
}

template <class Graph, class P>
inline vertexSubsetData<uintE> edgeMapFilter(Graph& G,
                                             vertexSubset& vs, P& p,
                                             const flags& fl = 0) {
  vs.toSparse();
  if (fl & pack_edges) {
    return packEdges<Graph, P>(G, vs, p, fl);
  }
  size_t m = vs.numNonzeros();
  size_t n = vs.numRows();
  using S = std::tuple<uintE, uintE>;
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  S* outV;
  if (should_output(fl)) {
    outV = pbbslib::new_array_no_init<S>(vs.size());
  }
  if (should_output(fl)) {
    par_for(0, m, 1, [&] (size_t i) {
      uintE v = vs.vtx(i);
      size_t ct = G.get_vertex(v).countOutNgh(v, p);
      outV[i] = std::make_tuple(v, ct);
    });
  } else {
    par_for(0, m, 1, [&] (size_t i) {
      uintE v = vs.vtx(i);
      G.get_vertex(v).countOutNgh(v, p);
    });
  }
  if (should_output(fl)) {
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    return vertexSubsetData<uintE>(n);
  }
}
