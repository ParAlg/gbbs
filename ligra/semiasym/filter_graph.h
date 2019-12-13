#pragma once

#include "graph.h"
#include "packed_graph.h"


// Used by MST and MaximalMatching
// Predicate returns three values:
// 0 : keep in graph, do not return in edge array
// 1 : remove from graph, do not return in edge array
// 2 : remove from graph, return in edge array
// Cost: O(n+m) work
template <template <class W> class vertex, class W, class P>
inline edge_array<W> filter_edges(graph<vertex<W>>& G, P& pred, const flags fl = 0) {
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
  parallel_for(0, n, [&] (size_t i) {
    auto res = G.get_vertex(i).reduceOutNgh(i, map_f, red_monoid);
    if (std::get<0>(res) > 0 || std::get<1>(res) > 0) {
      vtx_offs[i] = std::make_tuple(std::get<0>(res), std::get<1>(res),
                                    G.V[i].calculateOutTemporarySpaceBytes());
    } else {
      vtx_offs[i] = std::make_tuple(std::get<0>(res), std::get<1>(res), 0);
    }
  }, 1);
  reduce_t.stop(); reduce_t.reportTotal("reduce time");
  vtx_offs[n] = std::make_tuple(0, 0, 0);
  auto scan_f = [](const std::tuple<uintT, uintT, uintT>& l,
                   const std::tuple<uintT, uintT, uintT>& r) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r),
                           std::get<2>(l) + std::get<2>(r));
  };
  pbbslib::scan_inplace(vtx_offs.slice(), pbbslib::make_monoid(scan_f, std::make_tuple(0, 0, 0)));

  size_t total_space_bytes =
      std::get<2>(vtx_offs[n]);  // total space needed for all vertices
  size_t output_size = std::get<1>(vtx_offs[n]);
  std::cout << "# tmp bytes to allocate = " << total_space_bytes
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
      size_t deg = G.V[i].getOutDegree();
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
          G.V[i].filterOutNgh(i, pred_two, out_f, tmp_v);
        }
        // Pack out non-zero edges. This method updates the degree in G.
        if (n_to_pack < deg) {
          G.V[i].packOutNgh(i, pred_zero, tmp_v);
        } else {
          G.V[i].setOutDegree(0);
        }
      }
    });
  }
  auto degree_imap = pbbslib::make_sequence<size_t>(n,
      [&](size_t i) { return G.V[i].getOutDegree(); });

  G.m = pbbslib::reduce_add(degree_imap);
  std::cout << "# G.m is now = " << G.m << "\n";

  return edge_array<W>(arr.to_array(), n, n, arr.size());
}

// Used by MaximalMatching.
template <template <class W> class vertex, class W, class P>
inline edge_array<W> filter_all_edges(graph<vertex<W>>& G, P& p) {
  using edge = std::tuple<uintE, uintE, W>;
  size_t n = G.n;
  auto offs = sequence<std::tuple<uintT, uintT>>(n + 1);
  par_for(0, n, [&] (size_t i) {
    offs[i] = std::make_tuple(G.V[i].countOutNgh(i, p),
                              G.V[i].calculateOutTemporarySpace());
  });
  offs[n] = std::make_tuple(0, 0);
  auto scan_f = [](const std::tuple<uintT, uintT>& l,
                   const std::tuple<uintT, uintT>& r) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r));
  };
  pbbslib::scan_inplace(offs.slice(), pbbslib::make_monoid(scan_f, std::make_tuple(0, 0)));
  size_t total_space = std::get<1>(offs[n]);
  auto tmp = sequence<std::tuple<uintE, W>>(total_space);
  std::cout << "# tmp space allocated = " << total_space << "\n";

  size_t total_edges = std::get<0>(offs[n]);
  auto arr = sequence<edge>(total_edges);

  {
    auto for_inner = [&](size_t i) {
      size_t off = std::get<0>(offs[i]);
      if (G.V[i].getOutDegree() > 0) {
        std::tuple<uintE, W>* tmp_v = tmp.begin() + std::get<1>(offs[i]);
        auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
          arr[off + j] = std::make_tuple(i, std::get<0>(nw), std::get<1>(nw));
        };
        G.V[i].filterOutNgh(i, p, out_f, tmp_v);
        G.V[i].setOutDegree(0);
        G.V[i].setInDegree(0);
      }
    };
    par_for(0, n, [&] (size_t i) { for_inner(i); });
  }
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
    uintE ct = G.V[i].template reduceOutNgh<uintE>(i, map_f, red_monoid);
    if (ct > 0) {
      vtx_offs[i] = std::make_tuple(ct, G.V[i].calculateOutTemporarySpace());
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
        G.V[i].filterOutNgh(i, pred, out_f, tmp_v);
      }
    };
    par_for(0, n, [&] (size_t i) { for_inner(i); });
  }
  return edge_array<W>(output_arr.to_array(), n, n, output_arr.size());
}
