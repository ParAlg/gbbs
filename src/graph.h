// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <string>
#pragma once

#include <stdlib.h>
#include <fstream>
#include <functional>
#include <iostream>
#include "compressed_vertex.h"
#include "vertex.h"

#include "lib/get_time.h"
#include "lib/integer_sort.h"

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

template <class vertex>
struct graph {
  vertex* V;
  long n;
  long m;
  bool transposed;
  uintE* flags;
  std::function<void()> deletion_fn;
  std::function<graph<vertex>()> copy_fn;

  graph(vertex* _V, long _n, long _m, std::function<void()> _d,
        uintE* _flags = NULL)
      : V(_V), n(_n), m(_m), transposed(0), flags(_flags), deletion_fn(_d) {}

  graph(vertex* _V, long _n, long _m, std::function<void()> _d,
        std::function<graph<vertex>()> _c, uintE* _flags = NULL)
      : V(_V),
        n(_n),
        m(_m),
        flags(_flags),
        transposed(0),
        deletion_fn(_d),
        copy_fn(_c) {}

  auto copy() -> graph<vertex> { return copy_fn(); }

  void del() {
    if (flags != NULL) pbbs::free_array(flags);
    deletion_fn();
  }

  template <class F>
  void map_edges(F f) {
    par_for(0, n, [&] (size_t i) {
      V[i].mapOutNgh(i, f);
    });
  }

};

inline auto get_deletion_fn(void* V, void* edges) -> std::function<void()> {
  auto df = [&](void* V, void* edges) {
    pbbs::free_array(V);
    pbbs::free_array(edges);
  };
  return std::bind(df, V, edges);
}

inline auto get_deletion_fn(void* V, void* in_edges, void* out_edges)
    -> std::function<void()> {
  auto df = [&](void* V, void* in_edges, void* out_edges) {
    pbbs::free_array(V);
    pbbs::free_array(in_edges);
    pbbs::free_array(out_edges);
  };
  return std::bind(df, V, in_edges, out_edges);
}

template <class vertex, class E>
inline std::function<graph<vertex>()> get_copy_fn(vertex* V, E* edges, size_t n,
                                                  size_t m, size_t sizeofe) {
  auto df = [&](vertex* V, E* edges, size_t n, size_t m, size_t sizeofe) {
    auto NV = pbbs::new_array_no_init<vertex>(n);
    auto NE = pbbs::new_array_no_init<E>(sizeofe);
    par_for(0, sizeofe, pbbs::kSequentialForThreshold, [&] (size_t i)
                    { NE[i] = edges[i]; });
    par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i) {
      NV[i].setOutDegree(V[i].getOutDegree());
      size_t off = (V[i].getOutNeighbors() - edges);
      NV[i].setOutNeighbors(NE + off);
    });
    graph<vertex> G = graph<vertex>(NV, n, m, get_deletion_fn(NV, NE));
    G.copy_fn = get_copy_fn(NV, NE, n, m, sizeofe);
    //    std::cout << "Returning copied graph" << "\n";
    return G;
  };
  return std::bind(df, V, edges, n, m, sizeofe);
}

template <class vertex, class E>
inline std::function<graph<vertex>()> get_copy_fn(vertex* V, E* in_edges,
                                                  E* out_edges, size_t n,
                                                  size_t m, size_t m_in,
                                                  size_t m_out) {
  auto df = [&](vertex* V, E* in_edges, E* out_edges, size_t n, size_t m,
                size_t m_in, size_t m_out) {
    auto NV = pbbs::new_array_no_init<vertex>(n);
    auto Nin = pbbs::new_array_no_init<E>(m_in);
    auto Nout = pbbs::new_array_no_init<E>(m_out);
    par_for(0, m_in, pbbs::kSequentialForThreshold, [&] (size_t i)
                    { Nin[i] = in_edges[i]; });
    par_for(0, m_in, pbbs::kSequentialForThreshold, [&] (size_t i)
                    { Nout[i] = out_edges[i]; });
    par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i) {
      NV[i].setOutDegree(V[i].getOutDegree());
      NV[i].setInDegree(V[i].getInDegree());
      size_t out_off = (V[i].getOutNeighbors() - out_edges);
      NV[i].setOutNeighbors(Nout + out_off);
      size_t in_off = (V[i].getInNeighbors() - in_edges);
      NV[i].setOutNeighbors(Nin + in_off);
    });
    graph<vertex> G = graph<vertex>(
        V, n, m, get_deletion_fn((void*)NV, (void*)Nin, (void*)Nout));
    G.copy_fn = get_copy_fn(NV, Nin, Nout, n, m, m_in, m_out);
    return G;
  };
  return std::bind(df, V, in_edges, out_edges, n, m, m_in, m_out);
}


template <
    template <class W> class vertex, class W, typename P,
    typename std::enable_if<std::is_same<vertex<W>, symmetricVertex<W>>::value,
                            int>::type = 0>
inline graph<asymmetricVertex<W>> filter_graph(graph<vertex<W>>& G, P& pred) {
  using w_vertex = vertex<W>;
  size_t n = G.n;
  w_vertex* V = G.V;
  auto out_edge_sizes = sequence<uintT>(n + 1);
  auto in_edge_sizes = sequence<uintT>(n + 1);

  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i) {
    w_vertex u = V[i];
    auto out_f = [&](uintE j) {
      return static_cast<int>(pred(i, u.getOutNeighbor(j), u.getOutWeight(j)));
    };
    auto out_im = make_sequence<int>(u.getOutDegree(), out_f);
    auto in_f = [&](uintE j) {
      return static_cast<int>(pred(u.getInNeighbor(j), i, u.getInWeight(j)));
    };
    auto in_im = make_sequence<int>(u.getInDegree(), in_f);

    if (out_im.size() > 0)
      out_edge_sizes[i] = pbbs::reduce_add(out_im);
    else
      out_edge_sizes[i] = 0;
    if (in_im.size() > 0)
      in_edge_sizes[i] = pbbs::reduce_add(in_im);
    else
      in_edge_sizes[i] = 0;
  });

  out_edge_sizes[n] = 0;
  in_edge_sizes[n] = 0;
  long outEdgeCount = pbbs::scan_add(out_edge_sizes, out_edge_sizes);
  long inEdgeCount = pbbs::scan_add(in_edge_sizes, in_edge_sizes);

  assert(G.m / 2 == outEdgeCount);
  assert(G.m / 2 == inEdgeCount);

  using edge = std::tuple<uintE, W>;

  auto out_edges = sequence<edge>(outEdgeCount);
  auto in_edges = sequence<edge>(inEdgeCount);

  par_for(0, n, [&] (size_t i) {
    w_vertex u = V[i];
    size_t out_offset = out_edge_sizes[i];
    uintE d = u.getOutDegree();
    if (d > 0) {
      edge* nghs = u.getOutNeighbors();
      edge* dir_nghs = out_edges.start() + out_offset;
      auto pred_c = [&](const edge& e) {
        return pred(i, std::get<0>(e), std::get<1>(e));
      };
      auto n_im_f = [&](size_t i) { return nghs[i]; };
      auto n_im = make_sequence<edge>(d, n_im_f);
      auto res = pbbs::filter(n_im, pred_c, pbbs::no_flag, dir_nghs);
    }
  });

  par_for(0, n, [&] (size_t i) {
    w_vertex u = V[i];
    size_t in_offset = in_edge_sizes[i];
    uintE d = u.getInDegree();
    if (d > 0) {
      edge* nghs = u.getInNeighbors();
      edge* dir_nghs = in_edges.start() + in_offset;

      auto pred_c = [&](const edge& e) {
        return pred(std::get<0>(e), i, std::get<1>(e));
      };
      auto n_im_f = [&](size_t i) { return nghs[i]; };
      auto n_im = make_sequence<edge>(d, n_im_f);
      auto res = pbbs::filter(n_im, pred_c, pbbs::no_flag, dir_nghs);
    }
  });

  auto AV = pbbs::new_array_no_init<asymmetricVertex<W>>(n);
  par_for(0, n, [&] (size_t i) {
    uintT in_offset = in_edge_sizes[i];
    uintT out_offset = out_edge_sizes[i];
    AV[i] = asymmetricVertex<W>(
        in_edges.start() + in_offset, out_edges.start() + out_offset,
        in_edge_sizes[i + 1] - in_offset, out_edge_sizes[i + 1] - out_offset);
  });

  return graph<asymmetricVertex<W>>(
      AV, G.n, outEdgeCount,
      get_deletion_fn(AV, out_edges.get_array(), in_edges.get_array()));
}

struct print_t {
  bool srcTarg(uintE src, uintE ngh, uintE edgesRead) {
    std::cout << "src = " << src << " ngh = " << ngh << " edge# " << edgesRead
              << "\n";
    return true;
  }
};

// byte version
template <template <class W> class vertex, class W, typename P,
          typename std::enable_if<
              std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value ||
                  std::is_same<vertex<W>, csv_byte<W>>::value,
              int>::type = 0>
inline graph<cav_byte<W>> filter_graph(graph<vertex<W>>& G, P& pred) {
  using w_vertex = vertex<W>;
  size_t n = G.n;
  w_vertex* V = G.V;

  std::cout << "Filtering"
            << "\n";
  // 1. Calculate total size
  auto degrees = sequence<uintE>(n);
  auto byte_offsets = sequence<uintT>(n + 1);
  par_for(0, n, [&] (size_t i) {
    size_t total_bytes = 0;
    uintE last_ngh = 0;
    size_t deg = 0;
    uchar tmp[16];
    auto f = [&](uintE u, uintE v, W w) {
      if (pred(u, v, w)) {
        long bytes = 0;
        if (deg == 0) {
          bytes = encodings::byte::compressFirstEdge(tmp, bytes, u, v);
          bytes = encodings::byte::compressWeight<W>(tmp, bytes, w);
        } else {
          bytes = encodings::byte::compressEdge(tmp, bytes, v - last_ngh);
          bytes = encodings::byte::compressWeight<W>(tmp, bytes, w);
        }
        last_ngh = v;
        total_bytes += bytes;
        deg++;
      }
      return false;
    };
    V[i].mapOutNgh(i, f, false);

    degrees[i] = deg;
    byte_offsets[i] = total_bytes;
  });
  byte_offsets[n] = 0;
  size_t last_offset = pbbs::scan_add(byte_offsets, byte_offsets);
  std::cout << " size is: " << last_offset << "\n";

  auto edges = sequence<uchar>(last_offset);

  {
    auto for_inner = [&](size_t i) {
      uintE new_deg = degrees[i];
      if (new_deg > 0) {
        auto app_pred = [&](std::tuple<uintE, W> val) {
          return pred(i, std::get<0>(val), std::get<1>(val));
        };

        auto iter = V[i].getOutIter(i);
        auto f_it =
            ligra_utils::make_filter_iter<std::tuple<uintE, W>>(iter, app_pred);
        long nbytes = encodings::byte::sequentialCompressEdgeSet<W>(
            edges.start() + byte_offsets[i], 0, new_deg, i, f_it);
        if (nbytes != (byte_offsets[i + 1] - byte_offsets[i])) {
          std::cout << "degree is: " << new_deg << " nbytes should be: "
                    << (byte_offsets[i + 1] - byte_offsets[i])
                    << " but is: " << nbytes << "\n";
          assert(nbytes == (byte_offsets[i + 1] - byte_offsets[i]));
        }
      }
    };

    par_for(0, n, [&] (size_t i) { for_inner(i); });
  }

  auto AV = pbbs::new_array_no_init<cav_byte<W>>(n);
  par_for(0, n, [&] (size_t i) {
    size_t o = byte_offsets[i];
    uchar* our_edges = edges.start() + o;
    AV[i].inNeighbors = nullptr;
    AV[i].inDegree = 0;

    AV[i].outNeighbors = our_edges;
    AV[i].outDegree = degrees[i];
  });

  auto deg_f = [&](size_t i) { return degrees[i]; };
  auto deg_map = make_sequence<size_t>(n, deg_f);
  uintT total_deg = pbbs::reduce_add(deg_map);
  std::cout << "Filtered, total_deg = " << total_deg << "\n";
  return graph<cav_byte<W>>(AV, G.n, total_deg,
                            get_deletion_fn(AV, edges.get_array()));
}

template <
    template <class W> class vertex, class W, typename P,
    typename std::enable_if<std::is_same<vertex<W>, asymmetricVertex<W>>::value,
                            int>::type = 0>
inline auto filter_graph(graph<vertex<W>>& G, P& pred) -> decltype(G) {
  assert(false);  // Not implemented for directed graphs
  return G;
}

template <
    template <class W> class vertex, class W, typename P,
    typename std::enable_if<
        std::is_same<vertex<W>, cav_bytepd_amortized<W>>::value, int>::type = 0>
inline auto filter_graph(graph<vertex<W>>& G, P& pred) -> decltype(G) {
  assert(false);  // Not implemented for directed graphs
  return G;
}

// Edge Array Representation
template <class W>
struct edge_array {
  using edge = std::tuple<uintE, uintE, W>;
  edge* E;
  // for sq matrices, num_rows == num_cols
  size_t num_rows;
  size_t num_cols;
  // non_zeros is the #edges
  size_t non_zeros;
  void del() { pbbs::free_array(E); }
  edge_array(edge* _E, size_t r, size_t c, size_t nz)
      : E(_E), num_rows(r), num_cols(c), non_zeros(nz) {}
  edge_array() {}
  size_t size() { return non_zeros; }
};

// Used by MST and MaximalMatching
// Predicate returns three values:
// 0 : keep in graph, do not return in edge array
// 1 : remove from graph, do not return in edge array
// 2 : remove from graph, return in edge array
// Cost: O(n+m) work
template <template <class W> class vertex, class W, class P>
inline edge_array<W> filter_edges(graph<vertex<W>>& G, P& pred) {
  using edge = std::tuple<uintE, uintE, W>;
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
                  const std::tuple<uintT, uintT>& r) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r));
  };
  par_for(0, n, [&] (size_t i) {
    auto res = G.V[i].reduceOutNgh(i, id, map_f, red_f);
    if (std::get<0>(res) > 0 || std::get<1>(res) > 0) {
      vtx_offs[i] = std::make_tuple(std::get<0>(res), std::get<1>(res),
                                    G.V[i].calculateOutTemporarySpace());
    } else {
      vtx_offs[i] = std::make_tuple(std::get<0>(res), std::get<1>(res), 0);
    }
  });
  vtx_offs[n] = std::make_tuple(0, 0, 0);
  auto scan_f = [](const std::tuple<uintT, uintT, uintT>& l,
                   const std::tuple<uintT, uintT, uintT>& r) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r),
                           std::get<2>(l) + std::get<2>(r));
  };
  pbbs::scan(vtx_offs, vtx_offs, scan_f, std::make_tuple(0, 0, 0));

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

  std::cout << "starting pack my man " << "\n";
  // 2. pack and write out
  {
    auto for_inner = [&](size_t i) {
      size_t off = std::get<1>(vtx_offs[i]);
      size_t n_one = std::get<0>(vtx_offs[i + 1]) - std::get<0>(vtx_offs[i]);
      size_t n_two = std::get<1>(vtx_offs[i + 1]) - off;
      size_t n_to_pack = n_one + n_two;
      if (n_to_pack > 0) {
        std::tuple<uintE, W>* tmp_v = tmp.start() + std::get<2>(vtx_offs[i]);
        auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
          arr[off + j] = std::make_tuple(i, std::get<0>(nw), std::get<1>(nw));
        };
        // Filter out edges where pred == 2.
        if (n_two > 0) {
          G.V[i].filterOutNgh(i, pred_two, out_f, tmp_v);
        }
        // Pack out any non-zero edges
        G.V[i].packOutNgh(i, pred_zero, tmp_v);
      }
    };
    par_for(0, n, [&] (size_t i) { for_inner(i); });
  }
  std::cout << "packed up " << "\n";
  auto deg_f = [&](size_t i) { return G.V[i].getOutDegree(); };
  auto degree_imap =
      make_sequence<size_t>(n, deg_f);

  G.m = pbbs::reduce_add(degree_imap);
  std::cout << "G.m = " << G.m << "\n";

  return edge_array<W>(arr.get_array(), n, n, arr.size());
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
  //  std::cout << "fall e" << "\n";
  offs[n] = std::make_tuple(0, 0);
  auto scan_f = [](const std::tuple<uintT, uintT>& l,
                   const std::tuple<uintT, uintT>& r) {
    return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                           std::get<1>(l) + std::get<1>(r));
  };
  pbbs::scan(offs, offs, scan_f, std::make_tuple(0, 0));
  size_t total_space = std::get<1>(offs[n]);
  auto tmp = sequence<std::tuple<uintE, W>>(total_space);
  std::cout << "tmp space allocated = " << total_space << "\n";

  size_t total_edges = std::get<0>(offs[n]);
  auto arr = sequence<edge>(total_edges);

  {
    auto for_inner = [&](size_t i) {
      size_t off = std::get<0>(offs[i]);
      if (G.V[i].getOutDegree() > 0) {
        std::tuple<uintE, W>* tmp_v = tmp.start() + std::get<1>(offs[i]);
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
  //  std::cout << "G.m = " << G.m << "arr.size = " << arr.size() << "\n";
  G.m = 0;
  return edge_array<W>(arr.get_array(), n, n, arr.size());
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
  par_for(0, n, [&] (size_t i) {
    uintE ct = G.V[i].reduceOutNgh(i, id, map_f, red_f);
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
  pbbs::scan(vtx_offs, vtx_offs, scan_f, std::make_tuple(0, 0));

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
        std::tuple<uintE, W>* tmp_v = tmp.start() + std::get<1>(vtx_offs[i]);
        auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
          output_arr[off + j] =
              std::make_tuple(i, std::get<0>(nw), std::get<1>(nw));
        };
        G.V[i].filterOutNgh(i, pred, out_f, tmp_v);
      }
    };
    par_for(0, n, [&] (size_t i) { for_inner(i); });
  }
  return edge_array<W>(output_arr.get_array(), n, n, output_arr.size());
}

// Mutates (sorts) the underlying array
// Returns an unweighted, symmetric graph
template <class W>
inline graph<symmetricVertex<W>> sym_graph_from_edges(edge_array<W>& A,
                                                      bool is_sorted = false) {
  using V = symmetricVertex<W>;
  using edge = std::tuple<uintE, uintE, W>;
  size_t m = A.non_zeros;
  size_t n = std::max<size_t>(A.num_cols, A.num_rows);

  if (m == 0) {
    std::function<void()> del = []() {};
    if (n == 0) {
      return graph<V>(nullptr, 0, 0, del);
    } else {
      V* v = pbbs::new_array_no_init<V>(n);
      par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i) {
        v[i].degree = 0;
        v[i].neighbors = nullptr;
      });
      return graph<V>(v, n, 0, del);
    }
  }

  auto Am = make_sequence<edge>(A.E, m);
  if (!is_sorted) {
    auto first = [](std::tuple<uintE, uintE, W> a) { return std::get<0>(a); };
    size_t bits = pbbs::log2_up(n);
    pbbs::integer_sort(Am, Am, first, bits);
  }

  auto starts = sequence<uintT>(n);
  V* v = pbbs::new_array_no_init<V>(n);
  auto edges = sequence<uintE>(m, [&](size_t i) {
    // Fuse loops over edges (check if this helps)
    if (i == 0 || (std::get<0>(Am[i]) != std::get<0>(Am[i - 1]))) {
      starts[std::get<0>(Am[i])] = i;
    }
    return std::get<1>(Am[i]);
  });
  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i) {
    uintT o = starts[i];
    size_t degree = ((i == n - 1) ? m : starts[i + 1]) - o;
    v[i].degree = degree;
    v[i].neighbors = ((std::tuple<uintE, W>*)(edges.start() + o));
  });
  return graph<V>(v, n, m, get_deletion_fn(v, edges.get_array()));
}
