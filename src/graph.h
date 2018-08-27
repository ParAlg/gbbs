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

#pragma once

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "compressed_vertex.h"
#include "vertex.h"

#include "../lib/integer_sort.h"
using namespace std;

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
      : V(_V), n(_n), m(_m), flags(_flags), transposed(0), deletion_fn(_d) {}

  graph(vertex* _V, long _n, long _m, std::function<void()> _d,
        std::function<graph<vertex>()> _c, uintE* _flags = NULL)
      : V(_V),
        n(_n),
        m(_m),
        flags(_flags),
        transposed(0),
        deletion_fn(_d),
        copy_fn(_c) {}

  auto copy() { return copy_fn(); }

  void del() {
    if (flags != NULL) free(flags);
    deletion_fn();
  }
};

auto get_deletion_fn(void* V, void* edges) {
  auto df = [&](void* V, void* edges) {
    free(V);
    free(edges);
  };
  return std::bind(df, V, edges);
}

auto get_deletion_fn(void* V, void* in_edges, void* out_edges) {
  auto df = [&](void* V, void* in_edges, void* out_edges) {
    free(V);
    free(in_edges);
    free(out_edges);
  };
  return std::bind(df, V, in_edges, out_edges);
}

template <class vertex, class E>
std::function<graph<vertex>()> get_copy_fn(vertex* V, E* edges, size_t n,
                                           size_t m, size_t sizeofe) {
  auto df = [&](vertex* V, E* edges, size_t n, size_t m, size_t sizeofe) {
    auto NV = newA(vertex, n);
    auto NE = newA(E, sizeofe);
    parallel_for_bc(i, 0, sizeofe, (sizeofe > pbbs::kSequentialForThreshold),
                    { NE[i] = edges[i]; });
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      NV[i].setOutDegree(V[i].getOutDegree());
      size_t off = (V[i].getOutNeighbors() - edges);
      NV[i].setOutNeighbors(NE + off);
    });
    graph<vertex> G = graph<vertex>(NV, n, m, get_deletion_fn(NV, NE));
    G.copy_fn = get_copy_fn(NV, NE, n, m, sizeofe);
    //    cout << "Returning copied graph" << endl;
    return G;
  };
  return std::bind(df, V, edges, n, m, sizeofe);
}

template <class vertex, class E>
std::function<graph<vertex>()> get_copy_fn(vertex* V, E* in_edges, E* out_edges,
                                           size_t n, size_t m, size_t m_in,
                                           size_t m_out) {
  auto df = [&](vertex* V, E* in_edges, E* out_edges, size_t n, size_t m,
                size_t m_in, size_t m_out) {
    auto NV = newA(vertex, n);
    auto Nin = newA(E, m_in);
    auto Nout = newA(E, m_out);
    parallel_for_bc(i, 0, m_in, (m_in > pbbs::kSequentialForThreshold),
                    { Nin[i] = in_edges[i]; });
    parallel_for_bc(i, 0, m_in, (m_in > pbbs::kSequentialForThreshold),
                    { Nout[i] = out_edges[i]; });
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
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

auto get_copy_fn(void* V, void* in_edges, void* out_edges) {
  auto df = [&](void* V, void* in_edges, void* out_edges) {};
  return std::bind(df, V, in_edges, out_edges);
}

template <
    template <class W> typename vertex, class W, typename P,
    typename std::enable_if<std::is_same<vertex<W>, symmetricVertex<W>>::value,
                            int>::type = 0>
auto filter_graph(graph<vertex<W>>& G, P& pred) {
  using w_vertex = vertex<W>;
  size_t n = G.n;
  w_vertex* V = G.V;
  auto out_edge_sizes = array_imap<uintT>(n + 1);
  auto in_edge_sizes = array_imap<uintT>(n + 1);

  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    w_vertex u = V[i];
    auto out_im = make_in_imap<int>(u.getOutDegree(), [&](uintE j) {
      return static_cast<int>(pred(i, u.getOutNeighbor(j), u.getOutWeight(j)));
    });
    auto in_im = make_in_imap<int>(u.getInDegree(), [&](uintE j) {
      return static_cast<int>(pred(u.getInNeighbor(j), i, u.getInWeight(j)));
    });

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

  using edge = tuple<uintE, W>;

  auto out_edges = array_imap<edge>(outEdgeCount);
  auto in_edges = array_imap<edge>(inEdgeCount);

  parallel_for_bc(i, 0, n, true, {
    w_vertex u = V[i];
    size_t out_offset = out_edge_sizes[i];
    uintE d = u.getOutDegree();
    if (d > 0) {
      edge* nghs = u.getOutNeighbors();
      edge* dir_nghs = out_edges.start() + out_offset;
      auto pred_c = [&](const edge& e) {
        return pred(i, get<0>(e), get<1>(e));
      };
      auto n_im = make_in_imap<edge>(d, [&](size_t i) { return nghs[i]; });
      auto res = pbbs::filter(n_im, pred_c, pbbs::no_flag, dir_nghs);
      size_t k = res.size();
      assert(k == (out_edge_sizes[i + 1] - out_offset));
    }
  });

  parallel_for_bc(i, 0, n, true, {
    w_vertex u = V[i];
    size_t in_offset = in_edge_sizes[i];
    uintE d = u.getInDegree();
    if (d > 0) {
      edge* nghs = u.getInNeighbors();
      edge* dir_nghs = in_edges.start() + in_offset;

      auto pred_c = [&](const edge& e) {
        return pred(get<0>(e), i, get<1>(e));
      };
      auto n_im = make_in_imap<edge>(d, [&](size_t i) { return nghs[i]; });
      auto res = pbbs::filter(n_im, pred_c, pbbs::no_flag, dir_nghs);
      size_t k = res.size();
      assert(k == (in_edge_sizes[i + 1] - in_offset));
    }
  });

  auto AV = newA(asymmetricVertex<W>, n);
  parallel_for_bc(i, 0, n, true, {
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
    cout << "src = " << src << " ngh = " << ngh << " edge# " << edgesRead
         << endl;
    return true;
  }
};

// byte version
template <
    template <class W> typename vertex, class W, typename P,
    typename std::enable_if<
        std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value, int>::type = 0>
auto filter_graph(graph<vertex<W>>& G, P& pred) {
  using w_vertex = vertex<W>;
  size_t n = G.n;
  w_vertex* V = G.V;

  cout << "Filtering" << endl;
  // 1. Calculate total size
  auto degrees = array_imap<uintE>(n);
  auto byte_offsets = array_imap<uintT>(n + 1);
  parallel_for_bc(i, 0, n, true, {
    size_t total_bytes = 0;
    uintE last_ngh = 0;
    size_t deg = 0;
    uchar tmp[16];
    auto f = [&](uintE u, uintE v, W w) {
      if (pred(u, v, w)) {
        long bytes = 0;
        if (deg == 0) {
          bytes = compression::compressFirstEdge(tmp, bytes, u, v);
          bytes = compression::compressWeight<W>(tmp, bytes, w);
        } else {
          bytes = compression::compressEdge(tmp, bytes, v - last_ngh);
          bytes = compression::compressWeight<W>(tmp, bytes, w);
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
  cout << " size is: " << last_offset << endl;

  auto edges = array_imap<uchar>(last_offset);

  parallel_for_bc(i, 0, n, true, {
    uintE new_deg = degrees[i];
    if (new_deg > 0) {
      auto app_pred = [&](tuple<uintE COMMA W> val) {
        return pred(i, get<0>(val), get<1>(val));
      };

      auto iter = V[i].getOutIter(i);
      auto f_it =
          ligra_utils::make_filter_iter<tuple<uintE COMMA W>>(iter, app_pred);
      long nbytes = byte::sequentialCompressEdgeSet<W>(
          edges.start() + byte_offsets[i], 0, new_deg, i, f_it);
      if (nbytes != (byte_offsets[i + 1] - byte_offsets[i])) {
        cout << "degree is: " << new_deg
             << " nbytes should be: " << (byte_offsets[i + 1] - byte_offsets[i])
             << " but is: " << nbytes << endl;
        assert(nbytes == (byte_offsets[i + 1] - byte_offsets[i]));
      }
    }
  });

  auto AV = newA(cav_byte<W>, n);
  parallel_for_bc(i, 0, n, true, {
    size_t o = byte_offsets[i];
    uchar* our_edges = edges.start() + o;
    AV[i].inNeighbors = nullptr;
    AV[i].inDegree = 0;

    AV[i].outNeighbors = our_edges;
    AV[i].outDegree = degrees[i];
  });

  auto deg_map = make_in_imap<size_t>(n, [&](size_t i) { return degrees[i]; });
  uintT total_deg = pbbs::reduce_add(deg_map);
  cout << "Filtered, total_deg = " << total_deg << endl;
  return graph<cav_byte<W>>(AV, G.n, total_deg,
                            get_deletion_fn(AV, edges.get_array()));
}

template <
    template <class W> typename vertex, class W, typename P,
    typename std::enable_if<std::is_same<vertex<W>, asymmetricVertex<W>>::value,
                            int>::type = 0>
auto filter_graph(graph<vertex<W>>& G, P& pred) {
  exit(0);
  return G;
}

template <
    template <class W> typename vertex, class W, typename P,
    typename std::enable_if<
        std::is_same<vertex<W>, cav_bytepd_amortized<W>>::value, int>::type = 0>
auto filter_graph(graph<vertex<W>>& G, P& pred) {
  exit(0);
  return G;
}

// Edge Array Representation
template <class W>
struct edge_array {
  using edge = tuple<uintE, uintE, W>;
  edge* E;
  // for sq matrices, num_rows == num_cols
  size_t num_rows;
  size_t num_cols;
  // non_zeros is the #edges
  size_t non_zeros;
  void del() { free(E); }
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
edge_array<W> filter_edges(graph<vertex<W>>& G, P& pred) {
  using edge = tuple<uintE, uintE, W>;
  size_t n = G.n;
  auto vtx_offs = array_imap<tuple<size_t, size_t, size_t>>(n + 1);

  // 1. map and write the # filtered edges for each vtx into vtx_offs
  tuple<uintT, uintT> id =
      make_tuple(0, 0);  // #vals == 1, #vals == 2, space to allocate
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    uintE pr = pred(src, ngh, wgh);
    return make_tuple<uintT, uintT>(pr == 1, pr == 2);
  };
  auto red_f = [](const tuple<uintT, uintT>& l, const tuple<uintT, uintT>& r) {
    return make_tuple(get<0>(l) + get<0>(r), get<1>(l) + get<1>(r));
  };
  parallel_for_bc(i, 0, n, true, {
    auto res = G.V[i].reduceOutNgh(i, id, map_f, red_f);
    if (get<0>(res) > 0 || get<1>(res) > 0) {
      vtx_offs[i] = make_tuple(get<0>(res), get<1>(res),
                               G.V[i].calculateOutTemporarySpace());
    } else {
      vtx_offs[i] = make_tuple(get<0>(res), get<1>(res), 0);
    }
  });
  vtx_offs[n] = make_tuple(0, 0, 0);
  auto scan_f = [](const tuple<uintT, uintT, uintT>& l,
                   const tuple<uintT, uintT, uintT>& r) {
    return make_tuple(get<0>(l) + get<0>(r), get<1>(l) + get<1>(r),
                      get<2>(l) + get<2>(r));
  };
  pbbs::scan(vtx_offs, vtx_offs, scan_f, make_tuple(0, 0, 0));

  size_t total_space =
      get<2>(vtx_offs[n]);  // total space needed for all vertices
  size_t output_size = get<1>(vtx_offs[n]);
  cout << "tmp space to allocate = " << total_space
       << " output size = " << output_size << endl;
  auto arr = array_imap<edge>(output_size);
  auto tmp = array_imap<tuple<uintE, W>>(total_space);

  auto pred_two = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return pred(src, ngh, wgh) == 2;
  };

  // Only keeps zero degree vertices (pred == 1 and pred == 2 are packed out)
  auto pred_zero = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return pred(src, ngh, wgh) == 0;
  };

  // 2. pack and write out
  parallel_for_bc(i, 0, n, true, {
    size_t off = get<1>(vtx_offs[i]);
    size_t n_one = get<0>(vtx_offs[i + 1]) - get<0>(vtx_offs[i]);
    size_t n_two = get<1>(vtx_offs[i + 1]) - off;
    size_t n_to_pack = n_one + n_two;
    if (n_to_pack > 0) {
      tuple<uintE COMMA W>* tmp_v = tmp.start() + get<2>(vtx_offs[i]);
      auto out_f = [&](size_t j, const tuple<uintE, W>& nw) {
        arr[off + j] = make_tuple(i, get<0>(nw), get<1>(nw));
      };
      // Filter out edges where pred == 2.
      if (n_two > 0) {
        G.V[i].filterOutNgh(i, pred_two, out_f, tmp_v);
      }
      // Pack out any non-zero edges
      G.V[i].packOutNgh(i, pred_zero, tmp_v);
    }
  });
  auto degree_imap =
      make_in_imap<size_t>(n, [&](size_t i) { return G.V[i].getOutDegree(); });

  G.m = pbbs::reduce_add(degree_imap);
  return edge_array<W>(arr.get_array(), n, n, arr.size());
}

// Used by MaximalMatching.
template <template <class W> class vertex, class W, class P>
edge_array<W> filter_all_edges(graph<vertex<W>>& G, P& p) {
  using edge = tuple<uintE, uintE, W>;
  size_t n = G.n;
  auto offs = array_imap<tuple<uintT, uintT>>(n + 1);
  parallel_for_bc(i, 0, n, true, {
    offs[i] = make_tuple(G.V[i].countOutNgh(i, p),
                         G.V[i].calculateOutTemporarySpace());
  });
  //  cout << "fall e" << endl;
  offs[n] = make_tuple(0, 0);
  auto scan_f = [](const tuple<uintT, uintT>& l, const tuple<uintT, uintT>& r) {
    return make_tuple(get<0>(l) + get<0>(r), get<1>(l) + get<1>(r));
  };
  pbbs::scan(offs, offs, scan_f, make_tuple(0, 0));
  size_t total_space = get<1>(offs[n]);
  auto tmp = array_imap<tuple<uintE, W>>(total_space);
  cout << "tmp space allocated = " << total_space << endl;

  size_t total_edges = get<0>(offs[n]);
  auto arr = array_imap<edge>(total_edges);

  parallel_for_bc(i, 0, n, true, {
    size_t off = get<0>(offs[i]);
    if (G.V[i].getOutDegree() > 0) {
      tuple<uintE COMMA W>* tmp_v = tmp.start() + get<1>(offs[i]);
      auto out_f = [&](size_t j, const tuple<uintE, W>& nw) {
        arr[off + j] = make_tuple(i, get<0>(nw), get<1>(nw));
      };
      G.V[i].filterOutNgh(i, p, out_f, tmp_v);
      G.V[i].setOutDegree(0);
      G.V[i].setInDegree(0);
    }
  });
  //  cout << "G.m = " << G.m << "arr.size = " << arr.size() << endl;
  G.m = 0;
  return edge_array<W>(arr.get_array(), n, n, arr.size());
}

// Similar to filter_edges, except we only filter (no packing). Any edge s.t.
// pred(src, ngh, wgh) == 1 is returned in the output edge array.
template <template <class W> class vertex, class W, class P>
edge_array<W> sample_edges(graph<vertex<W>>& G, P& pred) {
  using edge = tuple<uintE, uintE, W>;
  size_t n = G.n;
  auto vtx_offs = array_imap<tuple<size_t, size_t>>(n + 1);

  // 1. Compute the # filtered edges and tmp-space required for each vtx.
  uintE id = 0;
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return pred(src, ngh, wgh);
  };
  auto red_f = [](size_t l, size_t r) { return l + r; };
  parallel_for_bc(i, 0, n, true, {
    uintE ct = G.V[i].reduceOutNgh(i, id, map_f, red_f);
    if (ct > 0) {
      vtx_offs[i] = make_tuple(ct, G.V[i].calculateOutTemporarySpace());
    } else {
      vtx_offs[i] = make_tuple(0, 0);
    }
  });
  vtx_offs[n] = make_tuple(0, 0);
  auto scan_f = [](const tuple<size_t, size_t>& l,
                   const tuple<size_t, size_t>& r) {
    return make_tuple(get<0>(l) + get<0>(r), get<1>(l) + get<1>(r));
  };
  pbbs::scan(vtx_offs, vtx_offs, scan_f, make_tuple(0, 0));

  size_t output_size = get<0>(vtx_offs[n]);
  auto output_arr = array_imap<edge>(output_size);

  size_t tmp_space = get<1>(vtx_offs[n]);
  auto tmp = array_imap<tuple<uintE, W>>(tmp_space);

  // 2. Filter edges into output arr.
  parallel_for_bc(i, 0, n, true, {
    size_t off = get<0>(vtx_offs[i]);
    size_t n_to_pack = get<0>(vtx_offs[i + 1]) - off;
    if (n_to_pack > 0) {
      tuple<uintE COMMA W>* tmp_v = tmp.start() + get<1>(vtx_offs[i]);
      auto out_f = [&](size_t j, const tuple<uintE, W>& nw) {
        output_arr[off + j] = make_tuple(i, get<0>(nw), get<1>(nw));
      };
      G.V[i].filterOutNgh(i, pred, out_f, tmp_v);
    }
  });
  return edge_array<W>(output_arr.get_array(), n, n, output_arr.size());
}

// Mutates (sorts) the underlying array
// Returns an unweighted, symmetric graph
template <class W>
auto sym_graph_from_edges(edge_array<W>& A, bool is_sorted = false) {
  using V = symmetricVertex<W>;
  using edge = tuple<uintE, uintE, W>;
  size_t m = A.non_zeros;
  size_t n = max<size_t>(A.num_cols, A.num_rows);

  if (m == 0) {
    std::function<void()> del = []() {};
    if (n == 0) {
      return graph<V>(nullptr, 0, 0, del);
    } else {
      V* v = pbbs::new_array_no_init<V>(n);
      parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
        v[i].degree = 0;
        v[i].neighbors = nullptr;
      });
      return graph<V>(v, n, 0, del);
    }
  }

  auto Am = make_array_imap<edge>(A.E, m);
  if (!is_sorted) {
    auto first = [](tuple<uintE, uintE, W> a) { return get<0>(a); };
    size_t bits = pbbs::log2_up(n);
    pbbs::integer_sort<uintE>(Am, Am, first, bits);
  }

  auto starts = array_imap<uintT>(n);
  V* v = pbbs::new_array_no_init<V>(n);
  auto edges = array_imap<uintE>(m, [&](size_t i) {
    // Fuse loops over edges (check if this helps)
    if (i == 0 || (get<0>(Am[i]) != get<0>(Am[i - 1]))) {
      starts[get<0>(Am[i])] = i;
    }
    return get<1>(Am[i]);
  });
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    uintT o = starts[i];
    size_t degree = ((i == n - 1) ? m : starts[i + 1]) - o;
    v[i].degree = degree;
    v[i].neighbors = ((tuple<uintE, W>*)(edges.start() + o));
  });
  return graph<V>(v, n, m, get_deletion_fn(v, edges.get_array()));
}
