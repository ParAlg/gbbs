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

#include <string>

#include <stdlib.h>
#include <fstream>
#include <functional>
#include <iostream>

#include "bridge.h"
#include "compressed_vertex.h"
#include "flags.h"
#include "vertex.h"

/* Compressed Sparse Row (CSR) based representation for symmetric graphs.
 * Takes two template parameters:
 * 1) vertex_type: vertex template, parametrized by the weight type associated with each edge
 * 2) W: the weight template
 * The graph is represented as an array of edges of type vertex_type::edge_type,
 * which is just a tuple<uintE, W>.*/
template <template <class W> class vertex_type, class W>
struct symmetric_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;

  vertex_data* v_data;

  /* Pointer to edges */
  edge_type* e0;
  /* Pointer to second copy of edges--relevant if using 2-socket NVM */
  edge_type* e1;

  /* number of vertices in G */
  size_t n;
  /* number of edges in G */
  size_t m;

  /* called to delete the graph */
  std::function<void()> deletion_fn;

symmetric_graph(vertex_data* v_data, size_t n, size_t m,
    std::function<void()> _deletion_fn, edge_type* _e0, edge_type* _e1=nullptr)
      : v_data(v_data),
        e0(_e0),
        e1(_e1),
        n(n),
        m(m),
        deletion_fn(_deletion_fn) {
    if (_e1 == nullptr) {
      e1 = e0; // handles NVM case when graph is stored in symmetric memory
    }
  }

  void del() {
    deletion_fn();
  }

#ifndef TWOSOCKETNVM
  vertex get_vertex(uintE i) {
    return vertex(e0, v_data[i]);
  }
#else
  vertex get_vertex(uintE i) {
    if (numanode() == 0) {
      return vertex(e0, v_data[i]);
    } else {
      return vertex(e1, v_data[i]);
    }
  }
#endif

  template <class P>
  uintE pack_neighbors(uintE id, P& p, std::tuple<uintE, W>* tmp) {
    uintE new_degree = get_vertex(id).packOutNgh(id, p, tmp);
    v_data[id].degree = new_degree; /* updates the degree */
    return new_degree;
  }

  /* degree must be <= old_degree */
  void decrease_vertex_degree(uintE id, uintE degree) {
    assert(degree <= v_data[id].degree);
    v_data[id].degree = degree;
  }

  void zero_vertex_degree(uintE id) {
    decrease_vertex_degree(id, 0);
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    parallel_for(0, n, [&](size_t i) {
      get_vertex(i).mapOutNgh(i, f, parallel_inner_map);
    }, 1);
  }

  // F : edge -> edge
  template <class F>
  void alter_edges(F f, bool parallel_inner_map = true) {
    abort(); /* unimplemented for CSR */
  }

  pbbs::sequence<std::tuple<uintE, uintE, W>> edges() {
    using g_edge = std::tuple<uintE, uintE, W>;
    auto degs = pbbs::sequence<size_t>(n, [&] (size_t i) {
      return get_vertex(i).getOutDegree();
    });
    size_t sum_degs = pbbslib::scan_add_inplace(degs.slice());
    assert(sum_degs == m);
    auto edges = pbbs::sequence<g_edge>(sum_degs);
    parallel_for(0, n, [&](size_t i) {
      size_t k = degs[i];
      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        edges[k++] = std::make_tuple(u, v, wgh);
      };
      get_vertex(i).mapOutNgh(i, map_f, false);
    }, 1);
    return edges;
  }
};

/* Compressed Sparse Row (CSR) based representation for asymmetric
 * graphs.  Note that the symmetric/asymmetric structures are pretty
 * similar, but defined separately. The purpose is to try and avoid
 * errors where an algorithm intended for symmetric graphs (e.g.,
 * biconnectivity) is not mistakenly called on a directed graph.
 *
 * Takes two template parameters:
 * 1) vertex_type: vertex template, parametrized by the weight type
 *    associated with each edge
 * 2) W: the weight template
 *
 * The graph is represented as an array of edges of type
 * vertex_type::edge_type, which is just a pair<uintE, W>.
 * */
template <template <class W> class vertex_type, class W>
struct asymmetric_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;

  /* number of vertices in G */
  size_t n;
  /* number of edges in G */
  size_t m;
  /* called to delete the graph */
  std::function<void()> deletion_fn;

  vertex_data* v_out_data;
  vertex_data* v_in_data;

  /* Pointer to out-edges */
  edge_type* out_edges_0;
  /* Pointer to second copy of out-edges--relevant if using 2-socket NVM */
  edge_type* out_edges_1;

  /* Pointer to in-edges */
  edge_type* in_edges_0;
  /* Pointer to second copy of in-edges--relevant if using 2-socket NVM */
  edge_type* in_edges_1;

#ifndef TWOSOCKETNVM
  vertex get_vertex(size_t i) {
    return vertex(out_edges_0, v_out_data[i], in_edges_0, v_in_data[i]);
  }
#else
  vertex get_vertex(size_t i) {
    if (numanode() == 0) {
      return vertex(out_edges_0, v_out_data[i], in_edges_0, v_in_data[i]);
    } else {
      return vertex(out_edges_1, v_out_data[i], in_edges_1, v_in_data[i]);
    }
  }
#endif

  asymmetric_graph(vertex_data* v_out_data,
      vertex_data* v_in_data,
      size_t n,
      size_t m,
      std::function<void()> _deletion_fn,
      edge_type* _out_edges_0,
      edge_type* _in_edges_0,
      edge_type* _out_edges_1=nullptr,
      edge_type* _in_edges_1=nullptr) :
    n(n),
    m(m),
    deletion_fn(_deletion_fn),
    v_out_data(v_out_data),
    v_in_data(v_in_data),
    out_edges_0(_out_edges_0),
    out_edges_1(_out_edges_1),
    in_edges_0(_in_edges_0),
    in_edges_1(_in_edges_1) {}

  void del() {
    deletion_fn();
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    parallel_for(0, n, [&](size_t i) {
      get_vertex(i).mapOutNgh(i, f, parallel_inner_map);
    }, 1);
  }

  /* degree must be <= old_degree */
  void decrease_out_degree(uintE i, uintE degree) {
    /* currently unimplemented --- implement when needed */
    assert(false); exit(0);
  }

  /* degree must be <= old_degree */
  void decrease_in_degree(uintE i, uintE degree) {
    /* currently unimplemented --- implement when needed */
    assert(false); exit(0);
  }
};

// Edge Array Representation
template <class W>
struct edge_array {
  using weight_type = W;
  using edge = std::tuple<uintE, uintE, W>;
  edge* E;
  // for sq matrices, num_rows == num_cols
  size_t num_rows; // n TODO: deprecate #rows/#cols
  size_t num_cols; // TODO deprecate

  size_t n;
  size_t m;

  // non_zeros is the #edges
  size_t non_zeros; // m TODO rename to "m"
  void del() { pbbslib::free_array(E); }
  edge_array(edge* _E, size_t r, size_t c, size_t nz)
      : E(_E), num_rows(r), num_cols(c), non_zeros(nz) {
    if (r != c) {
      std::cout << "# edge_array format currently expects square matrix" << std::endl;
      exit(0);
    }
    n = r;
    m = nz;
  }
  edge_array() {}
  size_t size() { return non_zeros; }

  pbbs::sequence<edge> to_seq() {
    auto ret = pbbs::sequence<edge>(E, non_zeros);
    non_zeros = 0; E = nullptr;
    return std::move(ret);
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    parallel_for(0, m, [&](size_t i) {
      uintE u, v; W w;
      std::tie(u,v, w) = E[i];
      f(u, v, w);
    }, 512);
  }

  // F : edge -> edge
  template <class F>
  void alter_edges(F f, bool parallel_inner_map = true) {
    parallel_for(0, m, [&](size_t i) {
      uintE u, v; W w;
      std::tie(u, v, w) = E[i];
      E[i] = f(u, v, w);
    }, 512);
  }

  template <class P>
  void filter_edges(P p) {
    auto in_seq = pbbslib::make_sequence<edge>(m, E);
    auto q = pbbs::filter(in_seq, p);
    size_t q_s = q.size();
    pbbs::free_array(E);
    E = q.to_array(); m = q_s;
  }

};

inline auto get_deletion_fn(void* a, void* b) -> std::function<void()> {
  auto df = [&](void* a, void* b) {
    pbbslib::free_array(a);
    pbbslib::free_array(b);
  };
  return std::bind(df, a, b);
}

inline auto get_deletion_fn(void* a, void* b, void* c) -> std::function<void()> {
  auto df = [&](void* a, void* b, void* c) {
    pbbslib::free_array(a);
    pbbslib::free_array(b);
    pbbslib::free_array(c);
  };
  return std::bind(df, a, b, c);
}

inline auto get_deletion_fn(void* a, void* b, void* c, void* d) -> std::function<void()> {
  auto df = [&](void* a, void* b, void* c, void* d) {
    pbbslib::free_array(a);
    pbbslib::free_array(b);
    pbbslib::free_array(c);
    pbbslib::free_array(d);
  };
  return std::bind(df, a, b, c, d);
}

template <class W, class Graph>
inline edge_array<W> to_edge_array(Graph& G) {
  using edge = std::tuple<uintE, uintE, W>;
  size_t n = G.n;
  auto sizes = pbbs::sequence<uintT>(n);
  parallel_for(0, n, [&] (size_t i) {
    sizes[i] = G.get_vertex(i).getOutDegree();
  });
  size_t m = pbbslib::scan_add_inplace(sizes.slice());
  assert(m == G.m);

  edge* arr = pbbs::new_array_no_init<edge>(m);
  parallel_for(0, n, [&] (size_t i) {
    size_t idx = 0;
    uintT offset = sizes[i];
    auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      arr[offset + idx] = std::make_tuple(u, v, wgh);
      idx++;
    };
    G.get_vertex(i).mapOutNgh(i, map_f, /* parallel = */false);
  });
  return edge_array<W>(arr, n, n, m);
}

// Mutates (sorts) the underlying array A containing a black-box description of
// an edge of typename A::value_type. The caller provides functions GetU, GetV, and GetW
// which extract the u, v, and weight of a (u,v,w) edge respective (if the edge
// is a std::tuple<uinte, uintE, W> this is just get<0>, ..<1>, ..<2>
// respectively.
// e.g.:
//   using edge = std::tuple<uintE, uintE, W>;
//   auto get_u = [&] (const edge& e) { return std::get<0>(e); };
//   auto get_v = [&] (const edge& e) { return std::get<1>(e); };
//   auto get_w = [&] (const edge& e) { return std::get<2>(e); };
//   auto G = sym_graph_from_edges<W>(coo1, get_u, get_v, get_w, 10, false);
template <class W, class EdgeSeq, class GetU, class GetV, class GetW>
inline symmetric_graph<symmetric_vertex, W> sym_graph_from_edges(
    EdgeSeq& A,
    size_t n,
    GetU& get_u,
    GetV& get_v,
    GetW& get_w,
    bool is_sorted = false) {
  using vertex = symmetric_vertex<W>;
  using edge_type = typename vertex::edge_type;
  size_t m = A.size();

  if (m == 0) {
    if (n == 0) {
      std::function<void()> del = []() {};
      return symmetric_graph<symmetric_vertex, W>(nullptr, 0, 0, del, nullptr);
    } else {
      auto v_data = pbbs::new_array_no_init<vertex_data>(n);
      par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
        v_data[i].offset = 0;
        v_data[i].degree = 0;
      });
      std::function<void()> del = get_deletion_fn(v_data, nullptr);
      return symmetric_graph<symmetric_vertex, W>(v_data, n, 0, del, nullptr);
    }
  }

  if (!is_sorted) {
    auto first = [](std::tuple<uintE, uintE, W> a) { return std::get<0>(a); };
    size_t bits = pbbslib::log2_up(n);
    pbbslib::integer_sort_inplace(A.slice(), first, bits);
  }

  auto starts = sequence<uintT>(n+1, (uintT) 0);

  using neighbor = std::tuple<uintE, W>;
  auto edges = sequence<neighbor>(m, [&](size_t i) {
    // Fuse loops over edges (check if this helps)
    if (i == 0 || (get_u(A[i]) != get_u(A[i - 1]))) {
      starts[get_u(A[i])] = i;
    }
    if (i != (m-1)) {
      uintE our_vtx = get_u(A[i]);
      uintE next_vtx = get_u(A[i+1]);
      if (our_vtx != next_vtx && (our_vtx + 1 != next_vtx)) {
        par_for(our_vtx+1, next_vtx, pbbslib::kSequentialForThreshold, [&] (size_t k) {
          starts[k] = i+1;
        });
      }
    }
    if (i == (m-1)) { /* last edge */
      starts[get_u(A[i]) + 1] = m;
    }
    return std::make_tuple(get_v(A[i]), get_w(A[i]));
  });

//  auto copy_f = [&] (const uintT& u, const uintT& v) -> uintT {
//    if (v == 0) { return u; }
//    return v;
//  };
//  auto copy_m = pbbs::make_monoid(copy_f, (uintT)0);
//  pbbs::scan_inplace(starts.slice(), copy_m, pbbs::fl_inplace);
//  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
//    uintT o = starts[i];
//    size_t degree = ((i == n - 1) ? m : starts[i + 1]) - o;
//    v[i].degree = degree;
//    v[i].neighbors = ((std::tuple<uintE, W>*)(edges.begin() + o));
//    /* sort each neighbor list if needed */
//  });
//  return graph<V>(v, n, m, get_deletion_fn(v, edges.to_array()));

  auto v_data = pbbs::new_array_no_init<vertex_data>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintT o = starts[i];
    v_data[i].offset = o;
    v_data[i].degree = (uintE)(((i == (n-1)) ? m : starts[i+1]) - o);
  });
  auto new_edge_arr = edges.to_array();
  return symmetric_graph<symmetric_vertex, W>(v_data, n, m, get_deletion_fn(v_data, new_edge_arr), (edge_type*)new_edge_arr);
}

template <class W>
inline symmetric_graph<symmetric_vertex, W> sym_graph_from_edges(
    pbbs::sequence<std::tuple<uintE, uintE, W>>& A,
    size_t n,
    bool is_sorted = false) {
  using edge = std::tuple<uintE, uintE, W>;
  auto get_u = [&] (const edge& e) { return std::get<0>(e); };
  auto get_v = [&] (const edge& e) { return std::get<1>(e); };
  auto get_w = [&] (const edge& e) { return std::get<2>(e); };
  return sym_graph_from_edges<W>(A, n, get_u, get_v, get_w, is_sorted);
}
