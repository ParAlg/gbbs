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

  /* degree must be <= old_degree */
  void decrease_degree(uintE i, uintE degree) {
    v_data[i].degree = degree;
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    parallel_for(0, n, [&](size_t i) {
      get_vertex(i).mapOutNgh(i, f, parallel_inner_map);
    }, 1);
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
    v_out_data(v_out_data),
    v_in_data(v_in_data),
    n(n),
    m(m),
    deletion_fn(_deletion_fn),
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
  using edge = std::tuple<uintE, uintE, W>;
  edge* E;
  // for sq matrices, num_rows == num_cols
  size_t num_rows; // n
  size_t num_cols;
  // non_zeros is the #edges
  size_t non_zeros; // m
  void del() { pbbslib::free_array(E); }
  edge_array(edge* _E, size_t r, size_t c, size_t nz)
      : E(_E), num_rows(r), num_cols(c), non_zeros(nz) {}
  edge_array() {}
  size_t size() { return non_zeros; }
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

// Mutates (sorts) the underlying array
// Returns an unweighted, symmetric graph
template <class W>
inline symmetric_graph<symmetric_vertex, W> sym_graph_from_edges(edge_array<W>& A,
                                                      bool is_sorted = false) {
  using wvertex = symmetric_vertex<W>;
  using edge = std::tuple<uintE, uintE, W>;
  using E = typename wvertex::E;
  size_t m = A.non_zeros;
  size_t n = std::max<size_t>(A.num_cols, A.num_rows);

  if (m == 0) {
    if (n == 0) {
      std::function<void()> del = []() {};
      return symmetric_graph<symmetric_vertex, W>(nullptr, 0, 0, del, nullptr);
    } else {
      uintT* offsets = pbbs::new_array_no_init<uintT>(n+1);
      par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
        offsets[i] = 0;
      });
      std::function<void()> del = get_deletion_fn(offsets, nullptr);
      return symmetric_graph<symmetric_vertex, W>(offsets, n, 0, del, nullptr);
    }
  }

  auto Am = pbbslib::make_sequence<edge>(A.E, m);
  if (!is_sorted) {
    auto first = [](std::tuple<uintE, uintE, W> a) { return std::get<0>(a); };
    size_t bits = pbbslib::log2_up(n);
    pbbslib::integer_sort_inplace(Am, first, bits);
  }

  auto starts = sequence<uintT>(n);
  uintT* offsets = pbbslib::new_array_no_init<uintT>(n+1);
  auto edges = sequence<uintE>(m, [&](size_t i) {
    // Fuse loops over edges (check if this helps)
    if (i == 0 || (std::get<0>(Am[i]) != std::get<0>(Am[i - 1]))) {
      starts[std::get<0>(Am[i])] = i;
    }
    return std::get<1>(Am[i]);
  });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintT o = starts[i];
    offsets[i]  = o;
  });
  offsets[n] = m;
  auto edge_arr = edges.to_array();
  return symmetric_graph<symmetric_vertex, W>(offsets, n, m, get_deletion_fn(offsets, edge_arr), (E*)edge_arr);
}
