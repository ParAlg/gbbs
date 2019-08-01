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
#include <functional>
#include <iostream>
#include <string>

#include "bridge.h"
#include "compressed_vertex.h"
#include "flags.h"
#include "pbbslib/parallel.h"
#include "vertex.h"

#ifdef NVM
#include <numa.h>
#include <utmpx.h>
#endif

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

template <template <class W> class vertex, class W>
struct symmetric_graph {
  using vtx_type = vertex<W>;
  using weight_type = W;
  using E = typename vtx_type::E;

  vertex_data* V;

  E* e0;
  E* e1;

  size_t n;
  size_t m;
  std::function<void()> deletion_fn;

#ifndef NVM
  vtx_type get_vertex(size_t i) {
    return vtx_type(V[i], e0);
  }
#else
  vtx_type get_vertex(size_t i) {
    if (numanode() == 0) {
      return vtx_type(V[i], e0);
    } else {
      return vtx_type(V[i], e1);
    }
  }
#endif

symmetric_graph(vertex_data* V, size_t n, size_t m, std::function<void()> _d,
        E* _e0, E* _e1=nullptr)
      : V(V),
        e0(_e0),
        e1(_e1),
        n(n),
        m(m),
        deletion_fn(_d) {
    if (_e1 == nullptr) {
      e1 = e0; // handles NVM case when graph is stored in symmetric memory
    }
  }

  void del() {
    deletion_fn();
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    par_for(0, n, 1,
            [&](size_t i) { get_vertex(i).mapOutNgh(i, f, parallel_inner_map); });
  }
};

template <template <class W> class vertex, class W>
struct asymmetric_graph {
  using vtx_type = vertex<W>;
  using weight_type = W;
  using E = typename vtx_type::E;

  vtx_type* V;

  size_t n;
  size_t m;
  std::function<void()> deletion_fn;

  vtx_type get_vertex(size_t i) {
    return V[i];
  }

  asymmetric_graph(vtx_type* V, size_t n, size_t m, std::function<void()> _d)
      : V(V),
        n(n),
        m(m),
        deletion_fn(_d) {
  }

  void del() {
    deletion_fn();
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    par_for(0, n, 1,
            [&](size_t i) { get_vertex(i).mapOutNgh(i, f, parallel_inner_map); });
  }
};



inline auto get_deletion_fn(void* V, void* edges) -> std::function<void()> {
  auto df = [&](void* V, void* edges) {
    pbbslib::free_array(V);
    pbbslib::free_array(edges);
  };
  return std::bind(df, V, edges);
}

inline auto get_deletion_fn(void* V, void* in_edges, void* out_edges)
    -> std::function<void()> {
  auto df = [&](void* V, void* in_edges, void* out_edges) {
    pbbslib::free_array(V);
    pbbslib::free_array(in_edges);
    pbbslib::free_array(out_edges);
  };
  return std::bind(df, V, in_edges, out_edges);
}

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

// Mutates (sorts) the underlying array
// Returns an unweighted, symmetric graph
template <class W>
inline symmetric_graph<symmetricVertex, W> sym_graph_from_edges(edge_array<W>& A,
                                                      bool is_sorted = false) {
  using wvertex = symmetricVertex<W>;
  using edge = std::tuple<uintE, uintE, W>;
  using E = typename wvertex::E;
  size_t m = A.non_zeros;
  size_t n = std::max<size_t>(A.num_cols, A.num_rows);

  if (m == 0) {
    std::function<void()> del = []() {};
    if (n == 0) {
      return symmetric_graph<symmetricVertex, W>(nullptr, 0, 0, del, nullptr);
    } else {
      vertex_data* v = pbbslib::new_array_no_init<vertex_data>(n);
      par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
        v[i].degree = 0;
        v[i].offset = 0;
      });
      return symmetric_graph<symmetricVertex, W>(v, n, 0, del, nullptr);
    }
  }

  auto Am = pbbslib::make_sequence<edge>(A.E, m);
  if (!is_sorted) {
    auto first = [](std::tuple<uintE, uintE, W> a) { return std::get<0>(a); };
    size_t bits = pbbslib::log2_up(n);
    pbbslib::integer_sort_inplace(Am, first, bits);
  }

  auto starts = sequence<uintT>(n);
  vertex_data* v = pbbslib::new_array_no_init<vertex_data>(n);
  auto edges = sequence<uintE>(m, [&](size_t i) {
    // Fuse loops over edges (check if this helps)
    if (i == 0 || (std::get<0>(Am[i]) != std::get<0>(Am[i - 1]))) {
      starts[std::get<0>(Am[i])] = i;
    }
    return std::get<1>(Am[i]);
  });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintT o = starts[i];
    size_t degree = ((i == n - 1) ? m : starts[i + 1]) - o;
    v[i].degree = degree;
    v[i].offset = o;
  });
  auto edge_arr = edges.to_array();
  return symmetric_graph<symmetricVertex, W>(v, n, m, get_deletion_fn(v, edge_arr), (E*)edge_arr);
}
