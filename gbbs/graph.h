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
#include <optional>
#include <string>

#include "bridge.h"
#include "compressed_vertex.h"
#include "edge_array.h"
#include "flags.h"
#include "macros.h"
#include "vertex.h"

namespace gbbs {

//  Compressed Sparse Row (CSR) based representation for symmetric graphs.
//  Takes two template parameters:
//  1) vertex_type: vertex template, parametrized by the weight type associated
//  with each edge
//  2) W: the weight template
//  The graph is represented as an array of edges of type
//  vertex_type::edge_type.
//  For uncompressed vertices, this type is equal to tuple<uintE, W>.
template <template <class W> class vertex_type, class W>
struct symmetric_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;
  using graph = symmetric_graph<vertex_type, W>;

  size_t num_vertices() { return n; }
  size_t num_edges() { return m; }

  // ======== Graph operators that perform packing ========
  template <class P>
  uintE packNeighbors(uintE id, P& p, uint8_t* tmp) {
    uintE new_degree = get_vertex(id).out_neighbors().pack(p, (std::tuple<uintE, W>*)tmp);
//    std::cout << "Packed id = " << id << " old_degree = " << v_data[id].degree << " new_degree = " << new_degree << std::endl;
    v_data[id].degree = new_degree;  // updates the degree
    return new_degree;
  }

  // degree must be <= old_degree
  void decreaseVertexDegree(uintE id, uintE degree) {
    assert(degree <= v_data[id].degree);
    v_data[id].degree = degree;
  }

  void zeroVertexDegree(uintE id) { decreaseVertexDegree(id, 0); }

  sequence<std::tuple<uintE, uintE, W>> edges() {
    using g_edge = std::tuple<uintE, uintE, W>;
    auto degs = sequence<size_t>(
        n, [&](size_t i) { return get_vertex(i).out_degree(); });
    size_t sum_degs = pbbslib::scan_inplace(make_slice(degs));
    assert(sum_degs == m);
    auto edges = sequence<g_edge>(sum_degs);
    parallel_for(0, n, [&](size_t i) {
      size_t k = degs[i];
      auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
       edges[k++] = std::make_tuple(u, v, wgh);
      };
      get_vertex(i).out_neighbors().map(map_f, false);
    }, 1);
    return edges;
  }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true, size_t granularity=1) {
    parallel_for(0, n, [&](size_t i) {
      get_vertex(i).out_neighbors().map(f, parallel_inner_map);
    }, granularity);
  }

  template <class M, class R>
  typename R::T reduceEdges(M map_f, R reduce_f) {
    using T = typename R::T;
    auto D = pbbslib::make_delayed<T>(n, [&] (size_t i) { return
      get_vertex(i).out_neighbors().reduce(map_f, reduce_f); });
    return pbbslib::reduce(D, reduce_f);
  }

  // ======================= Constructors and fields  ========================
  symmetric_graph()
      : v_data(nullptr),
        e0(nullptr),
        e1(nullptr),
        n(0),
        m(0),
        deletion_fn([]() {}) {}

  symmetric_graph(vertex_data* v_data, size_t n, size_t m,
                  std::function<void()> _deletion_fn, edge_type* _e0,
                  edge_type* _e1 = nullptr)
      : v_data(v_data),
        e0(_e0),
        e1(_e1),
        n(n),
        m(m),
        deletion_fn(_deletion_fn) {
    if (_e1 == nullptr) {
      e1 = e0;  // handles NVM case when graph is stored in symmetric memory
    }
  }
  void del() { deletion_fn(); }

  // creates an in-memory copy of the graph.
  graph copy() {
    auto vd = pbbslib::new_array_no_init<vertex_data>(n);
    auto ed = pbbslib::new_array_no_init<edge_type>(m);
    parallel_for(0, n, [&] (size_t i) {
      vd[i] = v_data[i];
    });
    parallel_for(0, m, [&] (size_t i) {
      ed[i] = e0[i];
    });
    return graph(vd, n, m, [=] () {
      pbbslib::free_array(vd, n);
      pbbslib::free_array(ed, m);
    }, ed);
  }

#ifndef SAGE
  vertex get_vertex(uintE i) { return vertex(e0, v_data[i], i); }
#else
  vertex get_vertex(uintE i) {
    if (pbbslib::numanode() == 0) {
      return vertex(e0, v_data[i], i);
    } else {
      return vertex(e1, v_data[i], i);
    }
  }
#endif

  // Graph Data
  vertex_data* v_data;
  // Pointer to edges
  edge_type* e0;
  // Pointer to second copy of edges--relevant if using 2-socket NVM
  edge_type* e1;

  // number of vertices in G
  size_t n;
  // number of edges in G
  size_t m;

  // called to delete the graph
  std::function<void()> deletion_fn;
};

// Similar to symmetric_graph, but edges are not necessarily allocated
// consecutively. The structure simply stores an array of vertex
// objects (which store an 8-byte pointer, and a uintE degree each).
template <template <class W> class vertex_type, class W>
struct symmetric_ptr_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;
  using graph = symmetric_ptr_graph<vertex_type, W>;

  size_t num_vertices() { return n; }
  size_t num_edges() { return m; }

  // ======== Graph operators that perform packing ========
  template <class P>
  uintE packNeighbors(uintE id, P& p, uint8_t* tmp) {
    uintE new_degree = get_vertex(id).out_neighbors().pack(p, (std::tuple<uintE, W>*)tmp);
    vertices[id].degree = new_degree;  // updates the degree
    return new_degree;
  }

  // degree must be <= old_degree
  void decreaseVertexDegree(uintE id, uintE degree) {
    assert(degree <= vertices[id].degree);
    vertices[id].degree = degree;
  }

  void zeroVertexDegree(uintE id) { decreaseVertexDegree(id, 0); }

  sequence<std::tuple<uintE, uintE, W>> edges() {
    using g_edge = std::tuple<uintE, uintE, W>;
    auto degs = sequence<size_t>(
        n, [&](size_t i) { return get_vertex(i).out_degree(); });
    size_t sum_degs = pbbslib::scan_inplace(make_slice(degs));
    assert(sum_degs == m);
    auto edges = sequence<g_edge>(sum_degs);
    parallel_for(0, n, [&](size_t i) {
      size_t k = degs[i];
      auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
       edges[k++] = std::make_tuple(u, v, wgh);
      };
      get_vertex(i).out_neighbors().map(i, map_f, false);
    }, 1);
    return edges;
  }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) {
    parallel_for(0, n, [&](size_t i) {
      get_vertex(i).out_neighbors().map(i, f, parallel_inner_map);
    }, 1);
  }

  template <class M, class R>
  typename R::T reduceEdges(M map_f, R reduce_f) {
    using T = typename R::T;
    auto D = pbbslib::make_delayed<T>(n, [&] (size_t i) { return
      get_vertex(i).out_neighbors().reduce(i, map_f, reduce_f); });
    return pbbslib::reduce(D, reduce_f);
  }

  // ======================= Constructors and fields  ========================
  symmetric_ptr_graph()
      : n(0),
        m(0),
        vertices(nullptr),
        edge_list_sizes(nullptr),
        deletion_fn([]() {}) {}

  symmetric_ptr_graph(size_t n, size_t m, vertex* _vertices, std::function<void()> _deletion_fn, uintE* _edge_list_sizes=nullptr)
      : n(n),
        m(m),
        vertices(_vertices),
        edge_list_sizes(_edge_list_sizes),
        deletion_fn(_deletion_fn) {
  }
  void del() { deletion_fn(); }

  // creates an in-memory copy of the graph.
  graph copy() {
    vertex* V = pbbslib::new_array_no_init<vertex>(n);
    auto offsets = sequence<size_t>(n+1);
    parallel_for(0, n, [&] (size_t i) {
      V[i] = vertices[i];
      offsets[i] = (edge_list_sizes == nullptr) ? V[i].out_degree() : edge_list_sizes[i];
    });
    offsets[n] = 0;
    size_t total_space = pbbslib::scan_inplace(make_slice(offsets));
    edge_type* E = pbbslib::new_array_no_init<edge_type>(total_space);
    parallel_for(0, n, [&] (size_t i) {
      size_t offset = offsets[i];
      if constexpr (std::is_same<vertex, symmetric_vertex<W>>::value) {
        auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh, size_t ind) {
          E[offset + ind] = std::make_tuple(v, wgh);
        };
        V[i].out_neighbors().map_with_index(map_f);
      } else {
        size_t next_offset = offsets[i+1];
        size_t to_copy = next_offset - offset;
        // memcpy?
        parallel_for(0, to_copy, [&] (size_t j) {
          E[offset + j] = V[i].neighbors[j];
        });
      }
    });
    return graph(n, m, V, [=] () {
        pbbslib::free_array(V, n);
        pbbslib::free_array(E, total_space);
    });
  }

  // Note that observers recieve a handle to a vertex object which is only valid
  // so long as this graph's memory is valid (i.e., before del() has been
  // called).
  vertex get_vertex(uintE i) { return vertices[i]; }

  // number of vertices in G
  size_t n;
  // number of edges in G
  size_t m;
  // pointer to array of vertex objects
  vertex* vertices;
  // pointer to array of vertex edge-list sizes---necessary if copying a
  // compressed graph in this representation.
  uintE* edge_list_sizes;

  // called to delete the graph
  std::function<void()> deletion_fn;
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

#ifndef SAGE
  vertex get_vertex(size_t i) {
    return vertex(out_edges_0, v_out_data[i], in_edges_0, v_in_data[i], i);
  }
#else
  vertex get_vertex(size_t i) {
    if (pbbslib::numanode() == 0) {
      return vertex(out_edges_0, v_out_data[i], in_edges_0, v_in_data[i], i);
    } else {
      return vertex(out_edges_1, v_out_data[i], in_edges_1, v_in_data[i], i);
    }
  }
#endif

  asymmetric_graph()
      : n(0),
        m(0),
        deletion_fn([]() {}),
        v_out_data(nullptr),
        v_in_data(nullptr),
        out_edges_0(nullptr),
        out_edges_1(nullptr),
        in_edges_0(nullptr),
        in_edges_1(nullptr) {}

  asymmetric_graph(vertex_data* v_out_data, vertex_data* v_in_data, size_t n,
                   size_t m, std::function<void()> _deletion_fn,
                   edge_type* _out_edges_0, edge_type* _in_edges_0,
                   edge_type* _out_edges_1 = nullptr,
                   edge_type* _in_edges_1 = nullptr)
      : n(n),
        m(m),
        deletion_fn(_deletion_fn),
        v_out_data(v_out_data),
        v_in_data(v_in_data),
        out_edges_0(_out_edges_0),
        out_edges_1(_out_edges_0),
        in_edges_0(_in_edges_0),
        in_edges_1(_in_edges_0) {
    if (_out_edges_1) {
      out_edges_1 = _out_edges_1;
    }
    if (_in_edges_1) {
      in_edges_1 = _in_edges_1;
    }
  }

  void del() { deletion_fn(); }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) {
    parallel_for(
        0, n,
        [&](size_t i) { get_vertex(i).out_neighbors().map(f, parallel_inner_map); },
        1);
  }
};

// Similar to asymmetric_graph, but edges are not necessarily allocated
// consecutively.
template <template <class W> class vertex_type, class W>
struct asymmetric_ptr_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using edge_type = typename vertex::edge_type;

  // number of vertices in G
  size_t n;
  // number of edges in G
  size_t m;
  // pointer to array of vertex object
  vertex* vertices;

  // called to delete the graph
  std::function<void()> deletion_fn;

  vertex get_vertex(size_t i) {
    return vertices[i];
  }

  asymmetric_ptr_graph()
      : n(0),
        m(0),
        vertices(nullptr),
        deletion_fn([]() {}) {}

  asymmetric_ptr_graph(size_t n, size_t m, vertex* _vertices, std::function<void()> _deletion_fn)
      : n(n),
        m(m),
        vertices(_vertices),
        deletion_fn(_deletion_fn) {}

  void del() { deletion_fn(); }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) {
    parallel_for(
        0, n,
        [&](size_t i) { get_vertex(i).out_neighbors().map(f, parallel_inner_map); },
        1);
  }
};

// Mutates (sorts) the underlying array A containing a black-box description of
// an edge of typename A::value_type. The caller provides functions GetU, GetV,
// and GetW which extract the u, v, and weight of a (u,v,w) edge respective (if
// the edge is a std::tuple<uinte, uintE, W> this is just get<0>, ..<1>, ..<2>
// respectively.
// e.g.:
//   using edge = std::tuple<uintE, uintE, W>;
//   auto get_u = [&] (const edge& e) { return std::get<0>(e); };
//   auto get_v = [&] (const edge& e) { return std::get<1>(e); };
//   auto get_w = [&] (const edge& e) { return std::get<2>(e); };
//   auto G = sym_graph_from_edges<W>(coo1, get_u, get_v, get_w, 10, false);
template <class Wgh, class EdgeSeq, class GetU, class GetV, class GetW>
static inline symmetric_graph<symmetric_vertex, Wgh> sym_graph_from_edges(
    EdgeSeq& A, size_t n, GetU&& get_u, GetV&& get_v, GetW&& get_w,
    bool is_sorted = false) {
  using vertex = symmetric_vertex<Wgh>;
  using edge_type = typename vertex::edge_type;
  size_t m = A.size();

  if (m == 0) {
    if (n == 0) {
      std::function<void()> del = []() {};
      return symmetric_graph<symmetric_vertex, Wgh>(nullptr, 0, 0, del,
                                                    nullptr);
    } else {
      auto v_data = pbbslib::new_array_no_init<vertex_data>(n);
      parallel_for(0, n, [&](size_t i) {
        v_data[i].offset = 0;
        v_data[i].degree = 0;
      });
      return symmetric_graph<symmetric_vertex, Wgh>(
          v_data, n, 0, [=]() { pbbslib::free_array(v_data, n); }, nullptr);
    }
  }

  if (!is_sorted) {
    pbbslib::integer_sort_inplace(make_slice(A), get_u);
  }

  auto starts = sequence<uintT>(n + 1, (uintT)0);

  using neighbor = std::tuple<uintE, Wgh>;
  auto edges = pbbslib::new_array_no_init<neighbor>(m);
  parallel_for(0, m, [&](size_t i) {
    if (i == 0 || (get_u(A[i]) != get_u(A[i - 1]))) {
      starts[get_u(A[i])] = i;
    }
    if (i != (m - 1)) {
      uintE our_vtx = get_u(A[i]);
      uintE next_vtx = get_u(A[i + 1]);
      if (our_vtx != next_vtx && (our_vtx + 1 != next_vtx)) {
        parallel_for(our_vtx + 1, next_vtx,
                [&](size_t k) { starts[k] = i + 1; });
      }
    }
    if (i == (m - 1)) { /* last edge */
      parallel_for(get_u(A[i]) + 1, starts.size(), [&](size_t j) { starts[j] = m; });
    }
    edges[i] = std::make_tuple(get_v(A[i]), get_w(A[i]));
  });

  auto v_data = pbbslib::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    uintT o = starts[i];
    v_data[i].offset = o;
    v_data[i].degree = (uintE)(((i == (n - 1)) ? m : starts[i + 1]) - o);
  });
  return symmetric_graph<symmetric_vertex, Wgh>(
      v_data, n, m, [=]() { pbbslib::free_array(v_data, n); pbbslib::free_array(edges, m); },
      (edge_type*)edges);
}

template <class Wgh>
static inline symmetric_graph<symmetric_vertex, Wgh> sym_graph_from_edges(
    sequence<std::tuple<uintE, uintE, Wgh>>& A, size_t n,
    bool is_sorted = false) {
  using edge = std::tuple<uintE, uintE, Wgh>;
  auto get_u = [&](const edge& e) { return std::get<0>(e); };
  auto get_v = [&](const edge& e) { return std::get<1>(e); };
  auto get_w = [&](const edge& e) { return std::get<2>(e); };
  return sym_graph_from_edges<Wgh>(A, n, get_u, get_v, get_w, is_sorted);
}



template <class Wgh, class EdgeSeq, class GetU, class GetV, class GetW>
std::tuple<uintE, Wgh>* get_edges(EdgeSeq& A, sequence<uintT>& starts, size_t m, const GetU& get_u, const GetV& get_v, const GetW& get_w) {
  using neighbor = std::tuple<uintE, Wgh>;
  auto edges = pbbslib::new_array_no_init<neighbor>(m);
  parallel_for(0, m, [&](size_t i) {
    if (i == 0 || (get_u(A[i]) != get_u(A[i - 1]))) {
      starts[get_u(A[i])] = i;
    }
    if (i != (m - 1)) {
      uintE our_vtx = get_u(A[i]);
      uintE next_vtx = get_u(A[i + 1]);
      if (our_vtx != next_vtx && (our_vtx + 1 != next_vtx)) {
        parallel_for(our_vtx + 1, next_vtx,
                [&](size_t k) { starts[k] = i + 1; });
      }
    }
    if (i == (m - 1)) { /* last edge */
      parallel_for(get_u(A[i]) + 1, starts.size(), [&](size_t j) { starts[j] = m; });
    }
    edges[i] = std::make_tuple(get_v(A[i]), get_w(A[i]));
  });
  return edges;
}


// Mutates (sorts) the underlying array A containing a black-box description of
// an edge of typename A::value_type. The caller provides functions GetU, GetV,
// and GetW which extract the u, v, and weight of a (u,v,w) edge respective (if
// the edge is a std::tuple<uinte, uintE, W> this is just get<0>, ..<1>, ..<2>
// respectively.
// e.g.:
//   using edge = std::tuple<uintE, uintE, W>;
//   auto get_u = [&] (const edge& e) { return std::get<0>(e); };
//   auto get_v = [&] (const edge& e) { return std::get<1>(e); };
//   auto get_w = [&] (const edge& e) { return std::get<2>(e); };
//   auto G = asym_graph_from_edges<W>(coo1, get_u, get_v, get_w, 10, false);
template <class Wgh, class EdgeSeq, class GetU, class GetV, class GetW>
static inline asymmetric_graph<asymmetric_vertex, Wgh> asym_graph_from_edges(
    EdgeSeq& A, size_t n, GetU&& get_u, GetV&& get_v, GetW&& get_w,
    bool is_sorted = false) {
  using vertex = asymmetric_vertex<Wgh>;
  using edge_type = typename vertex::edge_type;
  size_t m = A.size();

  if (m == 0) {
    if (n == 0) {
      std::function<void()> del = []() {};
      return asymmetric_graph<asymmetric_vertex, Wgh>(nullptr, nullptr, 0, 0, del,
                                                     nullptr, nullptr);
    } else {
      auto v_in_data = pbbslib::new_array_no_init<vertex_data>(n);
      auto v_out_data = pbbslib::new_array_no_init<vertex_data>(n);
      parallel_for(0, n, [&](size_t i) {
        v_in_data[i].offset = 0;
        v_in_data[i].degree = 0;

        v_out_data[i].offset = 0;
        v_out_data[i].degree = 0;
      });
      return asymmetric_graph<asymmetric_vertex, Wgh>(
          v_out_data, v_in_data, n, 0, [=]() { pbbslib::free_array(v_out_data, n); pbbslib::free_array(v_in_data, n); }, nullptr, nullptr);
    }
  }

  // flip to create the in-edges
  auto I = sequence<typename EdgeSeq::value_type>(A.size(), [&] (size_t i) {
    using T = typename EdgeSeq::value_type;
    auto e = A[i];
    return T(get_v(e), get_u(e), get_w(e));
  });

  if (!is_sorted) {
    size_t bits = pbbslib::log2_up(n);
    pbbslib::integer_sort_inplace(make_slice(A), get_u);
    pbbslib::integer_sort_inplace(make_slice(I), get_u);
  }

  auto in_starts = sequence<uintT>(n + 1, (uintT)0);
  auto out_starts = sequence<uintT>(n + 1, (uintT)0);

  auto in_edges = get_edges<Wgh>(I, in_starts, m, get_u, get_v, get_w);
  auto out_edges = get_edges<Wgh>(A, out_starts, m, get_u, get_v, get_w);

  auto in_v_data = pbbslib::new_array_no_init<vertex_data>(n);
  auto out_v_data = pbbslib::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    uintT in_o = in_starts[i];
    in_v_data[i].offset = in_o;
    in_v_data[i].degree = (uintE)(((i == (n - 1)) ? m : in_starts[i + 1]) - in_o);

    uintT out_o = out_starts[i];
    out_v_data[i].offset = out_o;
    out_v_data[i].degree = (uintE)(((i == (n - 1)) ? m : out_starts[i + 1]) - out_o);
  }, 1024);
  return asymmetric_graph<asymmetric_vertex, Wgh>(
      out_v_data, in_v_data, n, m, [=]() {
        pbbslib::free_array(in_v_data, n);
        pbbslib::free_array(out_v_data, n);
        pbbslib::free_array(in_edges, m);
        pbbslib::free_array(out_edges, m); },
      (edge_type*)out_edges, (edge_type*)in_edges);
}

template <class Wgh>
static inline asymmetric_graph<asymmetric_vertex, Wgh> asym_graph_from_edges(
    sequence<std::tuple<uintE, uintE, Wgh>>& A, size_t n,
    bool is_sorted = false) {
  using edge = std::tuple<uintE, uintE, Wgh>;
  auto get_u = [&](const edge& e) { return std::get<0>(e); };
  auto get_v = [&](const edge& e) { return std::get<1>(e); };
  auto get_w = [&](const edge& e) { return std::get<2>(e); };
  return asym_graph_from_edges<Wgh>(A, n, get_u, get_v, get_w, is_sorted);
}


}  // namespace gbbs
