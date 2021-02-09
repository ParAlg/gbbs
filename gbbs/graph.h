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
    v_data[id].degree = new_degree;  // updates the degree
    return new_degree;
  }

  // degree must be <= old_degree
  void decreaseVertexDegree(uintE id, uintE degree) {
    assert(degree <= v_data[id].degree);
    v_data[id].degree = degree;
  }

  void zeroVertexDegree(uintE id) { decreaseVertexDegree(id, 0); }

  parlay::sequence<std::tuple<uintE, uintE, W>> edges() {
    using g_edge = std::tuple<uintE, uintE, W>;
    auto degs = parlay::sequence<size_t>::from_function(
        n, [&](size_t i) { return get_vertex(i).out_degree(); });
    size_t sum_degs = pbbslib::scan_add_inplace(parlay::make_slice(degs));
    assert(sum_degs == m);
    auto edges = parlay::sequence<g_edge>::uninitialized(sum_degs);
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
    auto D = parlay::delayed_seq<T>(n, [&] (size_t i) { return
      get_vertex(i).out_neighbors().reduce(map_f, reduce_f); });
    return parlay::reduce(D, reduce_f);
  }

  // ======================= Constructors and fields  ========================
  symmetric_graph()
      : v_data(parlay::make_slice((vertex_data*) nullptr, (vertex_data*) nullptr)),
        e0(parlay::make_slice((edge_type*) nullptr, (edge_type*) nullptr)),
        e1(parlay::make_slice((edge_type*) nullptr, (edge_type*) nullptr)),
        n(0),
        m(0),
        deletion_fn([]() {})
        {}

  symmetric_graph(gbbs::slice<vertex_data> v_data, size_t n, size_t m,
                  std::function<void()> _deletion_fn, gbbs::slice<edge_type> _e0,
                  gbbs::slice<edge_type> _e1)
      : v_data(std::move(v_data)),
        e0(std::move(_e0)),
        e1(std::move(_e1)),
        n(n),
        m(m),
        deletion_fn(_deletion_fn) {}

  ~symmetric_graph() {
    deletion_fn();
  }

  // creates an in-memory copy of the graph.
  graph copy() {
    auto vd = gbbs::new_array_no_init<vertex_data>(n);
    auto ed = gbbs::new_array_no_init<edge_type>(m);
    parallel_for(0, n, [&] (size_t i) {
      vd[i] = v_data[i];
    });
    parallel_for(0, m, [&] (size_t i) {
      ed[i] = e0[i];
    });
    return graph(vd, n, m, [=] () {
      gbbs::free_array(vd, n);
      gbbs::free_array(ed, m);
    }, ed);
  }

#ifndef SAGE
  vertex get_vertex(uintE i) { return vertex(e0.begin(), v_data[i], i); }
#else
  vertex get_vertex(uintE i) {
    // TODO: fix numanode in sched
    assert(false); exit(-1);
    // if (pbbs::numanode() == 0) {
    //   return vertex(e0, v_data[i], i);
    // } else {
    //   return vertex(e1, v_data[i], i);
    // }
  }
#endif

  // Graph Data
  gbbs::slice<vertex_data> v_data;
  // Pointer to edges
  gbbs::slice<edge_type> e0;
  // Pointer to second copy of edges--relevant if using 2-socket NVM
  gbbs::slice<edge_type> e1;

  // number of vertices in G
  size_t n;
  // number of edges in G
  size_t m;

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

  gbbs::slice<vertex_data> v_out_data;
  gbbs::slice<vertex_data> v_in_data;

  /* Pointer to out-edges */
  gbbs::slice<edge_type> out_edges_0;
  /* Pointer to second copy of out-edges--relevant if using 2-socket NVM */
  gbbs::slice<edge_type> out_edges_1;

  /* Pointer to in-edges */
  gbbs::slice<edge_type> in_edges_0;
  /* Pointer to second copy of in-edges--relevant if using 2-socket NVM */
  gbbs::slice<edge_type> in_edges_1;

#ifndef SAGE
  vertex get_vertex(size_t i) {
    return vertex(out_edges_0.begin(), v_out_data[i], in_edges_0.begin(), v_in_data[i], i);
  }
#else
  vertex get_vertex(size_t i) {
    // TODO: fix numanode in sched
    assert(false); exit(-1);
    // if (pbbs::numanode() == 0) {
    //   return vertex(out_edges_0, v_out_data[i], in_edges_0, v_in_data[i], i);
    // } else {
    //   return vertex(out_edges_1, v_out_data[i], in_edges_1, v_in_data[i], i);
    // }
  }
#endif

  asymmetric_graph()
      : n(0),
        m(0),
        deletion_fn([]() {}) {}

  asymmetric_graph(gbbs::slice<vertex_data> v_out_data, gbbs::slice<vertex_data> v_in_data, size_t n,
                   size_t m, std::function<void()> _deletion_fn,
                   gbbs::slice<edge_type> _out_edges_0, gbbs::slice<edge_type> _in_edges_0,
                   gbbs::slice<edge_type> _out_edges_1,
                   gbbs::slice<edge_type> _in_edges_1)
      : n(n),
        m(m),
        deletion_fn(_deletion_fn),
        v_out_data(v_out_data),
        v_in_data(v_in_data),
        out_edges_0(_out_edges_0),
        out_edges_1(_out_edges_0),
        in_edges_0(_in_edges_0),
        in_edges_1(_in_edges_0) {}

  ~asymmetric_graph() {
    deletion_fn();
  }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) {
    parallel_for(0, n, [&](size_t i) {
      get_vertex(i).out_neighbors().map(f, parallel_inner_map);
    }, 1);
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
      return symmetric_graph<symmetric_vertex, Wgh>();
    } else {
      auto v_data = gbbs::new_array_no_init<vertex_data>(n);
      parallel_for(0, n, [&](size_t i) {
        v_data[i].offset = 0;
        v_data[i].degree = 0;
      });
      auto null_edges = parlay::make_slice((edge_type*) nullptr, (edge_type*) nullptr);

      return symmetric_graph<symmetric_vertex, Wgh>(
          parlay::make_slice(v_data, v_data + n), n, 0, [=]() { gbbs::free_array(v_data, n); },
          null_edges, null_edges);
    }
  }

  if (!is_sorted) {
    parlay::integer_sort_inplace(parlay::make_slice(A), get_u);
  }

  auto starts = parlay::sequence<uintT>(n + 1, (uintT)0);

//  for (size_t i=0; i<m; i++) {
//    std::cout << get_u(A[i]) << " " << get_v(A[i]) << std::endl;
//  }

  using neighbor = std::tuple<uintE, Wgh>;
  auto edges = gbbs::new_array_no_init<neighbor>(m);
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

  for (size_t i=0; i<starts.size(); i++) {
    std::cout << "Starts " << i << " " << starts[i] << std::endl;
  }
  std::cout << "Edges.size = " << m << std::endl;

  auto v_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    uintT o = starts[i];
    v_data[i].offset = o;
    v_data[i].degree = (uintE)(((i == (n - 1)) ? m : starts[i + 1]) - o);
  });

  auto G = symmetric_graph<symmetric_vertex, Wgh>(
      parlay::make_slice(v_data, v_data + n),
      n, m, [=]() { gbbs::free_array(v_data, n); gbbs::free_array(edges, m); },
      parlay::make_slice(edges, edges + m), parlay::make_slice(edges, edges + m));
  for (size_t i=0; i<G.n; i++) {
    auto map_fn = [&] (const uintE& u, const uintE& v, const Wgh& wgh) {
      std::cout << u << " " << v << std::endl;
    };
    G.get_vertex(i).out_neighbors().map(map_fn, false);
  }
  return G;
}

template <class Wgh>
static inline symmetric_graph<symmetric_vertex, Wgh> sym_graph_from_edges(
    parlay::sequence<std::tuple<uintE, uintE, Wgh>>& A, size_t n,
    bool is_sorted = false) {
  using edge = std::tuple<uintE, uintE, Wgh>;
  auto get_u = [&](const edge& e) { return std::get<0>(e); };
  auto get_v = [&](const edge& e) { return std::get<1>(e); };
  auto get_w = [&](const edge& e) { return std::get<2>(e); };
  return sym_graph_from_edges<Wgh>(A, n, get_u, get_v, get_w, is_sorted);
}


template <class Wgh, class EdgeSeq, class GetU, class GetV, class GetW>
std::tuple<uintE, Wgh>* get_edges(EdgeSeq& A, gbbs::slice<uintT> starts, size_t m,
    const GetU& get_u, const GetV& get_v, const GetW& get_w) {
  using neighbor = std::tuple<uintE, Wgh>;
  auto edges = gbbs::new_array_no_init<neighbor>(m);
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
      return asymmetric_graph<asymmetric_vertex, Wgh>();
    } else {
      auto v_in_data = gbbs::new_array_no_init<vertex_data>(n);
      auto v_out_data = gbbs::new_array_no_init<vertex_data>(n);
      parallel_for(0, n, [&](size_t i) {
        v_in_data[i].offset = 0;
        v_in_data[i].degree = 0;

        v_out_data[i].offset = 0;
        v_out_data[i].degree = 0;
      });
      return asymmetric_graph<asymmetric_vertex, Wgh>(
          v_out_data, v_in_data, n, 0, [=]() { gbbs::free_array(v_out_data, n); gbbs::free_array(v_in_data, n); }, gbbs::slice<edge_type>(), gbbs::slice<edge_type>());
    }
  }

  // flip to create the in-edges
  auto I = parlay::sequence<typename EdgeSeq::value_type>(A.size(), [&] (size_t i) {
    using T = typename EdgeSeq::value_type;
    auto e = A[i];
    return T(get_v(e), get_u(e), get_w(e));
  });

  if (!is_sorted) {
    parlay::integer_sort_inplace(parlay::make_slice(A), get_u);
    parlay::integer_sort_inplace(parlay::make_slice(I), get_u);
  }

  auto in_starts = parlay::sequence<uintT>(n + 1, (uintT)0);
  auto out_starts = parlay::sequence<uintT>(n + 1, (uintT)0);

  auto in_edges = get_edges<Wgh>(I, parlay::make_slice(in_starts), m, get_u, get_v, get_w);
  auto out_edges = get_edges<Wgh>(A, parlay::make_slice(out_starts), m, get_u, get_v, get_w);

  auto in_v_data = gbbs::new_array_no_init<vertex_data>(n);
  auto out_v_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    uintT in_o = in_starts[i];
    in_v_data[i].offset = in_o;
    in_v_data[i].degree = (uintE)(((i == (n - 1)) ? m : in_starts[i + 1]) - in_o);

    uintT out_o = out_starts[i];
    out_v_data[i].offset = out_o;
    out_v_data[i].degree = (uintE)(((i == (n - 1)) ? m : out_starts[i + 1]) - out_o);
  }, 1024);

  return asymmetric_graph<asymmetric_vertex, Wgh>(
      parlay::make_slice(out_v_data, out_v_data + n),
      parlay::make_slice(in_v_data, in_v_data + n),
      n,
      m,
      [=]() {
        gbbs::free_array(in_v_data, n);
        gbbs::free_array(out_v_data, n);
        gbbs::free_array(in_edges, m);
        gbbs::free_array(out_edges, m); },
      parlay::make_slice(out_edges, out_edges + m),
      parlay::make_slice(in_edges, in_edges + m));
}

template <class Wgh>
static inline asymmetric_graph<asymmetric_vertex, Wgh> asym_graph_from_edges(
    parlay::sequence<std::tuple<uintE, uintE, Wgh>>& A, size_t n,
    bool is_sorted = false) {
  using edge = std::tuple<uintE, uintE, Wgh>;
  auto get_u = [&](const edge& e) { return std::get<0>(e); };
  auto get_v = [&](const edge& e) { return std::get<1>(e); };
  auto get_w = [&](const edge& e) { return std::get<2>(e); };
  return asym_graph_from_edges<Wgh>(A, n, get_u, get_v, get_w, is_sorted);
}


}  // namespace gbbs
