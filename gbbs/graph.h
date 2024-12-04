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
//  2) W: the edge weight template
//  The graph is represented as an array of edges to neighbors of type
//  vertex_type::neighbor_type.
//  For uncompressed vertices, this type is equal to tuple<uintE, W>.
template <template <class W> class vertex_type, class W>
struct symmetric_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using neighbor_type = typename vertex::neighbor_type;
  using graph = symmetric_graph<vertex_type, W>;

  // Vertices can have an optional weight. Currently the type of this weight is
  // hard-coded as a double. We should consider generalizing this in the future
  // if the need arises.
  using vertex_weight_type = double;

  // The type of an edge. Used when building the graph from_edges, or returning
  // the graph as a list of edges.
  using edge = std::tuple<uintE, uintE, W>;

  size_t num_vertices() const { return n; }
  size_t num_edges() const { return m; }

  // ======== Graph operators that perform packing ========
  template <class P>
  uintE packNeighbors(uintE id, P &p, uint8_t *tmp) {
    uintE new_degree =
        get_vertex(id).out_neighbors().pack(p, (std::tuple<uintE, W> *)tmp);
    v_data[id].degree = new_degree;  // updates the degree
    return new_degree;
  }

  // degree must be <= old_degree
  void decreaseVertexDegree(uintE id, uintE degree) {
    assert(degree <= v_data[id].degree);
    v_data[id].degree = degree;
  }

  // Sets the provided vertex's degree to zero.
  void zeroVertexDegree(uintE id) { decreaseVertexDegree(id, 0); }

  // ======== Other useful graph operators ========

  // Apply the map operator f : (uintE * uintE * W) -> void
  // to each edge.
  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true,
                size_t granularity = 1) const {
    parlay::parallel_for(
        0, n,
        [&](size_t i) {
          get_vertex(i).out_neighbors().map(f, parallel_inner_map);
        },
        granularity);
  }

  // TODO(laxmand): update this function to take a parlay:::monoid.
  template <class M, class R>
  typename R::T reduceEdges(M map_f, R reduce_f) const {
    using T = typename R::T;
    auto D = parlay::delayed_seq<T>(n, [&](size_t i) {
      return get_vertex(i).out_neighbors().reduce(map_f, reduce_f);
    });
    return parlay::reduce(D, reduce_f);
  }

  // Returns the edge set of the graph. Each edge (u,v) will be output twice,
  // once as (u,v) and once as (v,u).
  sequence<edge> edges() const {
    using g_edge = std::tuple<uintE, uintE, W>;
    auto degs = sequence<size_t>::from_function(
        n, [&](size_t i) { return get_vertex(i).out_degree(); });
    size_t sum_degs = parlay::scan_inplace(make_slice(degs));
    assert(sum_degs == m);
    auto edges = sequence<g_edge>(sum_degs);
    parlay::parallel_for(
        0, n,
        [&](size_t i) {
          size_t k = degs[i];
          auto map_f = [&](const uintE &u, const uintE &v, const W &wgh) {
            edges[k++] = std::make_tuple(u, v, wgh);
          };
          get_vertex(i).out_neighbors().map(map_f, false);
        },
        1);
    return edges;
  }

  // Builds a symmetric graph from a sequence of edges. The input edges can be
  // asymmetric (this function will handle symmetrizing the edges).
  static symmetric_graph from_edges(
      const sequence<edge> &edges,
      size_t n = std::numeric_limits<size_t>::max()) {
    size_t m = edges.size();
    if (n == std::numeric_limits<size_t>::max()) {
      n = 1 +
          parlay::reduce(
              parlay::delayed_seq<size_t>(edges.size(), [&](size_t i) {
                return std::max(std::get<0>(edges[i]), std::get<1>(edges[i]));
              }));
    }
    if (m == 0) {
      if (n == 0) {
        std::function<void()> del = []() {};
        return symmetric_graph<vertex_type, W>(nullptr, 0, 0, std::move(del),
                                               nullptr);
      } else {
        auto v_data = gbbs::new_array_no_init<vertex_data>(n);
        parallel_for(0, n, [&](size_t i) {
          v_data[i].offset = 0;
          v_data[i].degree = 0;
        });
        return symmetric_graph<vertex_type, W>(
            v_data, n, 0, [=]() { gbbs::free_array(v_data, n); }, nullptr);
      }
    }

    auto symmetric_edges = EdgeUtils<W>::undirect_and_sort(edges);
    auto offsets = EdgeUtils<W>::compute_offsets(n, symmetric_edges);
    auto neighbors =
        EdgeUtils<W>::template get_neighbors<vertex>(symmetric_edges);
    size_t sym_m = symmetric_edges.size();

    auto v_data = gbbs::new_array_no_init<vertex_data>(n);
    parallel_for(0, n, [&](size_t i) {
      v_data[i].offset = offsets[i];
      v_data[i].degree = ((i == n - 1) ? sym_m : offsets[i + 1]) - offsets[i];
    });

    return symmetric_graph<vertex_type, W>(
        v_data, n, sym_m,
        [=]() {
          gbbs::free_array(v_data, n);
          gbbs::free_array(neighbors, sym_m);
        },
        (neighbor_type *)neighbors);
  }

  // ======================= Constructors and fields  ========================
  symmetric_graph()
      : v_data(nullptr),
        e0(nullptr),
        vertex_weights(nullptr),
        n(0),
        m(0),
        deletion_fn([]() {}) {}

  symmetric_graph(vertex_data *v_data, size_t n, size_t m,
                  std::function<void()> &&_deletion_fn, neighbor_type *_e0,
                  vertex_weight_type *_vertex_weights = nullptr)
      : v_data(v_data),
        e0(_e0),
        vertex_weights(_vertex_weights),
        n(n),
        m(m),
        deletion_fn(_deletion_fn) {}

  // Move constructor
  symmetric_graph(symmetric_graph &&other) noexcept {
    n = other.n;
    m = other.m;
    v_data = other.v_data;
    e0 = other.e0;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.v_data = nullptr;
    other.e0 = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
  }

  // Move assignment
  symmetric_graph &operator=(symmetric_graph &&other) noexcept {
    deletion_fn();
    n = other.n;
    m = other.m;
    v_data = other.v_data;
    e0 = other.e0;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.v_data = nullptr;
    other.e0 = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
    return *this;
  }

  // Copy constructor
  symmetric_graph(const symmetric_graph &other) {
    gbbs_debug(std::cout << "Copying symmetric graph." << std::endl;);
    n = other.n;
    m = other.m;
    v_data = gbbs::new_array_no_init<vertex_data>(n);
    e0 = gbbs::new_array_no_init<neighbor_type>(m);
    parallel_for(0, n, [&](size_t i) { v_data[i] = other.v_data[i]; });
    parallel_for(0, m, [&](size_t i) { e0[i] = other.e0[i]; });
    deletion_fn = [=]() {
      gbbs::free_array(v_data, n);
      gbbs::free_array(e0, m);
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, n);
      }
    };
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(n);
      parallel_for(
          0, n, [&](size_t i) { vertex_weights[i] = other.vertex_weights[i]; });
    }
  }

  ~symmetric_graph() { deletion_fn(); }

  vertex get_vertex(uintE i) const { return vertex(e0, v_data[i], i); }

  // Graph Data
  vertex_data *v_data;
  // Pointer to edges
  neighbor_type *e0;
  // Pointer to vertex weights
  vertex_weight_type *vertex_weights;

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
  using neighbor_type = typename vertex::neighbor_type;
  using graph = symmetric_ptr_graph<vertex_type, W>;
  using vertex_weight_type = double;
  using edge = std::tuple<uintE, uintE, W>;

  size_t num_vertices() const { return n; }
  size_t num_edges() const { return m; }

  // ======== Graph operators that perform packing ========
  template <class P>
  uintE packNeighbors(uintE id, P &p, uint8_t *tmp) {
    uintE new_degree =
        get_vertex(id).out_neighbors().pack(p, (std::tuple<uintE, W> *)tmp);
    vertices[id].degree = new_degree;  // updates the degree
    return new_degree;
  }

  // degree must be <= old_degree
  void decreaseVertexDegree(uintE id, uintE degree) {
    assert(degree <= vertices[id].degree);
    vertices[id].degree = degree;
  }

  void zeroVertexDegree(uintE id) { decreaseVertexDegree(id, 0); }

  sequence<std::tuple<uintE, uintE, W>> edges() const {
    using g_edge = std::tuple<uintE, uintE, W>;
    auto degs = sequence<size_t>::from_function(
        n, [&](size_t i) { return get_vertex(i).out_degree(); });
    size_t sum_degs = parlay::scan_inplace(make_slice(degs));
    assert(sum_degs == m);
    auto edges = sequence<g_edge>(sum_degs);
    parallel_for(
        0, n,
        [&](size_t i) {
          size_t k = degs[i];
          auto map_f = [&](const uintE &u, const uintE &v, const W &wgh) {
            edges[k++] = std::make_tuple(u, v, wgh);
          };
          get_vertex(i).out_neighbors().map(map_f, false);
        },
        1);
    return edges;
  }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) const {
    parallel_for(
        0, n,
        [&](size_t i) {
          get_vertex(i).out_neighbors().map(f, parallel_inner_map);
        },
        1);
  }

  template <class M, class R>
  typename R::T reduceEdges(M map_f, R reduce_f) const {
    using T = typename R::T;
    auto D = parlay::delayed_seq<T>(n, [&](size_t i) {
      return get_vertex(i).out_neighbors().reduce(map_f, reduce_f);
    });
    return parlay::reduce(D, reduce_f);
  }

  // Builds a symmetric_ptr_graph from a sequence of edges. The input edges can
  // be asymmetric (this function will handle symmetrizing the edges).
  static symmetric_ptr_graph from_edges(
      const sequence<edge> &edges,
      size_t n = std::numeric_limits<size_t>::max()) {
    size_t m = edges.size();
    if (n == std::numeric_limits<size_t>::max()) {
      n = 1 +
          parlay::reduce(
              parlay::delayed_seq<size_t>(edges.size(), [&](size_t i) {
                return std::max(std::get<0>(edges[i]), std::get<1>(edges[i]));
              }));
    }
    if (m == 0) {
      if (n == 0) {
        std::function<void()> del = []() {};
        return symmetric_ptr_graph<vertex_type, W>(0, 0, nullptr,
                                                   std::move(del));
      } else {
        auto vertices = gbbs::new_array_no_init<vertex>(n);
        parallel_for(0, n, [&](size_t i) {
          vertex_data vdata;
          vdata.offset = 0;
          vdata.degree = 0;
          vertices[i] = vertex(nullptr, vdata, i);
        });
        return symmetric_ptr_graph<vertex_type, W>(
            n, 0, vertices, [=]() { gbbs::free_array(vertices, n); });
      }
    }

    auto symmetric_edges = EdgeUtils<W>::undirect_and_sort(edges);
    auto offsets = EdgeUtils<W>::compute_offsets(n, symmetric_edges);
    auto neighbors =
        EdgeUtils<W>::template get_neighbors<vertex>(symmetric_edges);
    size_t sym_m = symmetric_edges.size();

    auto vertices = gbbs::new_array_no_init<vertex>(n);
    parallel_for(0, n, [&](size_t i) {
      vertex_data v_data;
      v_data.offset = offsets[i];
      v_data.degree = ((i == n - 1) ? sym_m : offsets[i + 1]) - offsets[i];
      vertices[i] = vertex(neighbors, v_data, i);
    });

    return symmetric_ptr_graph<vertex_type, W>(n, sym_m, vertices, [=]() {
      gbbs::free_array(vertices, n);
      gbbs::free_array(neighbors, sym_m);
    });
  }

  // ======================= Constructors and fields  ========================
  symmetric_ptr_graph()
      : n(0),
        m(0),
        vertices(nullptr),
        edge_list_sizes(nullptr),
        vertex_weights(nullptr),
        deletion_fn([]() {}) {}

  symmetric_ptr_graph(size_t n, size_t m, vertex *_vertices,
                      std::function<void()> _deletion_fn,
                      vertex_weight_type *_vertex_weights = nullptr,
                      uintE *_edge_list_sizes = nullptr)
      : n(n),
        m(m),
        vertices(_vertices),
        edge_list_sizes(_edge_list_sizes),
        vertex_weights(_vertex_weights),
        deletion_fn(_deletion_fn) {}

  // Move constructor
  symmetric_ptr_graph(symmetric_ptr_graph &&other) noexcept {
    n = other.n;
    m = other.m;
    vertices = other.vertices;
    edge_list_sizes = other.edge_list_sizes;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.vertices = nullptr;
    other.edge_list_sizes = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
  }

  // Move assignment
  symmetric_ptr_graph &operator=(symmetric_ptr_graph &&other) noexcept {
    deletion_fn();
    n = other.n;
    m = other.m;
    vertices = other.vertices;
    edge_list_sizes = other.edge_list_sizes;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.vertices = nullptr;
    other.edge_list_sizes = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
    return *this;
  }

  // Copy constructor
  symmetric_ptr_graph(const symmetric_ptr_graph &other) {
    n = other.n;
    m = other.m;
    vertices = gbbs::new_array_no_init<vertex>(n);
    auto offsets = sequence<size_t>(n + 1);
    parallel_for(0, n, [&](size_t i) {
      vertices[i] = other.vertices[i];
      offsets[i] = vertices[i].out_degree();
    });
    offsets[n] = 0;
    size_t total_space = parlay::scan_inplace(make_slice(offsets));
    neighbor_type *E = gbbs::new_array_no_init<neighbor_type>(total_space);

    parallel_for(0, n, [&](size_t i) {
      size_t offset = offsets[i];
      auto map_f = [&](const uintE &u, const uintE &v, const W &wgh,
                       size_t ind) {
        E[offset + ind] = std::make_tuple(v, wgh);
      };
      // Copy neighbor data into E.
      vertices[i].out_neighbors().map_with_index(map_f);
      // Update this vertex's pointer to point to E.
      vertices[i].neighbors = E + offset;
    });

    deletion_fn = [=]() {
      gbbs::free_array(vertices, n);
      gbbs::free_array(E, total_space);
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, n);
      }
    };
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(n);
      parallel_for(
          0, n, [&](size_t i) { vertex_weights[i] = other.vertex_weights[i]; });
    }
  }

  ~symmetric_ptr_graph() { deletion_fn(); }

  // Note that observers recieve a handle to a vertex object which is only
  // valid so long as this graph's memory is valid.
  vertex get_vertex(uintE i) const { return vertices[i]; }

  // number of vertices in G
  size_t n;
  // number of edges in G
  size_t m;
  // pointer to array of vertex objects
  vertex *vertices;
  // pointer to array of vertex edge-list sizes---necessary if copying a
  // compressed graph in this representation.
  uintE *edge_list_sizes;
  // pointer to array of vertex weights
  vertex_weight_type *vertex_weights;

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
 * 2) W: the edge weight template
 *
 * The graph is represented as an array of edges of type
 * vertex_type::neighbor_type, which is just a pair<uintE, W>.
 * */
template <template <class W> class vertex_type, class W>
struct asymmetric_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using neighbor_type = typename vertex::neighbor_type;
  using graph = asymmetric_graph<vertex_type, W>;
  using vertex_weight_type = double;
  using edge = std::tuple<uintE, uintE, W>;

  // number of vertices in G
  size_t n;
  // number of edges in G
  size_t m;
  // called to delete the graph
  std::function<void()> deletion_fn;

  vertex_data *v_out_data;
  vertex_data *v_in_data;

  // Pointer to out-edges
  neighbor_type *out_edges;

  // Pointer to in-edges
  neighbor_type *in_edges;

  // Pointer to vertex weights
  vertex_weight_type *vertex_weights;

  vertex get_vertex(size_t i) const {
    return vertex(out_edges, v_out_data[i], in_edges, v_in_data[i], i);
  }

  asymmetric_graph()
      : n(0),
        m(0),
        deletion_fn([]() {}),
        v_out_data(nullptr),
        v_in_data(nullptr),
        out_edges(nullptr),
        in_edges(nullptr),
        vertex_weights(nullptr) {}

  asymmetric_graph(vertex_data *v_out_data, vertex_data *v_in_data, size_t n,
                   size_t m, std::function<void()> _deletion_fn,
                   neighbor_type *_out_edges, neighbor_type *_in_edges,
                   vertex_weight_type *_vertex_weights = nullptr)
      : n(n),
        m(m),
        deletion_fn(_deletion_fn),
        v_out_data(v_out_data),
        v_in_data(v_in_data),
        out_edges(_out_edges),
        in_edges(_in_edges),
        vertex_weights(_vertex_weights) {}

  // Move constructor
  asymmetric_graph(asymmetric_graph &&other) noexcept {
    n = other.n;
    m = other.m;
    v_out_data = other.v_out_data;
    v_in_data = other.v_in_data;
    out_edges = other.out_edges;
    in_edges = other.in_edges;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.v_out_data = nullptr;
    other.v_in_data = nullptr;
    other.out_edges = nullptr;
    other.in_edges = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
  }

  // Move assignment
  asymmetric_graph &operator=(asymmetric_graph &&other) noexcept {
    deletion_fn();
    n = other.n;
    m = other.m;
    v_out_data = other.v_out_data;
    v_in_data = other.v_in_data;
    out_edges = other.out_edges;
    in_edges = other.in_edges;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.v_out_data = nullptr;
    other.v_in_data = nullptr;
    other.out_edges = nullptr;
    other.in_edges = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
    return *this;
  }

  // Copy constructor
  asymmetric_graph(const asymmetric_graph &other) {
    gbbs_debug(std::cout << "Copying asymmetric graph." << std::endl;);
    n = other.n;
    m = other.m;
    v_out_data = gbbs::new_array_no_init<vertex_data>(n);
    v_in_data = gbbs::new_array_no_init<vertex_data>(n);
    out_edges = gbbs::new_array_no_init<neighbor_type>(m);
    in_edges = gbbs::new_array_no_init<neighbor_type>(m);
    parallel_for(0, n, [&](size_t i) { v_out_data[i] = other.v_out_data[i]; });
    parallel_for(0, n, [&](size_t i) { v_in_data[i] = other.v_in_data[i]; });
    parallel_for(0, m, [&](size_t i) { out_edges[i] = other.out_edges[i]; });
    parallel_for(0, m, [&](size_t i) { in_edges[i] = other.in_edges[i]; });
    deletion_fn = [=]() {
      gbbs::free_array(v_out_data, n);
      gbbs::free_array(v_in_data, n);
      gbbs::free_array(out_edges, m);
      gbbs::free_array(in_edges, m);
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, n);
      }
    };
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(n);
      parallel_for(
          0, n, [&](size_t i) { vertex_weights[i] = other.vertex_weights[i]; });
    }
  }

  ~asymmetric_graph() { deletion_fn(); }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) const {
    parallel_for(
        0, n,
        [&](size_t i) {
          get_vertex(i).out_neighbors().map(f, parallel_inner_map);
        },
        1);
  }

  // Builds an asymmetric graph from a sequence of edges. Note that this
  // function will not remove duplicate edges if they exist.
  static asymmetric_graph from_edges(
      const sequence<edge> &edges,
      size_t n = std::numeric_limits<size_t>::max()) {
    size_t m = edges.size();

    if (n == std::numeric_limits<size_t>::max()) {
      n = 1 +
          parlay::reduce(
              parlay::delayed_seq<size_t>(edges.size(), [&](size_t i) {
                return std::max(std::get<0>(edges[i]), std::get<1>(edges[i]));
              }));
    }

    if (m == 0) {
      if (n == 0) {
        std::function<void()> del = []() {};
        return graph(nullptr, nullptr, 0, 0, del, nullptr, nullptr);
      } else {
        auto v_in_data = gbbs::new_array_no_init<vertex_data>(n);
        auto v_out_data = gbbs::new_array_no_init<vertex_data>(n);
        parallel_for(0, n, [&](size_t i) {
          v_in_data[i].offset = 0;
          v_in_data[i].degree = 0;
          v_out_data[i].offset = 0;
          v_out_data[i].degree = 0;
        });
        return graph(
            v_out_data, v_in_data, n, 0,
            [=]() {
              gbbs::free_array(v_out_data, n);
              gbbs::free_array(v_in_data, n);
            },
            nullptr, nullptr);
      }
    }

    auto all_out_edges = EdgeUtils<W>::sort_edges(edges);
    auto all_in_edges = EdgeUtils<W>::transpose(edges);
    EdgeUtils<W>::sort_edges_inplace(all_in_edges);

    auto out_offsets = EdgeUtils<W>::compute_offsets(n, all_out_edges);
    auto in_offsets = EdgeUtils<W>::compute_offsets(n, all_in_edges);

    auto out_edges =
        EdgeUtils<W>::template get_neighbors<vertex>(all_out_edges);
    auto in_edges = EdgeUtils<W>::template get_neighbors<vertex>(all_in_edges);

    auto in_v_data = gbbs::new_array_no_init<vertex_data>(n);
    auto out_v_data = gbbs::new_array_no_init<vertex_data>(n);

    parlay::parallel_for(0, n, [&](size_t i) {
      in_v_data[i].offset = in_offsets[i];
      in_v_data[i].degree =
          (uintE)(((i == (n - 1)) ? m : in_offsets[i + 1]) - in_offsets[i]);

      out_v_data[i].offset = out_offsets[i];
      out_v_data[i].degree =
          (uintE)(((i == (n - 1)) ? m : out_offsets[i + 1]) - out_offsets[i]);
    });
    return asymmetric_graph<asymmetric_vertex, W>(
        out_v_data, in_v_data, n, m,
        [=]() {
          gbbs::free_array(in_v_data, n);
          gbbs::free_array(out_v_data, n);
          gbbs::free_array(in_edges, m);
          gbbs::free_array(out_edges, m);
        },
        (neighbor_type *)out_edges, (neighbor_type *)in_edges);
  }
};

// Similar to asymmetric_graph, but edges are not necessarily allocated
// consecutively.
template <template <class W> class vertex_type, class W>
struct asymmetric_ptr_graph {
  using vertex = vertex_type<W>;
  using weight_type = W;
  using neighbor_type = typename vertex::neighbor_type;
  using vertex_weight_type = double;
  using edge = std::tuple<uintE, uintE, W>;
  using graph = asymmetric_ptr_graph<vertex_type, W>;

  // number of vertices in G
  size_t n;
  // number of edges in G
  size_t m;
  // pointer to array of vertex object
  vertex *vertices;
  // pointer to array of vertex weights
  vertex_weight_type *vertex_weights;

  // called to delete the graph
  std::function<void()> deletion_fn;

  vertex get_vertex(size_t i) const { return vertices[i]; }

  asymmetric_ptr_graph()
      : n(0),
        m(0),
        vertices(nullptr),
        deletion_fn([]() {}),
        vertex_weights(nullptr) {}

  asymmetric_ptr_graph(size_t n, size_t m, vertex *_vertices,
                       std::function<void()> _deletion_fn,
                       vertex_weight_type *_vertex_weights = nullptr)
      : n(n),
        m(m),
        vertices(_vertices),
        vertex_weights(_vertex_weights),
        deletion_fn(_deletion_fn) {}

  // Move constructor
  asymmetric_ptr_graph(asymmetric_ptr_graph &&other) noexcept {
    n = other.n;
    m = other.m;
    vertices = other.vertices;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.vertices = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
  }

  // Move assignment
  asymmetric_ptr_graph &operator=(asymmetric_ptr_graph &&other) noexcept {
    deletion_fn();
    n = other.n;
    m = other.m;
    vertices = other.vertices;
    vertex_weights = other.vertex_weights;
    deletion_fn = std::move(other.deletion_fn);
    other.vertices = nullptr;
    other.vertex_weights = nullptr;
    other.deletion_fn = []() {};
    return *this;
  }

  // Copy constructor
  asymmetric_ptr_graph(const asymmetric_ptr_graph &other) {
    n = other.n;
    m = other.m;
    vertices = gbbs::new_array_no_init<vertex>(n);
    auto in_offsets = sequence<size_t>(n + 1);
    auto out_offsets = sequence<size_t>(n + 1);
    parallel_for(0, n, [&](size_t i) {
      vertices[i] = other.vertices[i];
      out_offsets[i] = vertices[i].out_degree();
      in_offsets[i] = vertices[i].in_degree();
    });
    in_offsets[n] = 0;
    out_offsets[n] = 0;

    size_t in_space = parlay::scan_inplace(make_slice(in_offsets));
    size_t out_space = parlay::scan_inplace(make_slice(out_offsets));
    neighbor_type *inE = gbbs::new_array_no_init<neighbor_type>(in_space);
    neighbor_type *outE = gbbs::new_array_no_init<neighbor_type>(out_space);

    parallel_for(0, n, [&](size_t i) {
      size_t out_offset = out_offsets[i];
      if (out_offsets[i + 1] != out_offset) {
        auto map_f = [&](const uintE &u, const uintE &v, const W &wgh,
                         size_t ind) {
          outE[out_offset + ind] = std::make_tuple(v, wgh);
        };
        // Copy neighbor data into E.
        vertices[i].out_neighbors().map_with_index(map_f);
        // Update this vertex's pointer to point to E.
        vertices[i].out_nghs = outE + out_offset;
      }

      size_t in_offset = in_offsets[i];
      if (in_offsets[i + 1] != in_offset) {
        auto map_f = [&](const uintE &u, const uintE &v, const W &wgh,
                         size_t ind) {
          inE[in_offset + ind] = std::make_tuple(v, wgh);
        };
        // Copy neighbor data into E.
        vertices[i].in_neighbors().map_with_index(map_f);
        // Update this vertex's pointer to point to E.
        vertices[i].in_nghs = inE + in_offset;
      }
    });

    deletion_fn = [=]() {
      gbbs::free_array(vertices, n);
      gbbs::free_array(inE, in_space);
      gbbs::free_array(outE, out_space);
      if (vertex_weights != nullptr) {
        gbbs::free_array(vertex_weights, n);
      }
    };
    vertex_weights = nullptr;
    if (other.vertex_weights != nullptr) {
      vertex_weights = gbbs::new_array_no_init<vertex_weight_type>(n);
      parallel_for(
          0, n, [&](size_t i) { vertex_weights[i] = other.vertex_weights[i]; });
    }
  }

  ~asymmetric_ptr_graph() { deletion_fn(); }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) const {
    parallel_for(
        0, n,
        [&](size_t i) {
          get_vertex(i).out_neighbors().map(f, parallel_inner_map);
        },
        1);
  }

  // Builds an asymmetric_ptr_graph from a sequence of edges. Note that this
  // function will not remove duplicate edges if they exist.
  static asymmetric_ptr_graph from_edges(
      const sequence<edge> &edges,
      size_t n = std::numeric_limits<size_t>::max()) {
    size_t m = edges.size();

    if (n == std::numeric_limits<size_t>::max()) {
      n = 1 +
          parlay::reduce(
              parlay::delayed_seq<size_t>(edges.size(), [&](size_t i) {
                return std::max(std::get<0>(edges[i]), std::get<1>(edges[i]));
              }));
    }

    if (m == 0) {
      if (n == 0) {
        std::function<void()> del = []() {};
        return graph(0, 0, nullptr, del);
      } else {
        auto vertices = gbbs::new_array_no_init<vertex>(n);
        parallel_for(0, n, [&](size_t i) {
          vertex_data data;
          data.offset = 0;
          data.degree = 0;
          vertices[i] = vertex(nullptr, data, nullptr, data, i);
        });
        return graph(n, 0, vertices, [=]() { gbbs::free_array(vertices, n); });
      }
    }

    auto all_out_edges = EdgeUtils<W>::sort_edges(edges);
    auto all_in_edges = EdgeUtils<W>::transpose(edges);
    EdgeUtils<W>::sort_edges_inplace(all_in_edges);

    auto out_offsets = EdgeUtils<W>::compute_offsets(n, all_out_edges);
    auto in_offsets = EdgeUtils<W>::compute_offsets(n, all_in_edges);

    auto out_edges =
        EdgeUtils<W>::template get_neighbors<vertex>(all_out_edges);
    auto in_edges = EdgeUtils<W>::template get_neighbors<vertex>(all_in_edges);

    auto vertices = gbbs::new_array_no_init<vertex>(n);

    parlay::parallel_for(0, n, [&](size_t i) {
      vertex_data in_v_data;
      vertex_data out_v_data;
      in_v_data.offset = in_offsets[i];
      in_v_data.degree =
          (uintE)(((i == (n - 1)) ? m : in_offsets[i + 1]) - in_offsets[i]);
      out_v_data.offset = out_offsets[i];
      out_v_data.degree =
          (uintE)(((i == (n - 1)) ? m : out_offsets[i + 1]) - out_offsets[i]);

      vertices[i] = vertex(out_edges, out_v_data, in_edges, in_v_data, i);
    });
    return graph(n, m, vertices, [=]() {
      gbbs::free_array(vertices, n);
      gbbs::free_array(in_edges, m);
      gbbs::free_array(out_edges, m);
    });
  }
};

}  // namespace gbbs
