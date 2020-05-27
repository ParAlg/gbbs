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
#include "edge_map_data.h"
#include "edge_map_reduce.h"
#include "flags.h"
#include "graph_mutation.h"
#include "macros.h"
#include "vertex.h"
#include "vertex_subset.h"

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

  // =================== VertexSubset Operators: Mapping =====================
  // Applies the the supplied (update, cond) operators on edges incident to the
  // input vertex_subset, aggregating results at the neighbors of this vset.
  // This is the specialized version of the general nghMap function where the
  // output is a plain vertex_subset.
  // F is a struct providing the following functions:
  //   update : (uintE * uintE * W) -> bool
  //   updateAtomic : (uintE * uintE * W) -> bool
  //   cond : uintE -> bool
  template <class VS, class F>
  inline vertexSubsetData<pbbs::empty> nghMap(VS& vs, F f, intT threshold = -1,
                                              flags fl = 0) {
    return edgeMapData<pbbs::empty>(*this, vs, f, threshold, fl);
  }

  // The generalized version of edgeMap. Takes an input vertex_subset and
  // aggregates results at the neighbors of the vset. See above for the
  // requirements on F.
  template <class Data,  // data associated with vertices in the output
            class VS,    // input vertex_subset type
            class F>     // map function type
  inline vertexSubsetData<Data>
  nghMapData(VS& vs, F f, intT threshold = -1, flags fl = 0) {
    return edgeMapData<Data>(*this, vs, f, threshold, fl);
  }

  // srcMap. Applies
  template <class VS,  // input vertex_subset type
            class F>   // map function type
  inline vertexSubsetData<pbbs::empty>
  srcMap(VS& vs, F f, intT threshold = -1, flags fl = 0) {
    std::cout << "srcMap currently not implemented" << std::endl;
    exit(-1);  // currently unused by any benchmark
  }

  // The generalized version of srcMap. Takes an input vertex_subset and
  // aggregates results at the neighbors of the vset.
  template <class Data,  // data associated with vertices in the output
            class VS,    // input vertex_subset type
            class F>     // map function type
  inline vertexSubsetData<Data>
  srcMap(VS& vs, F f, intT threshold = -1, flags fl = 0) {
    std::cout << "srcMap currently not implemented" << std::endl;
    exit(-1);  // currently unused by any benchmark
  }

  // =================== VertexSubset Operators: Counting =====================
  template <
      class Data,  // data associated with vertices in the output vertex_subset
      class Cond,
      class Apply,  // function from std::tuple<uintE, uintE> ->
                    // std::optional<std::tuple<uintE, Data>>
      class VS>
  inline vertexSubsetData<Data> nghCount(VS& vs, Cond cond_f, Apply apply_f,
                                         pbbslib::hist_table<uintE, Data>& ht,
                                         flags fl = 0) {
    static_assert(std::is_same<Data, uintE>::value,
                  "Histogram code used in the implementation is specialized "
                  "for Data == counting_type (in this case uintE) for "
                  "performance.");
    return edgeMapCount<Data, Cond, Apply, VS>(*this, vs, cond_f, apply_f, ht, fl);
  }

  template <class Cond, class VS>
  inline vertexSubsetData<uintE> nghCount(VS& vs, Cond cond_f,
                                          pbbslib::hist_table<uintE, uintE>& ht,
                                          flags fl = 0) {
    auto apply_f = [&](const std::tuple<uintE, uintE>& ct) {
      return std::optional<std::tuple<uintE, uintE>>(ct);
    };
    return edgeMapCount<uintE, Cond, decltype(apply_f), VS>(*this, vs, cond_f,
                                                            apply_f, ht, fl);
  }

  template <class P, class VS>
  inline vertexSubsetData<uintE> srcCount(VS& vs, P p, flags fl = 0) {
    return edgeMapFilter(*this, vs, p, fl);
  }

  // =================== VertexSubset Operators: Reduction =====================
  template <class Data, class Apply, class Map, class Reduce, class VS>
  inline vertexSubsetData<Data> nghReduce(VS& vs,
                                          pbbslib::hist_table<uintE, Data>& ht,
                                          Apply apply_f, Map map_f,
                                          Reduce reduce_f, flags fl = 0) {
    std::cout
        << "TODO: map edgeMapReduce interface to the one in edge_map_reduce.h"
        << std::endl;
    exit(-1);
    return vertexSubset(n);
  }

  template <class VS, class Map, class Reduce>
  inline vertexSubsetData<uintE> nghReduce(VS& vs, Map map_f, Reduce reduce_f,
                                           flags fl = 0) {
    std::cout
        << "TODO: map edgeMapReduce interface to the one in edge_map_reduce.h"
        << std::endl;
    exit(-1);
    return vertexSubset(n);
  }

  // =================== VertexSubset Operators: Packing =====================
  template <class P>
  vertexSubsetData<uintE> srcPack(vertexSubset& vs, P p, flags fl = 0) {
    return packEdges(*this, vs, p, fl);
  }

  // Similar to srcPack, but returns the nghs as a vs.
  template <class P>
  vertexSubsetData<uintE> nghPack(vertexSubset& vs, P p, flags fl = 0) {
    assert(false);  // not implemented as this primitive is unused in the
                    // benchmark.
    exit(-1);
    return vertexSubsetData<uintE>();
  }

  // ======================= Graph Operators: Filtering ========================
  // Filters the symmetric graph, G, with a predicate function pred.  Note
  // that the predicate does not have to be symmetric, i.e. f(u,v) is
  // not necesssarily equal to f(v,u), but we only represent the out-edges of
  // this
  // (possibly) directed graph. For convenience in cases where the graph needed
  // is
  // symmetric, we coerce this to a symmetric_graph.
  template <template <class inner_wgh> class vtx_type, class wgh_type,
            typename P,
            typename std::enable_if<
                std::is_same<vtx_type<wgh_type>, symmetric_vertex<W>>::value,
                int>::type = 0>
  static inline symmetric_graph<symmetric_vertex, wgh_type> filterGraph(
      symmetric_graph<vtx_type, wgh_type>& G, P& pred) {
    auto[newN, newM, newVData, newEdges] = filter_graph<vtx_type, W>(G, pred);
    assert(newN == n);
    return symmetric_graph<symmetric_vertex, W>(
        newVData, newN, newM,
        [=]() { pbbslib::free_arrays(newVData, newEdges); }, newEdges);
  }

  template <
      template <class inner_wgh> class vtx_type, class wgh_type, typename P,
      typename std::enable_if<
          std::is_same<vtx_type<wgh_type>, csv_bytepd_amortized<W>>::value,
          int>::type = 0>
  static inline symmetric_graph<csv_byte, wgh_type> filterGraph(
      symmetric_graph<vtx_type, wgh_type>& G, P& pred) {
    auto[newN, newM, newVData, newEdges] = filter_graph<vtx_type, W>(G, pred);
    assert(newN == n);
    return symmetric_graph<csv_byte, W>(
        newVData, newN, newM,
        [=]() { pbbslib::free_arrays(newVData, newEdges); }, newEdges);
  }

  // Used by MST and MaximalMatching
  // Predicate returns three values:
  // 0 : keep in graph, do not return in edge array
  // 1 : remove from graph, do not return in edge array
  // 2 : remove from graph, return in edge array
  // Cost: O(n+m) work
  template <class P>
  edge_array<W> filterEdges(P& pred, flags fl = 0) {
    return filter_edges(*this, pred, fl);
  }

  template <class P>
  edge_array<W> filterAllEdges(P& pred, flags fl = 0) {
    return filter_all_edges(*this, pred, fl);
  }

  template <class P>
  edge_array<W> sampleEdges(P& pred) {
    return sample_edges(*this, pred);
  }

  // ======================= Graph Operators: Packing ========================
  template <class P>
  uintE packNeighbors(uintE id, P& p, std::tuple<uintE, W>* tmp) {
    uintE new_degree = get_vertex(id).packOutNgh(id, p, tmp);
    v_data[id].degree = new_degree;  // updates the degree
    return new_degree;
  }

  // degree must be <= old_degree
  void decreaseVertexDegree(uintE id, uintE degree) {
    assert(degree <= v_data[id].degree);
    v_data[id].degree = degree;
  }

  void zeroVertexDegree(uintE id) { decreaseVertexDegree(id, 0); }

  pbbs::sequence<std::tuple<uintE, uintE, W>> edges() {
    using g_edge = std::tuple<uintE, uintE, W>;
    auto degs = pbbs::sequence<size_t>(
        n, [&](size_t i) { return get_vertex(i).getOutDegree(); });
    size_t sum_degs = pbbslib::scan_add_inplace(degs.slice());
    assert(sum_degs == m);
    auto edges = pbbs::sequence<g_edge>(sum_degs);
    parallel_for(0, n,
                 [&](size_t i) {
                   size_t k = degs[i];
                   auto map_f = [&](const uintE& u, const uintE& v,
                                    const W& wgh) {
                     edges[k++] = std::make_tuple(u, v, wgh);
                   };
                   get_vertex(i).mapOutNgh(i, map_f, false);
                 },
                 1);
    return edges;
  }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) {
    parallel_for(
        0, n,
        [&](size_t i) { get_vertex(i).mapOutNgh(i, f, parallel_inner_map); },
        1);
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

#ifndef SAGE
  vertex get_vertex(uintE i) { return vertex(e0, v_data[i]); }
#else
  vertex get_vertex(uintE i) {
    if (numanode() == 0) {
      return vertex(e0, v_data[i]);
    } else {
      return vertex(e1, v_data[i]);
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

  // Applies the edgeMap operator on the input vertex_subset, aggregating
  // results
  // at the neighbors of this vset. This is the specialized version of the
  // general nghMap function where the output is a plain vertex_subset.
  // F must provide:
  // update : (uintE * uintE * W) -> bool
  // updateAtomic : (uintE * uintE * W) -> bool
  // cond : uintE -> bool
  template <class VS, /* input vertex_subset type */
            class F /* edgeMap function type */>
  inline vertexSubsetData<pbbs::empty> nghMap(VS& vs, F f, intT threshold = -1,
                                              flags fl = 0) {
    return edgeMapData<pbbs::empty>(*this, vs, f, threshold, fl);
  }

  // The generalized version of edgeMap. Takes an input vertex_subset and
  // aggregates results at the neighbors of the vset.
  template <class Data, /* data associated with vertices in the output
                           vertex_subset */
            class VS,   /* input vertex_subset type */
            class F /* edgeMap function type */>
  inline vertexSubsetData<Data> nghMap(VS& vs, F f, intT threshold = -1,
                                       flags fl = 0) {
    return edgeMapData(*this, vs, f, threshold, fl);
  }

  template <class F>
  void mapEdges(F f, bool parallel_inner_map = true) {
    parallel_for(
        0, n,
        [&](size_t i) { get_vertex(i).mapOutNgh(i, f, parallel_inner_map); },
        1);
  }
};

// Mutates (sorts) the underlying array A containing a black-box description
// of
// an edge of typename A::value_type. The caller provides functions GetU,
// GetV,
// and GetW
// which extract the u, v, and weight of a (u,v,w) edge respective (if the
// edge
// is a std::tuple<uinte, uintE, W> this is just get<0>, ..<1>, ..<2>
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
      auto v_data = pbbs::new_array_no_init<vertex_data>(n);
      par_for(0, n, pbbslib::kSequentialForThreshold, [&](size_t i) {
        v_data[i].offset = 0;
        v_data[i].degree = 0;
      });
      return symmetric_graph<symmetric_vertex, Wgh>(
          v_data, n, 0, [=]() { pbbs::free_array(v_data); }, nullptr);
    }
  }

  if (!is_sorted) {
    size_t bits = pbbslib::log2_up(n);
    pbbslib::integer_sort_inplace(A.slice(), get_u, bits);
  }

  auto starts = sequence<uintT>(n + 1, (uintT)0);

  using neighbor = std::tuple<uintE, Wgh>;
  auto edges = sequence<neighbor>(m, [&](size_t i) {
    if (i == 0 || (get_u(A[i]) != get_u(A[i - 1]))) {
      starts[get_u(A[i])] = i;
    }
    if (i != (m - 1)) {
      uintE our_vtx = get_u(A[i]);
      uintE next_vtx = get_u(A[i + 1]);
      if (our_vtx != next_vtx && (our_vtx + 1 != next_vtx)) {
        par_for(our_vtx + 1, next_vtx, pbbslib::kSequentialForThreshold,
                [&](size_t k) { starts[k] = i + 1; });
      }
    }
    if (i == (m - 1)) { /* last edge */
      par_for(get_u(A[i]) + 1, starts.size(), [&](size_t j) { starts[j] = m; });
    }
    return std::make_tuple(get_v(A[i]), get_w(A[i]));
  });

  auto v_data = pbbs::new_array_no_init<vertex_data>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&](size_t i) {
    uintT o = starts[i];
    v_data[i].offset = o;
    v_data[i].degree = (uintE)(((i == (n - 1)) ? m : starts[i + 1]) - o);
  });
  auto new_edge_arr = edges.to_array();
  return symmetric_graph<symmetric_vertex, Wgh>(
      v_data, n, m, [=]() { pbbslib::free_arrays(v_data, new_edge_arr); },
      (edge_type*)new_edge_arr);
}

template <class Wgh>
static inline symmetric_graph<symmetric_vertex, Wgh> sym_graph_from_edges(
    pbbs::sequence<std::tuple<uintE, uintE, Wgh>>& A, size_t n,
    bool is_sorted = false) {
  using edge = std::tuple<uintE, uintE, Wgh>;
  auto get_u = [&](const edge& e) { return std::get<0>(e); };
  auto get_v = [&](const edge& e) { return std::get<1>(e); };
  auto get_w = [&](const edge& e) { return std::get<2>(e); };
  return sym_graph_from_edges<Wgh>(A, n, get_u, get_v, get_w, is_sorted);
}
