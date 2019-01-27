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

#include "histogram.h"
#include "lib/macros.h"
#include "ligra.h"

#include <type_traits>

// edgeMapInduced
// Version of edgeMapSparse that maps over the one-hop frontier and returns it
// as a sparse array, without filtering.
template <class E, class vertex, class VS, class F>
inline vertexSubsetData<E> edgeMapInduced(graph<vertex>& GA, VS& V, F& f,
                                          bool out_ngh = true) {
  vertex* G = GA.V;
  uintT m = V.size();
  V.toSparse();
  auto degrees = sequence<uintT>(m);
  par_for(0, m, pbbs::kSequentialForThreshold, [&] (size_t i) {
    vertex v = G[V.vtx(i)];
    uintE degree = (out_ngh) ? v.getOutDegree() : v.getInDegree();
    degrees[i] = degree;
  });
  long edgeCount = pbbs::scan_add(degrees, degrees);
  if (edgeCount == 0) {
    return vertexSubsetData<E>(GA.n);
  }
  typedef std::tuple<uintE, E> VE;
  VE* edges = pbbs::new_array_no_init<VE>(edgeCount);

  auto gen = [&](const uintE& ngh, const uintE& offset,
                 const Maybe<E>& val = Maybe<E>()) {
    edges[offset] = std::make_tuple(ngh, val.t);
  };

  if (out_ngh) {
    par_for(0, m, 1, [&] (size_t i) { // TODO-granularity
      uintT o = degrees[i];
      auto v = V.vtx(i);
      G[v].template copyOutNgh(v, o, f, gen);
    });
  } else {
    par_for(0, m, 1, [&] (size_t i) { // TODO-granularity
      uintT o = degrees[i];
      auto v = V.vtx(i);
      G[v].template copyInNgh(v, o, f, gen);
    });
  }
  auto vs = vertexSubsetData<E>(GA.n, edgeCount, edges);
  return vs;
}

// edgeMapInducedFilter
// filter: (u, v, w) -> bool
template <class E, template <typename W> class vertex, class W, class VS, class F, class Filter>
inline vertexSubsetData<E> edgeMapInducedFilter(graph<vertex<W>>& GA, VS& V, F& f, Filter& filter_f,
                                                bool out_ngh = true) {
  using wvertex = vertex<W>;
  wvertex* G = GA.V;
  uintT m = V.size();
  V.toSparse();
  auto degrees = sequence<uintT>(m);
  par_for(0, m, pbbs::kSequentialForThreshold, [&] (size_t i) {
    uintE vtx_id = V.vtx(i);
    wvertex v = G[vtx_id];
    uintE degree = (out_ngh) ? v.countOutNgh(vtx_id, filter_f) : v.countInNgh(vtx_id, filter_f);
    degrees[i] = degree;
  });
  long edgeCount = pbbs::scan_add(degrees, degrees);
  if (edgeCount == 0) {
    return vertexSubsetData<E>(GA.n);
  }
  typedef std::tuple<uintE, E> VE;
  VE* edges = pbbs::new_array_no_init<VE>(edgeCount);

//  auto gen = [&](const uintE& ngh, const uintE& offset,
//                 const Maybe<E>& val = Maybe<E>()) {
//    edges[offset] = std::make_tuple(ngh, val.t);
//  };

  if (out_ngh) {
    par_for(0, m, 1, [&] (size_t i) { // TODO-granularity
      uintT o = degrees[i];
      auto v = V.vtx(i);
      auto out_f = [&] (size_t j, const std::tuple<uintE, W>& nw) {
        Maybe<E> val = f(v, std::get<0>(nw), std::get<1>(nw));
        edges[o + j] = std::make_tuple(std::get<0>(nw), val.t);
      };
      G[v].template filterOutNgh(v, filter_f, out_f, nullptr);
    });
  } else {
    exit(0); // unimplemented
    par_for(0, m, 1, [&] (size_t i) { // TODO-granularity
      uintT o = degrees[i];
      auto v = V.vtx(i);
//      G[v].template copyInNgh(v, o, f, gen);
    });
  }
  auto vs = vertexSubsetData<E>(GA.n, edgeCount, edges);
  return vs;
}


template <class V, template <typename W> class vertex, class W>
struct EdgeMap {
  using K = uintE;  // keys are always uintE's (vertex-identifiers)
  using KV = std::tuple<K, V>;
  using w_vertex = vertex<W>;
  graph<w_vertex>& G;
  pbbs::hist_table<K, V> ht;

  EdgeMap(graph<w_vertex>& _G, KV _empty,
          size_t ht_size = std::numeric_limits<size_t>::max())
      : G(_G) {
    if (ht_size == std::numeric_limits<size_t>::max()) {
      ht_size = G.m / 20;
    }
    ht = pbbs::hist_table<K, V>(_empty, ht_size);
  }

  // map_f: (uintE v, uintE ngh) -> M
  // reduce_f: (E, tuple(uintE ngh, E ngh_val)) -> E
  // apply_f: (uintE ngh, E reduced_val) -> O
  template <class O, class M, class Map, class Reduce, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapReduce(VS& vs, Map& map_f, Reduce& reduce_f,
                                           Apply& apply_f,
                                           bool out_ngh = true) {
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }

    auto wrapped_map_f = [&](const uintE& src, const uintE& ngh, const W& e) {
      return map_f(src, ngh);
    };

    auto oneHop =
        edgeMapInduced<M, w_vertex, VS>(G, vs, wrapped_map_f, out_ngh);
    oneHop.toSparse();

    auto elm_f = [&](size_t i) { return oneHop.vtxAndData(i); };
    auto get_elm = make_sequence<std::tuple<K, M> >(
        oneHop.size(), elm_f);
    auto key_f = [&](size_t i) -> uintE { return oneHop.vtx(i); };
    auto get_key = make_sequence<uintE>(
        oneHop.size(), key_f);

    auto q = [&](sequentialHT<K, V>& S, std::tuple<K, M> v) -> void {
      S.template insertF<M>(v, reduce_f);
    };
    auto res = pbbs::histogram_reduce<std::tuple<K, M>, std::tuple<K, O> >(
        get_elm, get_key, oneHop.size(), q, apply_f, ht);
    oneHop.del();
    return vertexSubsetData<O>(vs.n, res.first, res.second);
  }


  template <class O, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapCount(VS& vs, Apply& apply_f,
                                          bool out_ngh = true) {
    auto map_f = [](const uintE& i, const uintE& j, const W& wgh) {
      return pbbs::empty();
    };
    auto reduce_f = [&](const uintE& cur,
                        const std::tuple<uintE, pbbs::empty>& r) {
      return cur + 1;
    };
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    auto oneHop =
        edgeMapInduced<pbbs::empty, w_vertex, VS>(G, vs, map_f, out_ngh);
    oneHop.toSparse();

    auto key_f = [&](size_t i) -> uintE { return oneHop.vtx(i); };
    auto get_key = make_sequence<uintE>(
        oneHop.size(), key_f);
    auto res = pbbs::histogram<std::tuple<uintE, O> >(get_key, oneHop.size(),
                                                      apply_f, ht);
    oneHop.del();
    return vertexSubsetData<O>(vs.n, res.first, res.second);
  }


  template <class O, class Filter, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapCountFilter(VS& vs, Filter& filter_f, Apply& apply_f,
                                          bool out_ngh = true) {
    auto map_f = [](const uintE& i, const uintE& j, const W& wgh) {
      return pbbs::empty();
    };
    auto reduce_f = [&](const uintE& cur,
                        const std::tuple<uintE, pbbs::empty>& r) {
      return cur + 1;
    };
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    auto oneHop =
        edgeMapInducedFilter<pbbs::empty, vertex, W, VS>(G, vs, map_f, filter_f, out_ngh);
    oneHop.toSparse();


    auto key_f = [&](size_t i) -> uintE { return oneHop.vtx(i); };
    auto get_key = make_sequence<uintE>(
        oneHop.size(), key_f);
    auto res = pbbs::histogram<std::tuple<uintE, O> >(get_key, oneHop.size(),
                                                      apply_f, ht);
    oneHop.del();
    return vertexSubsetData<O>(vs.n, res.first, res.second);
  }

  ~EdgeMap() { ht.del(); }
};
