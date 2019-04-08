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

#include "bridge.h"
#include "histogram.h"
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
  par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    vertex v = G[V.vtx(i)];
    uintE degree = (out_ngh) ? v.getOutDegree() : v.getInDegree();
    degrees[i] = degree;
  });
  long edgeCount = pbbslib::scan_add_inplace(degrees);
  if (edgeCount == 0) {
    return vertexSubsetData<E>(GA.n);
  }
  typedef std::tuple<uintE, E> VE;
  VE* edges = pbbslib::new_array_no_init<VE>(edgeCount);

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

template <class ident_t, class V>
struct HistogramWrapper {
  using K = ident_t;
  using KV = std::tuple<K, V>;

  pbbslib::hist_table<K, V> ht;

  HistogramWrapper(size_t ht_size, KV empty) {
    ht = pbbslib::hist_table<K, V>(empty, ht_size);
  }

  // Wrapper for calling histogram with EM's hash table.
  // apply takes KV -> V
  template <class Apply>
  inline std::pair<size_t, KV*> edgeMapCount(sequence<ident_t>& values_seq, Apply&
      apply_f) {
    return pbbslib::histogram<std::tuple<ident_t, V>>(values_seq, values_seq.size(), apply_f, ht);
  }
};

template <class V, template <typename W> class vertex, class W>
struct EdgeMap {
  using K = uintE;  // keys are always uintE's (vertex-identifiers)
  using KV = std::tuple<K, V>;
  using w_vertex = vertex<W>;
  graph<w_vertex>& G;
  pbbslib::hist_table<K, V> ht;

  EdgeMap(graph<w_vertex>& _G, KV _empty,
          size_t ht_size = std::numeric_limits<size_t>::max())
      : G(_G) {
    if (ht_size == std::numeric_limits<size_t>::max()) {
      ht_size = G.m / 20;
    }
    ht = pbbslib::hist_table<K, V>(_empty, ht_size);
  }

  // sparse [write out neighbors]
  // id: M
  // map_f: (uintE v, uintE ngh) -> M
  // reduce_f: (M, tuple(uintE ngh, M ngh_val)) -> M
  // apply_f: (uintE ngh, M reduced_val) -> Maybe<O>
  template <class O, class M, class Map, class Reduce, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapReduce_sparse(VS& vs, Map& map_f, Reduce& reduce_f, Apply& apply_f, M id, bool out_ngh = true) {
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
    auto get_elm = pbbslib::make_sequence<std::tuple<K, M> >(
        oneHop.size(), elm_f);
    auto key_f = [&](size_t i) -> uintE { return oneHop.vtx(i); };
    auto get_key = pbbslib::make_sequence<uintE>(
        oneHop.size(), key_f);

    auto q = [&](sequentialHT<K, V>& S, std::tuple<K, M> v) -> void {
      S.template insertF<M>(v, reduce_f);
    };
    auto res = pbbslib::histogram_reduce<std::tuple<K, M>, std::tuple<K, O> >(
        get_elm, get_key, oneHop.size(), q, apply_f, ht);
    oneHop.del();
    return vertexSubsetData<O>(vs.n, res.first, res.second);
  }

  // dense [read all neighbors]
  // id: M
  // map_f: (uintE v, uintE ngh) -> M
  // reduce_f: (M, tuple(uintE ngh, M ngh_val)) -> M
  // apply_f: (uintE ngh, M reduced_val) -> Maybe<O>
  template <class O, class M, class Map, class Reduce, class Reduce2, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapReduce_dense(VS& vs, Map& map_f, Reduce& reduce_f, Reduce2& reduce_2, Apply& apply_f, M id, bool out_ngh = true) {
    size_t n = G.n;
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    using OT = std::tuple<bool, O>;
    auto out = pbbslib::new_array<OT>(n);
    auto mf = [&] (uintE u, uintE v, W& wgh) {
      return map_f(u, v);
    };

    par_for(0, n, [&] (size_t i) {
      M reduced_val = (out_ngh) ?
        G.V[i].reduceOutNgh(i, id, mf, reduce_2) :
        G.V[i].reduceInNgh(i, id, mf, reduce_2);
        auto tup = std::make_tuple(i, reduced_val);
      auto applied_val = apply_f(tup);
      if (applied_val.exists) {
        std::get<0>(out[i]) = true;
        std::get<1>(out[i]) = std::get<1>(applied_val.t);
      } else {
        std::get<0>(out[i]) = false;
      }
    });

    return vertexSubsetData<O>(n, out);
  }

  // id: M
  // map_f: (uintE v, uintE ngh) -> M
  // reduce_f: (M, tuple(uintE ngh, M ngh_val)) -> M
  // apply_f: (uintE ngh, M reduced_val) -> Maybe<O>
  template <class O, class M, class Map, class Reduce, class Reduce2, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapReduce(VS& vs, Map& map_f, Reduce& reduce_f, Reduce2& reduce_2, Apply& apply_f, M id, bool out_ngh = true, intT threshold = -1) {
    size_t n = G.n;
    vs.toSparse();
    auto degree_f = [&](size_t i) {
      return (!out_ngh) ? G.V[vs.vtx(i)].getInVirtualDegree()
                       : G.V[vs.vtx(i)].getOutVirtualDegree();
    };
    auto degree_imap = pbbslib::make_sequence<uintE>(vs.size(), degree_f);
    auto out_degrees = pbbslib::reduce_add(degree_imap);
    if (threshold == -1) threshold = G.m / 20;
    if (vs.size() + out_degrees > threshold) {
      // dense
      return edgeMapReduce_dense<O, M>(vs, map_f, reduce_f, reduce_2, apply_f, id, out_ngh);
    } else {
      // sparse
      return edgeMapReduce_sparse<O, M>(vs, map_f, reduce_f, apply_f, id, out_ngh);
    }
  }



  // sparse [write out neighbors]
  template <class O, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapCount_sparse(VS& vs, Apply& apply_f,
                                          bool out_ngh = true) {
    auto map_f = [](const uintE& i, const uintE& j, const W& wgh) {
      return pbbslib::empty();
    };
//    auto reduce_f = [&](const uintE& cur,
//                        const std::tuple<uintE, pbbslib::empty>& r) {
//      return cur + 1;
//    };
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    auto oneHop =
        edgeMapInduced<pbbslib::empty, w_vertex, VS>(G, vs, map_f, out_ngh);
    oneHop.toSparse();

    auto key_f = [&](size_t i) -> uintE { return oneHop.vtx(i); };
    auto get_key = pbbslib::make_sequence<uintE>(oneHop.size(), key_f);
    auto res = pbbslib::histogram<std::tuple<uintE, O> >(get_key, oneHop.size(),
                                                      apply_f, ht);
    oneHop.del();
    return vertexSubsetData<O>(vs.n, res.first, res.second);
  }

  // dense [read all neighbors]
  template <class O, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapCount_dense(VS& vs, Apply& apply_f,
                                                bool out_ngh = true) {
    size_t n = G.n;
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    using OT = std::tuple<bool, O>;
    auto out = pbbslib::new_array<OT>(n);
    vs.toDense();

    cout << "running dense" << endl << endl;

    auto count_f = [&] (uintE u, uintE v, W& wgh) -> size_t {
      return static_cast<size_t>(vs.isIn(v));
    };

    par_for(0, n, 1, [&] (size_t i) {
      size_t count = (out_ngh) ?
        G.V[i].countOutNgh(i, count_f) :
        G.V[i].countInNgh(i, count_f);
      auto tup = std::make_tuple(i, count);
      if (count > 0) {
        auto applied_val = apply_f(tup);
        if (applied_val.exists) {
          std::get<0>(out[i]) = true;
          std::get<1>(out[i]) = std::get<1>(applied_val.t);
        } else {
          std::get<0>(out[i]) = false;
        }
      } else {
        std::get<0>(out[i]) = false;
      }
    });

    return vertexSubsetData<O>(n, out);
  }

  template <class O, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapCount(VS& vs, Apply& apply_f, bool out_ngh = true, long threshold=-1) {
    vs.toSparse();
    auto degree_f = [&](size_t i) -> size_t {
      return (!out_ngh) ? G.V[vs.vtx(i)].getInVirtualDegree()
                       : G.V[vs.vtx(i)].getOutVirtualDegree();
    };
    auto degree_imap = pbbslib::make_sequence<size_t>(vs.size(), degree_f);
    auto out_degrees = pbbslib::reduce_add(degree_imap);
    size_t degree_threshold = threshold;
    if (threshold == -1) degree_threshold = G.m / 15;
    if (vs.size() + out_degrees > degree_threshold) {
      // dense
      return edgeMapCount_dense<O>(vs, apply_f, out_ngh);
    } else {
      // sparse
      return edgeMapCount_sparse<O>(vs, apply_f, out_ngh);
    }
  }

  ~EdgeMap() { ht.del(); }
};
