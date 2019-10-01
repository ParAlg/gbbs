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
#include "flags.h"

#include <type_traits>

// edgeMapInduced
// Version of edgeMapSparse that maps over the one-hop frontier and returns it
// as a sparse array, without filtering.
template <class E, class Graph, class VS, class F>
inline vertexSubsetData<E> edgeMapInduced(Graph& G, VS& V, F& f, const flags fl) {
  uintT m = V.size();
  V.toSparse();
  auto degrees = sequence<uintT>(m);
  par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    auto v = G.get_vertex(V.vtx(i));
    uintE degree = (fl & in_edges) ? v.getInDegree() : v.getOutDegree();
    degrees[i] = degree;
  });
  long edgeCount = pbbslib::scan_add_inplace(degrees);
  if (edgeCount == 0) {
    return vertexSubsetData<E>(G.n);
  }
  typedef std::tuple<uintE, E> VE;
  VE* edges = pbbslib::new_array_no_init<VE>(edgeCount);

  auto gen = [&](const uintE& ngh, const uintE& offset,
                 const Maybe<E>& val = Maybe<E>()) {
    edges[offset] = std::make_tuple(ngh, val.t);
  };

  if (fl & in_edges) {
    par_for(0, m, 1, [&] (size_t i) {
      uintT o = degrees[i];
      auto v = V.vtx(i);
      G.get_vertex(v).template copyInNgh(v, o, f, gen);
    });
  } else {
    par_for(0, m, 1, [&] (size_t i) {
      uintT o = degrees[i];
      auto v = V.vtx(i);
      G.get_vertex(v).template copyOutNgh(v, o, f, gen);
    });
  }
  auto vs = vertexSubsetData<E>(G.n, edgeCount, edges);
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

template <class V, class Graph>
struct EdgeMap {
  using K = uintE;  // keys are always uintE's (vertex-identifiers)
  using KV = std::tuple<K, V>;
  using W = typename Graph::weight_type;
  Graph& G;
  pbbslib::hist_table<K, V> ht;

  EdgeMap(Graph& _G, KV _empty,
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
  inline vertexSubsetData<O> edgeMapReduce_sparse(VS& vs, Map& map_f, Reduce& reduce_f, Apply& apply_f, M id, const flags fl) {
    size_t m = vs.size();
    size_t n = G.n;
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }

    auto oneHop =
        edgeMapInduced<M, Graph, VS>(G, vs, map_f, fl);
    oneHop.toSparse();

    auto elm_f = [&](size_t i) { return oneHop.vtxAndData(i); };
    auto get_elm = pbbslib::make_sequence<std::tuple<K, M> >(oneHop.size(), elm_f);
    auto key_f = [&](size_t i) -> uintE { return oneHop.vtx(i); };
    auto get_key = pbbslib::make_sequence<uintE>(oneHop.size(), key_f);

    auto q = [&](sequentialHT<K, V>& S, std::tuple<K, M> v) -> void {
      S.template insertF<M>(v, reduce_f);
    };
    auto res = pbbslib::histogram_reduce<std::tuple<K, M>, std::tuple<K, O> >(
        get_elm, get_key, oneHop.size(), q, apply_f, ht);
    oneHop.del();
    auto ret = vertexSubsetData<O>(vs.n, res.first, res.second);
    if (fl & no_output) {
      ret.del();
      return vertexSubsetData<O>(n);
    }
    return ret;
  }

  // dense [read all neighbors]
  // id: M
  // map_f: (uintE v, uintE ngh) -> M
  // reduce_f: M * M -> M
  // apply_f: (uintE ngh, M reduced_val) -> Maybe<O>
  template <class O, class M, class Cond, class Map, class Reduce, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapReduce_dense(VS& vs, Cond& cond_f, Map& map_f, Reduce& reduce_f, Apply& apply_f, M id, const flags fl) {
    size_t n = G.n;
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    using OT = std::tuple<bool, O>;

    auto red_monoid = pbbs::make_monoid(reduce_f, id);
    if (fl & no_output) {
      parallel_for(0, n, [&] (size_t i) {
        if (cond_f(i)) {
          M reduced_val = (fl & in_edges) ?
                           G.get_vertex(i).template reduceInNgh<M>(i, map_f, red_monoid) :
                           G.get_vertex(i).template reduceOutNgh<M>(i, map_f, red_monoid);
          auto tup = std::make_tuple(i, reduced_val);
          apply_f(tup);
        }
      }, 1);
      return vertexSubsetData<O>(n);
    } else {
      auto out = pbbslib::new_array<OT>(n);
      parallel_for(0, n, [&] (size_t i) {
        std::get<0>(out[i]) = false;
        if (cond_f(i)) {
          M reduced_val = (fl & in_edges) ?
                           G.get_vertex(i).template reduceInNgh<M>(i, map_f, red_monoid) :
                           G.get_vertex(i).template reduceOutNgh<M>(i, map_f, red_monoid);
          auto tup = std::make_tuple(i, reduced_val);
          auto applied_val = apply_f(tup);
          if (applied_val.exists) {
            std::get<0>(out[i]) = true;
            std::get<1>(out[i]) = std::get<1>(applied_val.t);
          }
        }
      }, 1);
      return vertexSubsetData<O>(n, out);
    }
  }

  // id: M
  // cond_f: uintE -> bool
  // map_f: (uintE v, uintE ngh) -> M
  // reduce_f: M * M -> M
  // apply_f: (uintE ngh, M reduced_val) -> Maybe<O>
  template <class O, class M, class Cond, class Map, class Reduce, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapReduce(VS& vs, Cond& cond_f, Map& map_f, Reduce& reduce_f, Apply& apply_f, M id, intT threshold = -1, const flags fl = 0) {
    assert(false);
    // Currently unimplemented: only using edgeMapReduce_dense. TODO: finish
    // implementation of this method (if needed).
    exit(0);
    size_t n = G.n;
    vs.toSparse();
    auto degree_f = [&](size_t i) {
      return (fl & in_edges) ? G.get_vertex(vs.vtx(i)).getInVirtualDegree() : G.get_vertex(vs.vtx(i)).getOutVirtualDegree();
    };
    auto degree_imap = pbbslib::make_sequence<uintE>(vs.size(), degree_f);
    auto out_degrees = pbbslib::reduce_add(degree_imap);
    if (threshold == -1) threshold = G.m / 20;
    if (vs.size() + out_degrees > threshold) {
      // dense
      return edgeMapReduce_dense<O, M>(vs, cond_f, map_f, reduce_f, apply_f, id, fl);
    } else {
      // sparse
      auto red_ht = [&] (M cur, std::tuple<uintE, M> nxt) -> M {
        return reduce_f(cur, std::get<1>(nxt));
      };
      return edgeMapReduce_sparse<O, M>(vs, map_f, red_ht, apply_f, id, fl);
    }
  }

  // ============================= Edge Map Count ===============================

  // sparse [write out neighbors]
  template <class O, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapCount_sparse(VS& vs, Apply& apply_f, const flags fl = 0) {
    auto map_f = [](const uintE& i, const uintE& j, const W& wgh) {
      return pbbslib::empty();
    };
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    auto oneHop =
        edgeMapInduced<pbbslib::empty, Graph, VS>(G, vs, map_f, fl);
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
  inline vertexSubsetData<O> edgeMapCount_dense(VS& vs, Apply& apply_f, const flags fl = 0) {
    size_t n = G.n;
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    using OT = std::tuple<bool, O>;
    vs.toDense();

    debug(cout << "running dense" << endl << endl;);

    auto count_f = [&] (uintE u, uintE v, W& wgh) -> size_t {
      return static_cast<size_t>(vs.isIn(v));
    };

    if (fl & no_output) {
      parallel_for(0, n, [&] (size_t i) {
        size_t count = (fl & in_edges) ?
          G.get_vertex(i).countInNgh(i, count_f) :
          G.get_vertex(i).countOutNgh(i, count_f);
        auto tup = std::make_tuple(i, count);
        if (count > 0) {
          apply_f(tup);
        }
      }, 1);
      return vertexSubsetData<O>(n);
    } else {
      auto out = pbbslib::new_array<OT>(n);
      parallel_for(0, n, [&] (size_t i) {
        size_t count = (fl & in_edges) ?
          G.get_vertex(i).countInNgh(i, count_f) :
          G.get_vertex(i).countOutNgh(i, count_f);
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
      }, 1);
      return vertexSubsetData<O>(n, out);
    }
  }

  template <class O, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapCount(VS& vs, Apply& apply_f, const flags fl = 0, long threshold=-1) {
    vs.toSparse();
    auto degree_f = [&](size_t i) -> size_t {
      return (fl & in_edges) ? G.get_vertex(vs.vtx(i)).getInVirtualDegree()
                             : G.get_vertex(vs.vtx(i)).getOutVirtualDegree();
    };
    auto degree_imap = pbbslib::make_sequence<size_t>(vs.size(), degree_f);
    auto out_degrees = pbbslib::reduce_add(degree_imap);
    size_t degree_threshold = threshold;
    if (threshold == -1) degree_threshold = G.m / 20;
    if (vs.size() + out_degrees > degree_threshold) {
      // dense
      return edgeMapCount_dense<O>(vs, apply_f, fl);
    } else {
      // sparse
      return edgeMapCount_sparse<O>(vs, apply_f, fl);
    }
  }

  ~EdgeMap() { ht.del(); }
};
