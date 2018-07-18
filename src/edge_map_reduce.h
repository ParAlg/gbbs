#pragma once

#include "ligra.h"
#include "histogram.h"

#include <type_traits>

// edgeMapInduced
// Version of edgeMapSparse that maps over the one-hop frontier and returns it as
// a sparse array, without filtering.
template <class E, class vertex, class VS, class F>
inline vertexSubsetData<E> edgeMapInduced(graph<vertex>& GA, VS& V, F& f, bool out_ngh=true) {
  vertex *G = GA.V;
  uintT m = V.size();
  V.toSparse();
  auto degrees = array_imap<uintT>(m);
  granular_for(i, 0, m, (m > 2000), {
    vertex v = G[V.vtx(i)];
    uintE degree = (out_ngh) ? v.getOutDegree() : v.getInDegree();
    degrees[i] = degree;
  });
  long edgeCount = pbbs::scan_add(degrees, degrees);
  if (edgeCount == 0) {
    return vertexSubsetData<E>(GA.n);
  }
  typedef tuple<uintE, E> VE;
  VE* edges = pbbs::new_array_no_init<VE>(edgeCount);

  auto gen = [&] (const uintE& ngh, const uintE& offset, const Maybe<E>& val = Maybe<E>()) {
    edges[offset] = make_tuple(ngh, val.t);
  };

  if (out_ngh) {
    granular_for(i, 0, m, (edgeCount > 2000), {
      uintT o = degrees[i];
      auto v = V.vtx(i);
      G[v].template copyOutNgh(v, o, f, gen);
    });
  } else {
    granular_for(i, 0, m, (edgeCount > 2000), {
      uintT o = degrees[i];
      auto v = V.vtx(i);
      G[v].template copyInNgh(v, o, f, gen);
    });
  }
  auto vs = vertexSubsetData<E>(GA.n, edgeCount, edges);
  return vs;
}

template <class W, class Map>
auto wrap_map_f(Map& map_f) {
  return [&] (const uintE& src, const uintE& ngh, const W& e) {
    return map_f(src, ngh);
  };
}

template <class V, template <typename W> class vertex, class W>
struct EdgeMap {
  using K = uintE; // keys are always uintE's (vertex-identifiers)
  using KV = tuple<K, V>;
  using w_vertex = vertex<W>;
  graph<w_vertex>& G;
  pbbs::hist_table<K, V> ht;

  EdgeMap(graph<w_vertex>& _G, KV _empty, size_t ht_size=numeric_limits<size_t>::max()) : G(_G) {
    if (ht_size == numeric_limits<size_t>::max()) {
      ht_size = G.m/20;
    }
    ht = pbbs::hist_table<K, V>(_empty, ht_size);
  }

  // map_f: (uintE v, uintE ngh) -> M
  // reduce_f: (E, tuple(uintE ngh, E ngh_val)) -> E
  // apply_f: (uintE ngh, E reduced_val) -> O
  template <class O, class M, class Map, class Reduce, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapReduce(VS& vs, Map& map_f, Reduce& reduce_f, Apply& apply_f, bool out_ngh=true) {
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }

    auto wrapped_map_f = wrap_map_f<W>(map_f);

    auto oneHop = edgeMapInduced<M, w_vertex, VS>(G, vs, wrapped_map_f, out_ngh);
    oneHop.toSparse();

    auto get_elm = make_in_imap<tuple<K, M> >(oneHop.size(), [&] (size_t i) { return oneHop.vtxAndData(i); });
    auto get_key = make_in_imap<uintE>(oneHop.size(), [&] (size_t i) -> uintE { return oneHop.vtx(i); });

    auto q = [&] (sequentialHT<K, V>& S, tuple<K, M> v) -> void { S.template insertF<M>(v, reduce_f); };
    auto res = pbbs::histogram_reduce<tuple<K, M>, tuple<K, O> >(get_elm, get_key, oneHop.size(), q, apply_f, ht);
    oneHop.del();
    return vertexSubsetData<O>(vs.n, res.first, res.second);
  }

  template <class O, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapCount(VS& vs, Apply& apply_f, bool out_ngh=true) {
    auto map_f = [] (const uintE& i, const uintE& j, const W& wgh) { return pbbs::empty(); };
    auto reduce_f = [&] (const uintE& cur, const tuple<uintE, pbbs::empty>& r) { return cur + 1; };
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    auto oneHop = edgeMapInduced<pbbs::empty, w_vertex, VS>(G, vs, map_f, out_ngh);
    oneHop.toSparse();

    auto get_key = make_in_imap<uintE>(oneHop.size(), [&] (size_t i) -> uintE { return oneHop.vtx(i); });
    auto res = pbbs::histogram<tuple<uintE, O> >(get_key, oneHop.size(), apply_f, ht);
    oneHop.del();
    return vertexSubsetData<O>(vs.n, res.first, res.second);
  }

  ~EdgeMap() {
    ht.del();
  }
};
