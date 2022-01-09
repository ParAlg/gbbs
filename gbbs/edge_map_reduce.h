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
#include "flags.h"
#include "helpers/histogram.h"
#include "vertex_subset.h"

#include <type_traits>

namespace gbbs {

// edgeMapInduced
// Version of edgeMapSparse that maps over the one-hop frontier and returns it
// as a sparse array, without filtering.
template <class E, class Graph, class VS, class Map, class Cond>
inline vertexSubsetData<E> edgeMapInduced(Graph& G, VS& V, Map& map_f,
                                          Cond& cond_f, uintE empty_key,
                                          const flags fl) {
  uintT m = V.size();
  V.toSparse();
  auto degrees = sequence<uintT>(m);
  parallel_for(0, m, [&](size_t i) {
    auto v = G.get_vertex(V.vtx(i));
    uintE degree = (fl & in_edges) ? v.in_degree() : v.out_degree();
    degrees[i] = degree;
  });
  long edgeCount = parlay::scan_inplace(make_slice(degrees));
  if (edgeCount == 0) {
    return vertexSubsetData<E>(G.n);
  }
  using S = typename vertexSubsetData<E>::S;
  auto edges = sequence<S>::uninitialized(edgeCount);

  auto gen = [&](const uintE& ngh, const uintE& offset,
                 const std::optional<E>& val = std::nullopt) {
    if (cond_f(ngh)) {
      if
        constexpr(!std::is_same<E, gbbs::empty>()) {
          edges[offset] = std::make_tuple(ngh, *val);
        }
      else {
        edges[offset] = ngh;
      }
    } else {
      if
        constexpr(!std::is_same<E, gbbs::empty>()) {
          edges[offset] = std::make_tuple(empty_key, *val);
        }
      else {
        edges[offset] = empty_key;
      }
    }
  };

  if (fl & in_edges) {
    parallel_for(0, m,
                 [&](size_t i) {
                   uintT o = degrees[i];
                   auto v = V.vtx(i);
                   G.get_vertex(v).in_neighbors().copy(o, map_f, gen);
                 },
                 1);
  } else {
    parallel_for(0, m,
                 [&](size_t i) {
                   uintT o = degrees[i];
                   auto v = V.vtx(i);
                   G.get_vertex(v).out_neighbors().copy(o, map_f, gen);
                 },
                 1);
  }
  auto vs = vertexSubsetData<E>(G.n, std::move(edges));
  return vs;
}

// ============================= Edge Map Count ===============================
// sparse [write out neighbors]
template <class O, class Cond, class Apply, class VS, class Graph>
inline vertexSubsetData<O> edgeMapCount_sparse(Graph& GA, VS& vs,
                                               hist_table<uintE, O>& ht,
                                               Cond& cond_f, Apply& apply_f,
                                               const flags fl = 0) {
  static_assert(
      std::is_same<O, uintE>::value,
      "Currently apply_f must emit the same type as the count-type (uintE)");
  using W = typename Graph::weight_type;
  auto map_f = [](const uintE& i, const uintE& j, const W& wgh) {
    return gbbs::empty();
  };
  size_t m = vs.size();
  if (m == 0) {
    return vertexSubsetData<O>(vs.numNonzeros());
  }
  uintE empty_key = std::get<0>(ht.empty);
  auto oneHop = edgeMapInduced<gbbs::empty, Graph, VS>(GA, vs, map_f, cond_f,
                                                       empty_key, fl);
  oneHop.toSparse();

  auto key_f = [&](size_t i) -> uintE { return oneHop.vtx(i); };
  auto get_key = parlay::delayed_seq<uintE>(oneHop.size(), key_f);
  auto res =
      histogram<std::tuple<uintE, O> >(get_key, oneHop.size(), apply_f, ht);
  return vertexSubsetData<O>(vs.n, std::move(res));
}

// dense [read all neighbors]
template <class O, class Cond, class Apply, class VS, class Graph>
inline vertexSubsetData<O> edgeMapCount_dense(Graph& GA, VS& vs, Cond& cond_f,
                                              Apply& apply_f,
                                              const flags fl = 0) {
  using W = typename Graph::weight_type;
  size_t n = GA.n;
  size_t m = vs.size();
  if (m == 0) {
    return vertexSubsetData<O>(vs.numNonzeros());
  }
  using OT = std::tuple<bool, O>;
  std::cout << "Making vs dense!" << std::endl;
  vs.toDense();
  std::cout << "Made it dense!" << std::endl;

  debug(std::cout << "running dense" << std::endl << std::endl;);

  auto count_f = [&](uintE u, uintE v, W& wgh) -> size_t {
    return static_cast<size_t>(vs.isIn(v));
  };

  if (fl & no_output) {
    parallel_for(0, n,
                 [&](size_t i) {
                   if (cond_f(i)) {
                     auto neighbors = (fl & in_edges)
                                          ? GA.get_vertex(i).in_neighbors()
                                          : GA.get_vertex(i).out_neighbors();
                     size_t count = neighbors.count(count_f);

                     auto tup = std::make_tuple(i, count);
                     if (count > 0) {
                       apply_f(tup);
                     }
                   }
                 },
                 1);
    return vertexSubsetData<O>(n);
  } else {
    auto out = sequence<OT>::uninitialized(n);
    parallel_for(0, n,
                 [&](size_t i) {
                   if (cond_f(i)) {
                     auto neighbors = (fl & in_edges)
                                          ? GA.get_vertex(i).in_neighbors()
                                          : GA.get_vertex(i).out_neighbors();
                     size_t count = neighbors.count(count_f);
                     auto tup = std::make_tuple(i, count);
                     if (count > 0) {
                       auto applied_val = apply_f(tup);
                       std::get<0>(out[i]) = applied_val.has_value();
                       std::get<1>(out[i]) = std::get<1>(*applied_val);
                     } else {
                       std::get<0>(out[i]) = false;
                     }
                   }
                 },
                 1);
    return vertexSubsetData<O>(n, std::move(out));
  }
}

template <class O, class Cond, class Apply, class VS, class Graph>
inline vertexSubsetData<O> edgeMapCount(Graph& GA, VS& vs, Cond& cond_f,
                                        Apply& apply_f,
                                        hist_table<uintE, O>& ht,
                                        const flags fl = 0,
                                        long threshold = -1) {
  if (fl & no_dense) {
    return edgeMapCount_sparse<O>(GA, vs, ht, cond_f, apply_f, fl);
  } else if (fl & dense_only) {
    return edgeMapCount_dense<O>(GA, vs, cond_f, apply_f, fl);
  }
  vs.toSparse();
  auto degree_f = [&](size_t i) -> size_t {
    auto neighbors = (fl & in_edges) ? GA.get_vertex(i).in_neighbors()
                                     : GA.get_vertex(i).out_neighbors();
    return neighbors.get_virtual_degree();
  };
  auto degree_imap = parlay::delayed_seq<size_t>(vs.size(), degree_f);
  auto out_degrees = parlay::reduce(degree_imap);
  size_t degree_threshold = threshold;
  if (threshold == -1) degree_threshold = GA.m / 20;
  if (vs.size() + out_degrees > degree_threshold) {
    // dense
    return edgeMapCount_dense<O>(GA, vs, cond_f, apply_f, fl);
  } else {
    // sparse
    return edgeMapCount_sparse<O>(GA, vs, ht, cond_f, apply_f, fl);
  }
}

template <class O, class Apply, class VS, class Graph>
inline vertexSubsetData<O> edgeMapCount(Graph& GA, VS& vs, Apply& apply_f,
                                        hist_table<uintE, O>& ht,
                                        const flags fl = 0,
                                        long threshold = -1) {
  auto cond_true = [&](const uintE& u) { return true; };
  return edgeMapCount(GA, vs, cond_true, apply_f, ht, fl, threshold);
}

template <class O, class Cond, class Apply, class VS, class Graph>
inline vertexSubsetData<O> srcCount(Graph& GA, VS& vs, Cond cond_f,
                                    Apply apply_f, const flags fl = 0) {
  size_t n = GA.n;
  if (vs.dense()) {
    using OT = std::tuple<bool, O>;
    auto out = sequence<OT>(n);
    parallel_for(0, n, [&](size_t i) {
      if (vs.isIn(i)) {
        if (cond_f(i)) {
          auto tup = {true, GA.get_vertex(i).out_degree()};
          out[i] = apply_f(tup);
        }
      } else {
        std::get<0>(out[i]) = false;
      }
    });
    return vertexSubsetData<O>(n, std::move(out));
  } else {
    using OT = std::tuple<uintE, O>;
    auto out = sequence<OT>(vs.size());
    parallel_for(0, vs.size(), [&](size_t i) {
      uintE v = vs.vtx(i);
      out[i] = apply_f({v, GA.get_vertex(v).out_degree()});
    });
    return vertexSubsetData<O>(n, std::move(out));
  }
}

template <class V, class Graph>
struct EdgeMap {
  using K = uintE;  // keys are always uintE's (vertex-identifiers)
  using KV = std::tuple<K, V>;
  using W = typename Graph::weight_type;
  Graph& G;
  hist_table<K, V> ht;

  EdgeMap(Graph& _G, KV _empty,
          size_t ht_size = std::numeric_limits<size_t>::max())
      : G(_G) {
    if (ht_size == std::numeric_limits<size_t>::max()) {
      ht_size = G.m / 20;
    }
    ht = hist_table<K, V>(_empty, ht_size);
  }

  // sparse [write out neighbors]
  // id: M
  // map_f: (uintE v, uintE ngh) -> M
  // reduce_f: (M, tuple(uintE ngh, M ngh_val)) -> M
  // apply_f: (uintE ngh, M reduced_val) -> std::optional<O>
  template <class O, class M, class Map, class Reduce, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapReduce_sparse(VS& vs, Map& map_f,
                                                  Reduce& reduce_f,
                                                  Apply& apply_f, M id,
                                                  const flags fl) {
    size_t m = vs.size();
    size_t n = G.n;
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }

    auto cond_f = [&](const uintE& u) { return true; };

    uintE empty_key = std::get<0>(ht.empty);
    auto oneHop =
        edgeMapInduced<M, Graph, VS>(G, vs, map_f, cond_f, empty_key, fl);
    oneHop.toSparse();

    auto elm_f = [&](size_t i) { return oneHop.vtxAndData(i); };
    auto get_elm = parlay::delayed_seq<std::tuple<K, M> >(oneHop.size(), elm_f);
    auto key_f = [&](size_t i) -> uintE { return oneHop.vtx(i); };
    auto get_key = parlay::delayed_seq<uintE>(oneHop.size(), key_f);

    auto q = [&](sequentialHT<K, V>& S, std::tuple<K, M> v) -> void {
      S.template insertF<M>(v, reduce_f);
    };
    auto res = histogram_reduce<std::tuple<K, M>, std::tuple<K, O> >(
        get_elm, get_key, oneHop.size(), q, apply_f, ht);
    auto ret = vertexSubsetData<O>(vs.n, std::move(res));
    if (fl & no_output) {
      return vertexSubsetData<O>(n);
    }
    return ret;
  }

  // dense [read all neighbors]
  // id: M
  // map_f: (uintE v, uintE ngh) -> M
  // reduce_f: M * M -> M
  // apply_f: (uintE ngh, M reduced_val) -> std::optional<O>
  template <class O, class M, class Cond, class Map, class Reduce, class Apply,
            class VS>
  inline vertexSubsetData<O> edgeMapReduce_dense(VS& vs, Cond& cond_f,
                                                 Map& map_f, Reduce& reduce_f,
                                                 Apply& apply_f, M id,
                                                 const flags fl) {
    size_t n = G.n;
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    using OT = typename vertexSubsetData<O>::D;

    auto red_monoid = parlay::make_monoid(reduce_f, id);
    if (fl & no_output) {
      parallel_for(0, n,
                   [&](size_t i) {
                     if (cond_f(i)) {
                       auto neighbors = (fl & in_edges)
                                            ? G.get_vertex(i).in_neighbors()
                                            : G.get_vertex(i).out_neighbors();
                       M reduced_val = neighbors.reduce(map_f, red_monoid);
                       auto tup = std::make_tuple(i, reduced_val);
                       apply_f(tup);
                     }
                   },
                   1);
      return vertexSubsetData<O>(n);
    } else {
      auto out = sequence<OT>(n);
      parallel_for(0, n,
                   [&](size_t i) {
                     if
                       constexpr(!std::is_same<O, gbbs::empty>()) {
                         std::get<0>(out[i]) = false;
                       }
                     else {
                       out[i] = false;
                     }
                     if (cond_f(i)) {
                       auto neighbors = (fl & in_edges)
                                            ? G.get_vertex(i).in_neighbors()
                                            : G.get_vertex(i).out_neighbors();
                       M reduced_val = neighbors.reduce(map_f, red_monoid);
                       auto tup = std::make_tuple(i, reduced_val);
                       auto applied_val = apply_f(tup);
                       if (applied_val.has_value()) {
                         if
                           constexpr(!std::is_same<O, gbbs::empty>()) {
                             std::get<0>(out[i]) = true;
                             std::get<1>(out[i]) = std::get<1>(*applied_val);
                           }
                         else {
                           out[i] = true;
                         }
                       }
                     }
                   },
                   1);
      return vertexSubsetData<O>(n, std::move(out));
    }
  }

  // id: M
  // cond_f: uintE -> bool
  // map_f: (uintE v, uintE ngh) -> M
  // reduce_f: M * M -> M
  // apply_f: (uintE ngh, M reduced_val) -> std::optional<O>
  template <class O, class M, class Cond, class Map, class Reduce, class Apply,
            class VS>
  inline vertexSubsetData<O> edgeMapReduce(VS& vs, Cond& cond_f, Map& map_f,
                                           Reduce& reduce_f, Apply& apply_f,
                                           M id, intT threshold = -1,
                                           const flags fl = 0) {
    assert(false);
    // Currently unimplemented: only using edgeMapReduce_dense. TODO: finish
    // implementation of this method (if needed).
    exit(0);
    vs.toSparse();
    auto degree_f = [&](size_t i) {
      auto neighbors = (fl & in_edges)
                           ? G.get_vertex(vs.vtx(i)).in_neighbors()
                           : G.get_vertex(vs.vtx(i)).out_neighbors();
      return neighbors.get_virtual_degree();
    };
    auto degree_imap = parlay::delayed_seq<uintE>(vs.size(), degree_f);
    auto out_degrees = parlay::reduce(degree_imap);
    if (threshold == -1) threshold = G.m / 20;
    if (vs.size() + out_degrees > threshold) {
      // dense
      return edgeMapReduce_dense<O, M>(vs, cond_f, map_f, reduce_f, apply_f, id,
                                       fl);
    } else {
      // sparse
      auto red_ht = [&](M cur, std::tuple<uintE, M> nxt) -> M {
        return reduce_f(cur, std::get<1>(nxt));
      };
      return edgeMapReduce_sparse<O, M>(vs, map_f, red_ht, apply_f, id, fl);
    }
  }

  // ============================= Edge Map Count
  // ===============================

  // sparse [write out neighbors]
  template <class O, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapCount_sparse(VS& vs, Apply& apply_f,
                                                 const flags fl = 0) {
    auto map_f = [](const uintE& i, const uintE& j, const W& wgh) {
      return gbbs::empty();
    };
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    auto cond_f = [&](const uintE& u) { return true; };
    uintE empty_key = std::get<0>(ht.empty);
    auto oneHop = edgeMapInduced<gbbs::empty, Graph, VS>(G, vs, map_f, cond_f,
                                                         empty_key, fl);
    oneHop.toSparse();

    auto key_f = [&](size_t i) -> uintE { return oneHop.vtx(i); };
    auto get_key = parlay::delayed_seq<uintE>(oneHop.size(), key_f);
    auto res =
        histogram<std::tuple<uintE, O> >(get_key, oneHop.size(), apply_f, ht);
    return vertexSubsetData<O>(vs.n, std::move(res));
  }

  // dense [read all neighbors]
  template <class O, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapCount_dense(VS& vs, Apply& apply_f,
                                                const flags fl = 0) {
    size_t n = G.n;
    size_t m = vs.size();
    if (m == 0) {
      return vertexSubsetData<O>(vs.numNonzeros());
    }
    using OT = std::tuple<bool, O>;
    vs.toDense();

    debug(std::cout << "running dense" << std::endl << std::endl;);

    auto count_f = [&](uintE u, uintE v, W& wgh) -> size_t {
      return static_cast<size_t>(vs.isIn(v));
    };

    if (fl & no_output) {
      parallel_for(0, n,
                   [&](size_t i) {
                     auto neighbors = (fl & in_edges)
                                          ? G.get_vertex(i).in_neighbors()
                                          : G.get_vertex(i).out_neighbors();
                     size_t count = neighbors.count(count_f);
                     auto tup = std::make_tuple(i, count);
                     if (count > 0) {
                       apply_f(tup);
                     }
                   },
                   1);
      return vertexSubsetData<O>(n);
    } else {
      auto out = sequence<OT>(n);
      parallel_for(0, n,
                   [&](size_t i) {
                     auto neighbors = (fl & in_edges)
                                          ? G.get_vertex(i).in_neighbors()
                                          : G.get_vertex(i).out_neighbors();
                     size_t count = neighbors.count(count_f);
                     auto tup = std::make_tuple(i, count);
                     if (count > 0) {
                       auto applied_val = apply_f(tup);
                       if (applied_val.has_value()) {
                         std::get<0>(out[i]) = true;
                         std::get<1>(out[i]) = std::get<1>(*applied_val);
                       } else {
                         std::get<0>(out[i]) = false;
                       }
                     } else {
                       std::get<0>(out[i]) = false;
                     }
                   },
                   1);
      return vertexSubsetData<O>(n, std::move(out));
    }
  }

  template <class O, class Apply, class VS>
  inline vertexSubsetData<O> edgeMapCount(VS& vs, Apply& apply_f,
                                          const flags fl = 0,
                                          long threshold = -1) {
    vs.toSparse();
    auto degree_f = [&](size_t i) -> size_t {
      auto neighbors = (fl & in_edges)
                           ? G.get_vertex(vs.vtx(i)).in_neighbors()
                           : G.get_vertex(vs.vtx(i)).out_neighbors();
      return neighbors.get_virtual_degree();
    };
    auto degree_imap = parlay::delayed_seq<size_t>(vs.size(), degree_f);
    auto out_degrees = parlay::reduce(degree_imap);
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
};

}  // namespace gbbs
