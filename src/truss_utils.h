// Copyright (c) 2019 Laxman Dhulipala, Guy Blelloch, and Julian Shun
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

#include "lib/sample_sort.h"
#include "lib/utilities.h"

namespace truss_utils {

  struct edge_id {
    uintE a;
    uintE b;
    edge_id() : a(std::numeric_limits<uintE>::max()), b(std::numeric_limits<uintE>::max()) {}
    edge_id(uintE _a, uintE _b) : a(std::min(_a, _b)), b(std::max(_a, _b)) { }
    bool operator==(const edge_id &other) const {
      return (a == other.a) && (b == other.b);
    }
  };

  // A table that stores a collection of elements, including duplicates.
  template <class Val>
  class valueHT {
  public:
    typedef Val V;
    size_t m;
    size_t mask;
    V empty;
    V* table;
    bool alloc;

    uint32_t thomas_hash(uint32_t a) {
        a = (a ^ 61) ^ (a >> 16);
        a = a + (a << 3);
        a = a ^ (a >> 4);
        a = a * 0x27d4eb2d;
        a = a ^ (a >> 15);
        return a;
    }

    inline size_t hashToRangeValue(size_t h) {return h & mask;}
    inline size_t firstIndex(V& v) {return hashToRangeValue(pbbs::hash32(v));}
    inline size_t incrementIndex(uintT h) {return hashToRangeValue(h+1);}

    size_t size() {
      return m;
    }

    static void clearA(V* A, long n, V kv) {
      par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i)
                      { A[i] = kv; });
    }

    // Size is the maximum number of values the hash table will hold.
    // Overfilling the table could put it into an infinite loop.
    valueHT(size_t _m, V _empty) :
      m((size_t) 1 << pbbs::log2_up((size_t)(_m))),
      mask(m-1),
      empty(_empty) {
        table = pbbs::new_array_no_init<V>(m);
        clearA(table, m, empty);
        alloc=true;
    }

    // Size is the maximum number of values the hash table will hold.
    // Overfilling the table could put it into an infinite loop.
    valueHT(size_t _m, V _empty, V* _table) :
      m(_m), mask(m-1), empty(_empty), table(_table), alloc(false) { }

    void del() {
      if (alloc) {
        pbbs::free_array(table);
      }
      alloc=false;
    }

    bool insert(V v) {
      size_t h = firstIndex(v);
      while (1) {
        if(table[h] == empty && pbbs::CAS(&table[h],empty,v)) {
          return 1;
        }
        h = incrementIndex(h);
      }
      return 0; // should never get here
    }

    sequence<V> entries() {
      auto pred = [&](V t) { return t != empty; };
      auto table_seq = make_sequence<V>(table, m);
      return pbbs::filter(table_seq, pred);
    }

    void clear() {
      par_for(0, m, [&] (size_t i) { table[i] = empty; });
    }

  };

  template <class F, template <class W> class vertex, class W>
  struct countF {
    using w_vertex = vertex<W>;
    w_vertex* V;
    F f;
    countF(w_vertex* _V, F _f) : V(_V), f(_f) {}

    inline bool update(uintE s, uintE d) {
      V[s].intersect(&V[d], s, d, f);
      return 1;
    }

    inline bool updateAtomic(uintE s, uintE d) {
      V[s].intersect_f(&V[d], s, d, f);
      return 1;
    }

    inline bool cond(uintE d) { return cond_true(d); }
  };

  template <class vertex, class VS, class F>
  vertexSubset emdf(graph<vertex> GA, VS& vs, F f, const flags& fl = 0) {
    return edgeMapDenseForward<pbbs::empty>(GA, vs, f, fl);
  }

  template <class F, template <class W> class vertex, class W>
  void TCDirected(graph<vertex<W>>& DG, F f) {
    size_t n = DG.n;
    auto frontier = sequence<bool>(n);
    par_for(0, n, [&] (size_t i) { frontier[i] = 1; });
    vertexSubset Frontier(n, n, frontier.start());
    emdf(DG, Frontier, wrap_em_f<W>(countF<F, vertex, W>(DG.V, f)), no_output);
  }

  template <class vertex>
  sequence<uintE> rankNodes(vertex* V, size_t n) {
    auto r = sequence<uintE>(n);
    auto o = sequence<uintE>(n);

    par_for(0, n, [&] (size_t i) { o[i] = i; });
    pbbs::sample_sort(o.start(), n, [&](const uintE u, const uintE v) {
      return V[u].getOutDegree() < V[v].getOutDegree();
    });
    par_for(0, n, [&] (size_t i) { r[o[i]] = i; });
    return r;
  }

  template <class TT>
  void increment_trussness(TT& ht, const uintE& u, const uintE& v) {
    size_t loc_uv = ht.idx(std::make_tuple(std::min(u, v), std::max(u, v)));
    auto val_loc = &std::get<1>(ht.table[loc_uv]);
    pbbs::write_add(val_loc, (uintE)1);
  }

//  template <class trussness_t, class edge_t>
//  inline bool should_remove(uintE k, trussness_t trussness_uv, trussness_t trussness_uw, trussness_t trussness_vw,
//                            edge_t uv_id, edge_t uw_id, edge_t vw_id) {
//    constexpr trussness_t TOP_BIT = static_cast<size_t>(std::numeric_limits<long>::max()) + 1;
//
//    // this triangle was removed in a previous round.
//    if (k > std::min(trussness_uw, trussness_vw)) { return false; }
//
//    if (trussness_uw & TOP_BIT) { // uw also removed on this round
//      return (uv_id < uw_id) && !(trussness_vw & TOP_BIT); // sym break, make sure vw also not removed.
//    } else if (trussness_vw & TOP_BIT) {
//      return (uv_id < vw_id); // sym break, we know uw does not have top bit set.
//    }
//
//    // neither other edge is removed this round. can safely remove this
//    // triangle as uv.
//    return true;
//  }

  template <class trussness_t, class edge_t>
  inline bool should_remove(uintE k, trussness_t trussness_uv, trussness_t trussness_uw, trussness_t trussness_vw,
                            edge_t uv_id, edge_t uw_id, edge_t vw_id) {

    // this triangle was removed in a previous round.
    if (k > std::min(trussness_uw, trussness_vw)) { return false; }

    if (trussness_uw == k) { // uw also removed on this round
      return (uv_id < uw_id) && !(trussness_vw == k); // sym break, make sure vw also not removed.
    } else if (trussness_vw == k) {
      return (uv_id < vw_id); // sym break, we know uw does not have top bit set.
    }

    // neither other edge is removed this round. can safely remove this
    // triangle as uv.
    return true;
  }


  // get_trussness_and_id: (uintE, uintE) -> (trussness, id)
  // Intersect u and v's neighbors. Check if we should remove an edge, and if
  // so, insert it into decrement_tab (concurrent write)
  template <class edge_t, class trussness_t, class Trussness, class Graph>
  void decrement_trussness(Graph& G, edge_t id, uintE u, uintE v, valueHT<edge_t>& decrement_tab, Trussness& get_trussness_and_id, uintE k) {

    trussness_t trussness_uv = k; edge_t uv_id = id;

    auto f = [&] (uintE __u, uintE __v, uintE w) { // w in both N(u), N(v)
      trussness_t trussness_uw, trussness_vw;
      edge_t uw_id, vw_id;
      std::tie(trussness_uw, uw_id) = get_trussness_and_id(u, w);
      std::tie(trussness_vw, vw_id) = get_trussness_and_id(v, w);

      if (should_remove(k, trussness_uv, trussness_uw, trussness_vw, uv_id, uw_id, vw_id)) {
        decrement_tab.insert(uw_id);
        decrement_tab.insert(vw_id);
      }
    };
    G.V[u].intersect_f_par(&(G.V[v]), u, v, f);
  }

}; // namespace truss_utils

