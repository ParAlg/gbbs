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

#include "pbbslib/sample_sort.h"
#include "pbbslib/utilities.h"

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
    inline size_t firstIndex(V& v) {return hashToRangeValue(pbbslib::hash32_3(v));}
    inline size_t incrementIndex(uintT h) {return hashToRangeValue(h+1);}

    size_t size() {
      return m;
    }

    static void clearA(V* A, long n, V kv) {
      par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                      { A[i] = kv; });
    }

    // Size is the maximum number of values the hash table will hold.
    // Overfilling the table could put it into an infinite loop.
    valueHT(size_t _m, V _empty) :
      m((size_t) 1 << pbbslib::log2_up((size_t)(_m))),
      mask(m-1),
      empty(_empty) {
        table = pbbslib::new_array_no_init<V>(m);
        clearA(table, m, empty);
        alloc=true;
    }

    // Size is the maximum number of values the hash table will hold.
    // Overfilling the table could put it into an infinite loop.
    valueHT(size_t _m, V _empty, V* _table) :
      m(_m), mask(m-1), empty(_empty), table(_table), alloc(false) { }

    void del() {
      if (alloc) {
        pbbslib::free_array(table);
      }
      alloc=false;
    }

    bool insert(V v) {
      size_t h = firstIndex(v);
      size_t iters = 0;
      while (1) {
        if(table[h] == empty && pbbslib::CAS(&table[h],empty,v)) {
          return 1;
        }
        h = incrementIndex(h);
        iters++;
      }
      return 0; // should never get here
    }

    sequence<V> entries() {
      auto pred = [&](V t) { return t != empty; };
      auto table_seq = pbbslib::make_sequence<V>(table, m);
      return pbbslib::filter(table_seq, pred);
    }

    void clear() {
      par_for(0, m, [&] (size_t i) { table[i] = empty; });
    }

  };


  template <class K, class V, class GS>
  class multi_table {
   public:
    using T = std::tuple<K, V>;

    size_t n; // number of internal tables
    bool alloc;

    T* big_table; // start of big table
    size_t big_size; // size of the big table
    V empty_val; // value used for empty cells

    sequence<size_t> offsets;

    class inner_table {
      public:
        size_t mask; // pow2 - 1
        T* table;
        K empty_key;

        inline size_t hashToRange(size_t h) { return h & mask; }
        inline size_t firstIndex(K& k) { return hashToRange(pbbslib::hash32(k)); } // hacks for now
        inline size_t incrementIndex(size_t h) { return hashToRange(h + 1); }

        // Precondition: k is not present in the table
        inline void insert(std::tuple<K, V>& kv) {
          K& k = std::get<0>(kv);
          size_t h = firstIndex(k);
          while (true) {
            if (std::get<0>(table[h]) == empty_key) {
              if (pbbslib::CAS(&std::get<0>(table[h]), empty_key, k)) {
                std::get<1>(table[h]) = std::get<1>(kv); // insert value
                return;
              }
            }
            h = incrementIndex(h);
          }
        }

        // Precondition: k is present in the table
        inline void increment(K k) {
          size_t h = firstIndex(k);
          while (true) {
            if (std::get<0>(table[h]) == k) {
              pbbslib::write_add(&std::get<1>(table[h]), static_cast<V>(1));
              return;
            }
            h = incrementIndex(h);
          }
        }

        // Precondition: k must be present in the table
        inline size_t idx(K k) {
          size_t h = firstIndex(k);
          while (true) {
            if (std::get<0>(table[h]) == k) {
              return h;
            }
            h = incrementIndex(h);
          }
        }
    };

    inner_table* tables; // one object per n (16 bytes) (note we can probably red to 8 bytes each)

    multi_table() : n(0) {
      alloc = false;
    }

    // Size is the maximum number of values the hash table will hold.
    // Overfilling the table could put it into an infinite loop.
    multi_table(size_t n, V empty_val, GS size_func) : n(n), empty_val(empty_val) {
      // compute offsets
      offsets = sequence<size_t>(n+1);
      par_for(0, n, [&] (size_t i) {
        size_t table_elms = size_func(i);
        offsets[i] = (1 << pbbslib::log2_up((size_t)(table_elms*1.2))) + 2; // 2 cell padding (l, r)
      });
      offsets[n] = 0;
      size_t total_space = pbbslib::scan_add_inplace(offsets);
      std::cout << "total space = " << total_space << std::endl;
      std::cout << "empty val is " << empty_val << std::endl;

      big_table = pbbslib::new_array_no_init<T>(total_space);
      big_size = total_space;

      tables = pbbslib::new_array_no_init<inner_table>(n);
      par_for(0, n, [&] (size_t i) {
        size_t off = offsets[i];
        size_t sz = offsets[i+1] - off;
        sz -= 2; // 2 cell's padding
        auto table_loc = big_table + offsets[i];

        tables[i].table = table_loc + 1; // 1 cell padding (l
        tables[i].mask = sz-1;
        tables[i].empty_key = i;

        // clear i's table
        V val = empty_val;
        auto empty_i = std::make_tuple((K)i, val);
        par_for(0, sz+2, [&] (size_t j) { // sz + 2 cell padding
          table_loc[j] = empty_i;
        });
      });

  //    par_for(0, big_size, [&] (size_t i) {
  //      assert(std::get<1>(big_table[i]) == empty_val);
  //    });

    }

    inline void check_consistency() {

      par_for(0, n, [&] (size_t i) {
        size_t off = offsets[i];
        size_t sz = offsets[i+1] - off;
        sz -= 2; // 2 cell's padding
        auto table_loc = big_table + offsets[i];
        auto table_start = table_loc + 1;
        auto table_end = table_loc + sz + 1;
      });

      par_for(0, big_size, [&] (size_t i) {
        if (std::get<1>(big_table[i]) != empty_val) {
          size_t idx = u_for_id(i);
        }
      });
    }

    inline size_t size() {
      return big_size;
    }

    inline void insert(K a, std::tuple<K, V> kv) {
      tables[a].insert(kv);
    }

    inline size_t idx(K a, K b) {
      K u = std::min(a,b);
      return offsets[u] + 1 + tables[u].idx(std::max(a, b));
    }

    inline void increment(K a, K b) {
      tables[std::min(a,b)].increment(std::max(a, b));
    }

    inline K u_for_id (size_t id) {
      // walk right until we hit an empty cell with val == empty_val
      // assumption is that id is not empty
      size_t idx = id + 1;
  //    assert(std::get<1>(big_table[id]) != empty_val);
      while (true) {
        if (std::get<1>(big_table[idx]) == empty_val) {
          return std::get<0>(big_table[idx]); // return key
        }
        idx++;
  //      if (idx > big_size) {
  //        assert(false);
  //      }
      }
    }

  };

  template <class K, class V, class GS>
  inline multi_table<K, V, GS> make_multi_table(size_t n, V empty_val, GS gs) {
    return multi_table<K, V, GS>(n, empty_val, gs);
  }




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
    return edgeMapDenseForward<pbbslib::empty>(GA, vs, f, fl);
  }

  template <class F, template <class W> class vertex, class W>
  void TCDirected(graph<vertex<W>>& DG, F f) {
    size_t n = DG.n;
    auto frontier = sequence<bool>(n);
    par_for(0, n, [&] (size_t i) { frontier[i] = 1; });
    vertexSubset Frontier(n, n, frontier.begin());
    emdf(DG, Frontier, wrap_em_f<W>(countF<F, vertex, W>(DG.V, f)), no_output);
  }

  template <class vertex>
  sequence<uintE> rankNodes(vertex* V, size_t n) {
    auto r = sequence<uintE>(n);
    auto o = sequence<uintE>(n);

    par_for(0, n, [&] (size_t i) { o[i] = i; });
    pbbslib::sample_sort(o.begin(), n, [&](const uintE u, const uintE v) {
      return V[u].getOutDegree() < V[v].getOutDegree();
    });
    par_for(0, n, [&] (size_t i) { r[o[i]] = i; });
    return r;
  }

//  template <class TT>
//  void increment_trussness(TT& ht, const uintE& u, const uintE& v) {
//    size_t loc_uv = ht.idx(std::make_tuple(std::min(u, v), std::max(u, v)));
//    auto val_loc = &std::get<1>(ht.table[loc_uv]);
//    pbbslib::write_add(val_loc, (uintE)1);
//  }


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

    size_t ctr = 0;
    auto f = [&] (uintE __u, uintE __v, uintE w) { // w in both N(u), N(v)
      trussness_t trussness_uw, trussness_vw;
      edge_t uw_id, vw_id;
      std::tie(trussness_uw, uw_id) = get_trussness_and_id(u, w);
      std::tie(trussness_vw, vw_id) = get_trussness_and_id(v, w);

      if (should_remove(k, trussness_uv, trussness_uw, trussness_vw, uv_id, uw_id, vw_id)) {
//        if (u == 14073 && v == 744491) {
//          cout << "trussness = " << trussness_uv << " " << trussness_uw << " " << trussness_vw << endl;
//          cout << "ids: " << uv_id << " " << uw_id << " " << vw_id << endl;
//        }
        ctr++;
        decrement_tab.insert(uw_id);
        decrement_tab.insert(vw_id);
//        if (uw_id == 2099851) {
//          std::cout << "decrementing our guy from " << uv_id << " " << uw_id << " " << vw_id << std::endl;
//        }
      }
    };
    G.V[u].intersect_f_par(&(G.V[v]), u, v, f);
  }

}; // namespace truss_utils

