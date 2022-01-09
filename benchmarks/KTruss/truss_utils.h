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

namespace gbbs {
namespace truss_utils {

struct edge_id {
  uintE a;
  uintE b;
  edge_id()
      : a(std::numeric_limits<uintE>::max()),
        b(std::numeric_limits<uintE>::max()) {}
  edge_id(uintE _a, uintE _b) : a(std::min(_a, _b)), b(std::max(_a, _b)) {}
  bool operator==(const edge_id& other) const {
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

  inline size_t hashToRangeValue(size_t h) { return h & mask; }
  inline size_t firstIndex(V& v) {
    return hashToRangeValue(parlay::hash32_3(v));
  }
  inline size_t incrementIndex(uintT h) { return hashToRangeValue(h + 1); }

  size_t size() { return m; }

  static void clearA(V* A, long n, V kv) {
    parallel_for(0, n, kDefaultGranularity, [&](size_t i) { A[i] = kv; });
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  valueHT(size_t _m, V _empty)
      : m((size_t)1 << parlay::log2_up((size_t)(_m))),
        mask(m - 1),
        empty(_empty) {
    table = gbbs::new_array_no_init<V>(m);
    clearA(table, m, empty);
    alloc = true;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  valueHT(size_t _m, V _empty, V* _table)
      : m(_m), mask(m - 1), empty(_empty), table(_table), alloc(false) {}

  void del() {
    if (alloc) {
      gbbs::free_array(table);
    }
    alloc = false;
  }

  bool insert(V v) {
    size_t h = firstIndex(v);
    size_t iters = 0;
    while (1) {
      if (table[h] == empty &&
          gbbs::atomic_compare_and_swap(&table[h], empty, v)) {
        return 1;
      }
      h = incrementIndex(h);
      iters++;
    }
    return 0;  // should never get here
  }

  sequence<V> entries() {
    auto pred = [&](V t) { return t != empty; };
    auto table_seq = parlay::delayed_seq<V>(table, m);
    return parlay::filter(table_seq, pred);
  }

  void clear() {
    parallel_for(0, m, [&](size_t i) { table[i] = empty; });
  }
};

template <class K, class V, class GS>
class multi_table {
 public:
  using T = std::tuple<K, V>;

  size_t n;  // number of internal tables
  bool alloc;

  T* big_table;     // start of big table
  size_t big_size;  // size of the big table
  V empty_val;      // value used for empty cells

  sequence<size_t> offsets;

  class inner_table {
   public:
    size_t mask;  // pow2 - 1
    T* table;
    K empty_key;

    inline size_t hashToRange(size_t h) { return h & mask; }
    inline size_t firstIndex(K& k) {
      return hashToRange(parlay::hash32(k));
    }  // hacks for now
    inline size_t incrementIndex(size_t h) { return hashToRange(h + 1); }

    // Precondition: k is not present in the table
    inline void insert(std::tuple<K, V>& kv) {
      K& k = std::get<0>(kv);
      size_t h = firstIndex(k);
      while (true) {
        if (std::get<0>(table[h]) == empty_key) {
          if (gbbs::atomic_compare_and_swap(&std::get<0>(table[h]), empty_key,
                                            k)) {
            std::get<1>(table[h]) = std::get<1>(kv);  // insert value
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
          gbbs::write_add(&std::get<1>(table[h]), static_cast<V>(1));
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

  inner_table* tables;  // one object per n (16 bytes) (note we can probably red
                        // to 8 bytes each)

  multi_table() : n(0) { alloc = false; }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  multi_table(size_t n, V empty_val, GS size_func)
      : n(n), empty_val(empty_val) {
    // compute offsets
    offsets = sequence<size_t>(n + 1);
    parallel_for(0, n, [&](size_t i) {
      size_t table_elms = size_func(i);
      offsets[i] = (1 << parlay::log2_up((size_t)(table_elms * 1.2))) +
                   2;  // 2 cell padding (l, r)
    });
    offsets[n] = 0;
    size_t total_space = parlay::scan_inplace(offsets);
    std::cout << "total space = " << total_space << std::endl;
    std::cout << "empty val is " << empty_val << std::endl;

    big_table = gbbs::new_array_no_init<T>(total_space);
    big_size = total_space;

    tables = gbbs::new_array_no_init<inner_table>(n);
    parallel_for(0, n, [&](size_t i) {
      size_t off = offsets[i];
      size_t sz = offsets[i + 1] - off;
      sz -= 2;  // 2 cell's padding
      auto table_loc = big_table + offsets[i];

      tables[i].table = table_loc + 1;  // 1 cell padding (l
      tables[i].mask = sz - 1;
      tables[i].empty_key = i;

      // clear i's table
      V val = empty_val;
      auto empty_i = std::make_tuple((K)i, val);
      parallel_for(0, sz + 2, [&](size_t j) {  // sz + 2 cell padding
        table_loc[j] = empty_i;
      });
    });

    //    parallel_for(0, big_size, [&] (size_t i) {
    //      assert(std::get<1>(big_table[i]) == empty_val);
    //    });
  }

  inline size_t size() { return big_size; }

  inline void insert(K a, std::tuple<K, V> kv) { tables[a].insert(kv); }

  inline size_t idx(K a, K b) {
    K u = std::min(a, b);
    return offsets[u] + 1 + tables[u].idx(std::max(a, b));
  }

  inline void increment(K a, K b) {
    tables[std::min(a, b)].increment(std::max(a, b));
  }

  inline K u_for_id(size_t id) {
    // walk right until we hit an empty cell with val == empty_val
    // assumption is that id is not empty
    size_t idx = id + 1;
    //    assert(std::get<1>(big_table[id]) != empty_val);
    while (true) {
      if (std::get<1>(big_table[idx]) == empty_val) {
        return std::get<0>(big_table[idx]);  // return key
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

template <class F, class Graph>
struct countF {
  Graph& G;
  F f;
  using W = typename Graph::weight_type;
  countF(Graph& G, F f) : G(G), f(f) {}

  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    auto v_s = G.get_vertex(s).out_neighbors();
    auto v_d = G.get_vertex(d).out_neighbors();
    v_s.intersect_f(&v_d, f);
    return 1;
  }

  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    auto v_s = G.get_vertex(s).out_neighbors();
    auto v_d = G.get_vertex(d).out_neighbors();
    v_s.intersect_f(&v_d, f);
    return 1;
  }

  inline bool cond(uintE d) { return cond_true(d); }
};

template <class Graph, class VS, class F>
vertexSubset emdf(Graph& GA, VS& vs, F f, const flags& fl = 0) {
  return edgeMapDenseForward<gbbs::empty>(GA, vs, f, fl);
}

template <class F, class Graph>
void TCDirected(Graph& DG, F f) {
  size_t n = DG.n;
  auto frontier = sequence<bool>::uninitialized(n);
  parallel_for(0, n, [&](size_t i) { frontier[i] = 1; });
  vertexSubset Frontier(n, n, std::move(frontier));
  emdf(DG, Frontier, countF<F, Graph>(DG, f), no_output);
}

template <class Graph>
sequence<uintE> rankNodes(Graph& G) {
  size_t n = G.n;
  auto r = sequence<uintE>(n);
  auto o = sequence<uintE>(n);

  parallel_for(0, n, [&](size_t i) { o[i] = i; });
  parlay::sample_sort_inplace(make_slice(o), [&](const uintE u, const uintE v) {
    return G.get_vertex(u).out_degree() < G.get_vertex(v).out_degree();
  });
  parallel_for(0, n, [&](size_t i) { r[o[i]] = i; });
  return r;
}

//  template <class TT>
//  void increment_trussness(TT& ht, const uintE& u, const uintE& v) {
//    size_t loc_uv = ht.idx(std::make_tuple(std::min(u, v), std::max(u, v)));
//    auto val_loc = &std::get<1>(ht.table[loc_uv]);
//    gbbs::write_add(val_loc, (uintE)1);
//  }

//  template <class trussness_t, class edge_t>
//  inline bool should_remove(uintE k, trussness_t trussness_uv, trussness_t
//  trussness_uw, trussness_t trussness_vw,
//                            edge_t uv_id, edge_t uw_id, edge_t vw_id) {
//    constexpr trussness_t TOP_BIT =
//    static_cast<size_t>(std::numeric_limits<long>::max()) + 1;
//
//    // this triangle was removed in a previous round.
//    if (k > std::min(trussness_uw, trussness_vw)) { return false; }
//
//    if (trussness_uw & TOP_BIT) { // uw also removed on this round
//      return (uv_id < uw_id) && !(trussness_vw & TOP_BIT); // sym break, make
//      sure vw also not removed.
//    } else if (trussness_vw & TOP_BIT) {
//      return (uv_id < vw_id); // sym break, we know uw does not have top bit
//      set.
//    }
//
//    // neither other edge is removed this round. can safely remove this
//    // triangle as uv.
//    return true;
//  }

template <class trussness_t, class edge_t>
inline bool should_remove(uintE k, trussness_t trussness_uv,
                          trussness_t trussness_uw, trussness_t trussness_vw,
                          edge_t uv_id, edge_t uw_id, edge_t vw_id) {
  // this triangle was removed in a previous round.
  if (k > std::min(trussness_uw, trussness_vw)) {
    return false;
  }

  if (trussness_uw == k) {  // uw also removed on this round
    return (uv_id < uw_id) &&
           !(trussness_vw == k);  // sym break, make sure vw also not removed.
  } else if (trussness_vw == k) {
    return (uv_id < vw_id);  // sym break, we know uw does not have top bit set.
  }

  // neither other edge is removed this round. can safely remove this
  // triangle as uv.
  return true;
}

// get_trussness_and_id: (uintE, uintE) -> (trussness, id)
// Intersect u and v's neighbors. Check if we should remove an edge, and if
// so, insert it into decrement_tab (concurrent write)
template <class edge_t, class trussness_t, class HT, class Trussness,
          class Graph>
void decrement_trussness(Graph& G, edge_t id, uintE u, uintE v,
                         HT& decrement_tab, Trussness& get_trussness_and_id,
                         uintE k) {
  trussness_t trussness_uv = k;
  edge_t uv_id = id;

  auto add_f = [&](uintE* ct, const std::tuple<uintE, uintE>& tup) {
    gbbs::fetch_and_add(ct, (uintE)1);
  };

  size_t ctr = 0;
  auto f = [&](uintE __u, uintE __v, uintE w) {  // w in both N(u), N(v)
    trussness_t trussness_uw, trussness_vw;
    edge_t uw_id, vw_id;
    std::tie(trussness_uw, uw_id) = get_trussness_and_id(u, w);
    std::tie(trussness_vw, vw_id) = get_trussness_and_id(v, w);

    if (should_remove(k, trussness_uv, trussness_uw, trussness_vw, uv_id, uw_id,
                      vw_id)) {
      ctr++;
      if (trussness_uw > k) {
        decrement_tab.insert_f(std::make_tuple(uw_id, (uintE)1), add_f);
        //          decrement_tab.insert(uw_id);
      }
      if (trussness_vw > k) {
        decrement_tab.insert_f(std::make_tuple(vw_id, (uintE)1), add_f);
        //          decrement_tab.insert(vw_id);
      }
    }
  };
  auto v_v = G.get_vertex(v).out_neighbors();
  G.get_vertex(u).out_neighbors().intersect_f_par(&v_v, f);
}

}  // namespace truss_utils
}  // namespace gbbs
