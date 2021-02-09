// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010-2016 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include <tuple>
#include <cassert>

#include "gbbs/bridge.h"

namespace pbbslib {

template <class K, class V, class KeyHash>
class sparse_table {
 public:
  using T = std::tuple<K, V>;

  size_t m;
  size_t mask;
  T empty;
  K empty_key;
  T* table;
  bool alloc;
  KeyHash key_hash;

  size_t size() const {
    return m;
  }

  static void clearA(T* A, long n, T kv) {
    parallel_for(0, n, [&] (size_t i) { A[i] = kv; });
  }

  inline size_t hashToRange(size_t h) const { return h & mask; }
  inline size_t firstIndex(K& k) const { return hashToRange(key_hash(k)); }
  inline size_t incrementIndex(size_t h) const { return hashToRange(h + 1); }

  void del() {
    if (alloc) {
      gbbs::free_array(table, m);
      alloc = false;
    }
  }

  // Does not copy elements in the current table; just resize the underlying
  // table.
  // Incoming must be a power of two
  void resize_no_copy(size_t incoming) {
    if (incoming > m) {
      if (alloc) {
        gbbs::free_array(table, m);
      }
      std::cout << "# Resizing decrement table, was: " << m;
      m = incoming;
      mask = m - 1;
      table = gbbs::new_array_no_init<T>(m);
      clearA(table, m, empty);
      alloc = true;
      std::cout << "#  is now: " << m << std::endl;
    }
  }

  // Incoming must be a power of two
  void resize(size_t incoming) {
    if (incoming > m) {
      bool old_alloc = alloc;
      auto old_table = table;
      size_t old_m = m;
      std::cout << "# Resizing table, was: " << m;
      size_t new_size = 1 << parlay::log2_up(2*(m + incoming));
      m = new_size;
      mask = m - 1;
      table = gbbs::new_array_no_init<T>(m);
      clearA(table, m, empty);
      alloc = true;
      std::cout << "#  is now: " << m << std::endl;
      parallel_for(0, old_m, [&] (size_t i) {
        if (std::get<0>(old_table[i]) != empty_key) {
          insert(old_table[i]);
        }
      });

      if (old_alloc) {
        gbbs::free_array(old_table, old_m);
      }
    }
  }


  sparse_table() : m(0) {
    mask = 0;
    alloc = false;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_table(size_t _m, T _empty, KeyHash _key_hash, long inp_space_mult=-1)
      : empty(_empty),
        empty_key(std::get<0>(empty)),
        key_hash(_key_hash) {
    double space_mult = 1.1;
    if (inp_space_mult != -1) space_mult = inp_space_mult;
    m = (size_t)1 << parlay::log2_up((size_t)(space_mult * _m) + 1);
    mask = m - 1;
    table = gbbs::new_array_no_init<T>(m);
    clearA(table, m, empty);
    alloc = true;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_table(size_t _m, T _empty, KeyHash _key_hash, T* _tab, bool clear=true)
      : m(_m),
        mask(m - 1),
        empty(_empty),
        empty_key(std::get<0>(empty)),
        table(_tab),
        key_hash(_key_hash) {
    if (clear) {
      clearA(table, m, empty);
    }
    alloc = false;
  }

  // Pre-condition: k must be present in T.
  inline size_t idx(K k) const {
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == k) {
        return h;
      }
      h = incrementIndex(h);
    }
  }

  bool insert(std::tuple<K, V> kv) {
    K k = std::get<0>(kv);
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == empty_key) {
        if (pbbslib::CAS(&std::get<0>(table[h]), empty_key, k)) {
          std::get<1>(table[h]) = std::get<1>(kv);
          return true;
        }
      }
      if (std::get<0>(table[h]) == k) {
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  template <class F>
  bool insert_f(std::tuple<K, V> kv, const F& f) {
    K k = std::get<0>(kv);
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == empty_key) {
        if (pbbslib::CAS(&std::get<0>(table[h]), empty_key, k)) {
//          std::get<1>(table[h]) = std::get<1>(kv);
          f(&std::get<1>(table[h]), kv);
          return true;
        }
      }
      if (std::get<0>(table[h]) == k) {
        f(&std::get<1>(table[h]), kv);
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  bool insert_seq(std::tuple<K, V> kv) {
    K k = std::get<0>(kv);
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == empty_key) {
        table[h] = kv;
        return 1;
      }
      if (std::get<0>(table[h]) == k) {
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  bool insert_check(std::tuple<K, V> kv, bool* abort) {
    if (*abort) { return false; }
    K k = std::get<0>(kv);
    size_t h = firstIndex(k);
    size_t n_probes = 0;
    while (true) {
      if (std::get<0>(table[h]) == empty_key) {
        if (pbbslib::CAS(&std::get<0>(table[h]), empty_key, k)) {
          std::get<1>(table[h]) = std::get<1>(kv);
          return true;
        }
      }
      if (std::get<0>(table[h]) == k) {
        return false;
      }
      h = incrementIndex(h);
      n_probes++;
      if (n_probes > 10000) {
        *abort = true;
        break;
      }
    }
    return false;
  }

  void mark_seq(K k) {
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == empty_key) {
        return;
      }
      if (std::get<0>(table[h]) == k) {
        std::get<0>(table[h]) = empty_key - 1;
        return;
      }
      h = incrementIndex(h);
    }
  }

  bool contains(K k) const {
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == k) {
        return true;
      } else if (std::get<0>(table[h]) == empty_key) {
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  V find(K k, V default_value) const {
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == k) {
        return std::get<1>(table[h]);
      } else if (std::get<0>(table[h]) == empty_key) {
        return default_value;
      }
      h = incrementIndex(h);
    }
    return default_value;
  }

  template <class F>
  void map(F f) {
    parallel_for(0, m, [&] (size_t i) {
      if (std::get<0>(table[i]) != empty_key) {
        f(table[i]);
      }
    });
  }

  sequence<T> entries() const {
    auto pred = [&](const T& t) { return std::get<0>(t) != empty_key; };
    auto table_seq = parlay::make_slice(table, table + m);
    return parlay::filter(table_seq, pred);
  }

  void clear() {
    parallel_for(0, m, [&] (size_t i) { table[i] = empty; });
  }
};

template <class K, class V, class KeyHash>
inline sparse_table<K, V, KeyHash> make_sparse_table(size_t m,
                                                     std::tuple<K, V> empty,
                                                     KeyHash key_hash,
                                                     long space_mult=-1) {
  return sparse_table<K, V, KeyHash>(m, empty, key_hash, space_mult);
}

template <class K, class V, class KeyHash>
inline sparse_table<K, V, KeyHash> make_sparse_table(std::tuple<K, V>* tab,
                                                     size_t m,
                                                     std::tuple<K, V> empty,
                                                     KeyHash key_hash, bool clear=true) {
  return sparse_table<K, V, KeyHash>(m, empty, key_hash, tab, clear);
}

}  // namespace pbbslib
