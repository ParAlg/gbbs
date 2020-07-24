  
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

template <class K, class KeyHash>
class sparse_set {
 public:
  using T = K;
  using V = K;

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
      pbbslib::free_array(table);
      alloc = false;
    }
  }

  // Does not copy elements in the current table; just resize the underlying
  // table.
  // Incoming must be a power of two
  void resize_no_copy(size_t incoming) {
    if (incoming > m) {
      if (alloc) {
        pbbslib::free_array(table);
      }
      std::cout << "# Resizing decrement table, was: " << m;
      m = incoming;
      mask = m - 1;
      table = pbbslib::new_array_no_init<T>(m);
      clearA(table, m, empty);
      alloc = true;
      std::cout << "#  is now: " << m << std::endl;
    }
  }

  sparse_set() : m(0) {
    mask = 0;
    alloc = false;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_set(size_t _m, T _empty, KeyHash _key_hash, long inp_space_mult=-1)
      : empty(_empty),
        empty_key(empty),
        key_hash(_key_hash) {
    double space_mult = 1.1;
    if (inp_space_mult != -1) space_mult = inp_space_mult;
    m = (size_t)1 << pbbslib::log2_up((size_t)(space_mult * _m) + 1);
    mask = m - 1;
    table = pbbslib::new_array_no_init<T>(m);
    clearA(table, m, empty);
    alloc = true;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_set(size_t _m, T _empty, KeyHash _key_hash, T* _tab, bool clear=true)
      : m(_m),
        mask(m - 1),
        empty(_empty),
        empty_key(empty),
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
      if (table[h] == k) {
        return h;
      }
      h = incrementIndex(h);
    }
  }

  bool insert(T kv) {
    K k = kv;
    size_t h = firstIndex(k);
    while (true) {
      if (table[h] == empty_key) {
        if (pbbslib::CAS(&table[h], empty_key, k)) {
        //   std::get<1>(table[h]) = std::get<1>(kv);
          return true;
        }
      }
      if (table[h] == k) {
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  template <class F>
  bool insert_f(T kv, const F& f) {
    K k = kv;
    size_t h = firstIndex(k);
    while (true) {
      if (table[h] == empty_key) {
        if (pbbslib::CAS(&table[h], empty_key, k)) {
//          std::get<1>(table[h]) = std::get<1>(kv);
          f(&table[h], kv);
          return true;
        }
      }
      if (table[h] == k) {
        f(&table[h], kv);
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  bool insert_seq(T kv) {
    K k = kv;
    size_t h = firstIndex(k);
    while (true) {
      if (table[h] == empty_key) {
        table[h] = kv;
        return 1;
      }
      if (table[h] == k) {
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  bool insert_check(T kv, bool* abort) {
    if (*abort) { return false; }
    K k = kv;
    size_t h = firstIndex(k);
    size_t n_probes = 0;
    while (true) {
      if (table[h] == empty_key) {
        if (pbbslib::CAS(&table[h], empty_key, k)) {
        //   std::get<1>(table[h]) = std::get<1>(kv);
          return true;
        }
      }
      if (table[h] == k) {
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
      if (table[h] == empty_key) {
        return;
      }
      if (table[h] == k) {
        table[h] = empty_key - 1;
        return;
      }
      h = incrementIndex(h);
    }
  }

  bool contains(K k) const {
    size_t h = firstIndex(k);
    while (true) {
      if (table[h] == k) {
        return true;
      } else if (table[h] == empty_key) {
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  V find(K k, V default_value) const {
    size_t h = firstIndex(k);
    while (true) {
      if (table[h] == k) {
        return table[h];
      } else if (table[h] == empty_key) {
        return default_value;
      }
      h = incrementIndex(h);
    }
    return default_value;
  }

  sequence<T> entries() const {
    auto pred = [&](const T& t) { return t != empty_key; };
    auto table_seq = pbbslib::make_sequence<T>(table, m);
    return pbbslib::filter(table_seq, pred);
  }

  void clear() {
    parallel_for(0, m, [&] (size_t i) { table[i] = empty; });
  }
};

template <class K, class KeyHash>
inline sparse_set<K, KeyHash> make_sparse_set(size_t m,
                                                     K empty,
                                                     KeyHash key_hash,
                                                     long space_mult=-1) {
  return sparse_set<K, KeyHash>(m, empty, key_hash, space_mult);
}

template <class K, class KeyHash>
inline sparse_set<K, KeyHash> make_sparse_set(K* tab,
                                                     size_t m,
                                                     K empty,
                                                     KeyHash key_hash, bool clear=true) {
  return sparse_set<K, KeyHash>(m, empty, key_hash, tab, clear);
}




}  // namespace pbbslib