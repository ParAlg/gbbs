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

#include "lib/sequence_ops.h"
#include "lib/extra_sequence_ops.h"
#include "lib/utilities.h"
#include "lib/macros.h"
#include <cassert>

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
  KeyHash& key_hash;

  size_t size() {
    return m;
  }

  static void clearA(T* A, long n, T kv) {
    par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i)
                    { A[i] = kv; });
  }

  inline size_t hashToRange(size_t h) { return h & mask; }
  inline size_t firstIndex(K& k) { return hashToRange(key_hash(k)); }
  inline size_t incrementIndex(size_t h) { return hashToRange(h + 1); }

  void del() {
    if (alloc) {
      pbbs::free_array(table);
      alloc = false;
    }
  }

  sparse_table() : m(0) {
    mask = 0;
    alloc = false;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_table(size_t _m, T _empty, KeyHash _key_hash, long space_mult=-1)
      : empty(_empty),
        empty_key(std::get<0>(empty)),
        key_hash(_key_hash) {
    if (space_mult == -1) space_mult = 1.1;
    m = (size_t)1 << pbbs::log2_up((size_t)(space_mult * _m));
    mask = m - 1;
    table = pbbs::new_array_no_init<T>(m);
    std::cout << "T size is " << sizeof(T) << std::endl;
    clearA(table, m, empty);
    alloc = true;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_table(size_t _m, T _empty, KeyHash _key_hash, T* _tab)
      : m(_m),
        mask(m - 1),
        table(_tab),
        empty(_empty),
        empty_key(std::get<0>(empty)),
        key_hash(_key_hash) {
    clearA(table, m, empty);
    alloc = false;
  }

  // Pre-condition: k must be present in T.
  inline size_t idx(K k) {
    size_t h = firstIndex(k);
    while (1) {
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
        if (pbbs::CAS(&std::get<0>(table[h]), empty_key, k)) {
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
        if (pbbs::CAS(&std::get<0>(table[h]), empty_key, k)) {
          std::get<1>(table[h]) = std::get<1>(kv);
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

  bool contains(K k) {
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

  V find(K k, V default_value) {
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

  sequence<T> entries() {
    auto pred = [&](T& t) { return std::get<0>(t) != empty_key; };
    auto table_seq = make_sequence<T>(table, m);
    return pbbs::filter(table_seq, pred);
  }

  void clear() {
    par_for(0, m, 2048, [&] (size_t i) { table[i] = empty; });
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
                                                     KeyHash key_hash) {
  return sparse_table<K, V, KeyHash>(m, empty, key_hash, tab);
}
