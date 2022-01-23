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

#include "gbbs/bridge.h"

namespace gbbs {

template <class K, class V>
class sparse_additive_map {
 public:
  using T = std::tuple<K, V>;

  size_t m;
  size_t mask;
  T empty;
  K empty_key;
  T* table;
  bool alloc;

  static void clearA(T* A, long n, T kv) {
    parallel_for(0, n, [&](size_t i) { A[i] = kv; });
  }

  inline size_t hashToRange(size_t h) { return h & mask; }
  inline size_t firstIndex(K& k) { return hashToRange(parlay::hash64(k)); }
  inline size_t incrementIndex(size_t h) { return hashToRange(h + 1); }

  void del() {
    if (alloc) {
      gbbs::free_array(table);
      alloc = false;
    }
  }

  sparse_additive_map() : m(0) {
    mask = 0;
    alloc = false;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_additive_map(size_t _m, T _empty)
      : m((size_t)1 << parlay::log2_up((size_t)(1.1 * _m))),
        mask(m - 1),
        empty(_empty),
        empty_key(std::get<0>(empty)) {
    size_t line_size = 64;
    size_t bytes = ((m * sizeof(T)) / line_size + 1) * line_size;
    table = (T*)aligned_alloc(line_size, bytes);
    clearA(table, m, empty);
    alloc = true;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_additive_map(size_t _m, T _empty, T* _tab)
      : m(_m),
        mask(m - 1),
        table(_tab),
        empty(_empty),
        empty_key(std::get<0>(empty)) {
    clearA(table, m, empty);
    alloc = false;
  }

  bool insert(std::tuple<K, V> kv) {
    K k = std::get<0>(kv);
    V v = std::get<1>(kv);
    size_t h = firstIndex(k);
    while (1) {
      if (std::get<0>(table[h]) == empty_key) {
        if (gbbs::atomic_compare_and_swap(&std::get<0>(table[h]), empty_key,
                                          k)) {
          gbbs::write_add(&std::get<1>(table[h]), v);
          return 1;
        }
      }
      if (std::get<0>(table[h]) == k) {
        gbbs::write_add(&std::get<1>(table[h]), v);
        return false;
      }
      h = incrementIndex(h);
    }
    return 0;
  }

  bool contains(K k) {
    size_t h = firstIndex(k);
    while (1) {
      if (std::get<0>(table[h]) == k) {
        return 1;
      } else if (std::get<0>(table[h]) == empty_key) {
        return 0;
      }
      h = incrementIndex(h);
    }
    return 0;
  }

  auto entries() {
    auto pred = [&](const T& t) { return std::get<0>(t) != empty_key; };
    auto table_seq = gbbs::make_slice<T>(table, m);
    return parlay::filter(table_seq, pred);
  }

  void clear() {
    parallel_for(0, m, [&](size_t i) { table[i] = empty; });
  }
};

}  // namespace gbbs
