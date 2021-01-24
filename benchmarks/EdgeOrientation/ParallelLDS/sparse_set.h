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

namespace gbbs {

template <class K>
class sparse_set {
 public:

  static constexpr K kEmptyKey = std::numeric_limits<K>::max();
  static constexpr double kSpaceMult = 1.25;

  uintE mask;  // table.size() - 1
  uintE num_elms;
  parlay::sequence<K> table;  // table.size() is a power of two

  size_t size() const {
    return table.size();
  }

  static void clearA(K* A, size_t n, K k) {
    parallel_for(0, n, [&] (size_t i) { A[i] = k; });
  }

  inline size_t hashToRange(size_t h) const { return h & mask; }
  inline size_t firstIndex(K& k) const { return hashToRange(parlay::hash32(k)); }
  inline size_t incrementIndex(size_t h) const { return hashToRange(h + 1); }

  // Updates num_elms by the number of incoming elements (guaranteed not to be
  // in the table).
  void resize(size_t incoming) {
    size_t total = num_elms + incoming;
    // std::cout << "total = " << total << " n_elms = " << num_elms << " incoming = " << incoming << std::endl;
    if (total * kSpaceMult >= table.size()) {
      size_t new_size = (1 << parlay::log2_up((size_t)(kSpaceMult * total) + 1));
      // std::cout << "new_size = " << new_size << std::endl;
      auto new_table = parlay::sequence<K>(new_size, kEmptyKey);
      auto old_table = std::move(table);
      table = std::move(new_table);
      // std::cout << "old_table_size = " << old_table.size() << std::endl;
      parallel_for(0, old_table.size(), [&] (size_t i) {
        if (old_table[i] != kEmptyKey) {
          insert(old_table[i]);
        }
      });
    }
    num_elms += incoming;
  }

  sparse_set() : mask(0), num_elms(0) {}

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_set(size_t _m, long inp_space_mult=-1) {
    double space_mult = 1.1;
    if (inp_space_mult != -1) space_mult = inp_space_mult;
    auto m = (size_t)1 << parlay::log2_up((size_t)(space_mult * _m) + 1);
    mask = m - 1;
    table = parlay::sequence<K>(m, kEmptyKey);
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

  bool insert(K k) {
    size_t h = firstIndex(k);
    while (true) {
      if (table[h] == kEmptyKey) {
        if (pbbslib::CAS(&table[h], kEmptyKey, k)) {
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

  bool contains(K k) const {
    if (table.size() > 0) {
      size_t h = firstIndex(k);
      while (true) {
        if (table[h] == k) {
          return true;
        } else if (table[h] == kEmptyKey) {
          return false;
        }
        h = incrementIndex(h);
      }
    }
    return false;
  }

  sequence<K> entries() const {
    auto pred = [&](const K& k) { return k != kEmptyKey; };
    return parlay::filter(parlay::make_slice(table), pred);
  }

  void clear() {
    parallel_for(0, table.size(), [&] (size_t i) { table[i] = kEmptyKey; });
  }
};

}  // namespace gbbs
