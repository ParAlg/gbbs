// This code is based on the paper "Phase-Concurrent Hash Tables for
// Determinism" by Julian Shun and Guy Blelloch from SPAA 2014.
// Copyright (c) 2014 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights (to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANBILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#pragma once

#include <optional>

#include "gbbs/bridge.h"
#include "gbbs/macros.h"

namespace gbbs {

template <class K, class V>
class sequentialHT {
  // Note that if compiling without compiler-llvm-libcxx, std::tuple<..> must be
  // changed to absl::container_internal::CompressedTuple. Otherwise the empty
  // base class optimization that eliminates storage for empty structs will not
  // get performed, causing poor performance, and possibly bugs.
  typedef std::tuple<K, V> T;

 public:
  size_t m;
  intT mask;
  T empty;
  K max_key;
  T* table;

  inline size_t toRange(size_t h) { return h & mask; }
  inline size_t firstIndex(K v) { return toRange(parlay::hash64(v)); }
  inline size_t incrementIndex(size_t h) { return toRange(h + 1); }

  sequentialHT(T* _table, size_t size, float loadFactor,
               std::tuple<K, V> _empty)
      : m((size_t)1 << parlay::log2_up((size_t)(loadFactor * size + 1))),
        mask(m - 1),
        empty(_empty),
        max_key(std::get<0>(_empty)),
        table(_table) {}

  // m must be a power of two
  sequentialHT(T* _table, size_t _m, std::tuple<K, V> _empty)
      : m((size_t)_m),
        mask(m - 1),
        empty(_empty),
        max_key(std::get<0>(_empty)),
        table(_table) {}

  template <class M, class F>
  inline void insertF(std::tuple<K, M>& v, F& f) {
    K vKey = std::get<0>(v);
    size_t h = firstIndex(vKey);
    while (1) {
      auto k = std::get<0>(table[h]);
      if (k == max_key) {
        std::get<0>(table[h]) = vKey;
        V cur = std::get<1>(table[h]);
        std::get<1>(table[h]) = f(cur, v);
        return;
      } else if (k == vKey) {
        V cur = std::get<1>(table[h]);
        std::get<1>(table[h]) = f(cur, v);
        return;
      }
      h = incrementIndex(h);
    }
  }

  // V must support ++
  inline bool insertAdd(K& vKey) {
    size_t h = firstIndex(vKey);
    while (1) {
      auto k = std::get<0>(table[h]);
      if (k == max_key) {
        table[h] = std::make_tuple(vKey, 1);
        return true;
      } else if (k == vKey) {
        std::get<1>(table[h])++;
        return false;
      }
      h = incrementIndex(h);
    }
  }

  // V must support ++, T<1> must be numeric
  inline bool insertAdd(T& v) {
    const K& vKey = std::get<0>(v);
    size_t h = firstIndex(vKey);
    while (1) {
      auto k = std::get<0>(table[h]);
      if (k == max_key) {
        table[h] = std::make_tuple(vKey, 1);
        return true;
      } else if (k == vKey) {
        std::get<1>(table[h]) += std::get<1>(v);
        return false;
      }
      h = incrementIndex(h);
    }
  }

  inline T find(K& v) {
    size_t h = firstIndex(v);
    T c = table[h];
    while (1) {
      if (std::get<0>(c) == max_key) {
        return empty;
      } else if (std::get<0>(c) == v) {
        return c;
      }
      h = incrementIndex(h);
      c = table[h];
    }
  }

  // F : KV -> E
  template <class E, class F>
  inline size_t compactInto(F& f, E* Out) {
    size_t k = 0;
    for (size_t i = 0; i < m; i++) {
      auto kv = table[i];
      auto key = std::get<0>(kv);
      if (key != max_key) {
        table[i] = empty;
        std::optional<E> value = f(kv);
        if (value.has_value()) {
          Out[k++] = *value;
        }
      }
    }
    return k;
  }

  // F : KV -> std::optional<std::tuple<uintE,?>>
  template <class F>
  inline size_t compactIntoSelf(F& f) {
    size_t k = 0;
    for (size_t i = 0; i < m; i++) {
      auto kv = table[i];
      auto key = std::get<0>(kv);
      if (key != max_key) {
        table[i] = empty;
        auto value = f(kv);
        if (value.has_value()) {
          table[k++] = *value;
        }
      }
    }
    return k;
  }
};

}  // namespace gbbs
