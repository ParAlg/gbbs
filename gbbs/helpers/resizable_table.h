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
#include <unordered_set>
#include <vector>

#include "gbbs/bridge.h"

namespace gbbs {
// TODO: see if striding by an entire page improves times further.
constexpr size_t kResizableTableCacheLineSz = 128;

inline size_t hashToRange(const size_t& h, const size_t& mask) {
  return h & mask;
}
inline size_t incrementIndex(const size_t& h, const size_t& mask) {
  return hashToRange(h + 1, mask);
}

template <class K, class V>
struct iter_kv {
  using T = std::tuple<K, V>;
  K k;
  size_t h;
  size_t mask;
  // size_t num_probes;  // For debugging.
  T* table;
  K empty_key;
  iter_kv(K _k, size_t _h, size_t _mask, T* _table, K _empty_key)
      : k(_k),
        h(_h),
        mask(_mask),
        // num_probes(0),
        table(_table),
        empty_key(_empty_key) {}

  // Finds the location of the first key
  bool init() {
    while (1) {
      if (std::get<0>(table[h]) == empty_key) {
        return false;
      } else if (std::get<0>(table[h]) == k) {
        return true;
      }
      h = incrementIndex(h, mask);
      // num_probes++;
    }
    return false;
  }

  // Probes until we've found another key
  bool has_next() {
    h = incrementIndex(h, mask);
    while (1) {
      if (std::get<0>(table[h]) == empty_key) {
        return false;
      } else if (std::get<0>(table[h]) == k) {
        return true;
      }
      h = incrementIndex(h, mask);
      // num_probes++;
    }
    return false;
  }

  T next() { return table[h]; }
};

template <class K, class V, class KeyHash>
class resizable_table {
 public:
  using T = std::tuple<K, V>;

  size_t m;
  size_t mask;
  size_t ne;
  T empty;
  K empty_key;
  sequence<T> table;
  KeyHash key_hash;
  sequence<size_t> cts;

  inline size_t firstIndex(K& k) { return hashToRange(key_hash(k), mask); }

  void init_counts() {
    size_t workers = num_workers();
    cts = sequence<size_t>::uninitialized(kResizableTableCacheLineSz * workers);
    for (size_t i = 0; i < workers; i++) {
      cts[i * kResizableTableCacheLineSz] = 0;
    }
  }

  void update_nelms() {
    size_t workers = num_workers();
    for (size_t i = 0; i < workers; i++) {
      ne += cts[i * kResizableTableCacheLineSz];
      cts[i * kResizableTableCacheLineSz] = 0;
    }
  }

  resizable_table() : m(0), ne(0) {
    mask = 0;
    init_counts();
  }

  resizable_table(size_t _m, T _empty, KeyHash _key_hash)
      : m((size_t)1 << parlay::log2_up((size_t)(1.1 * _m))),
        mask(m - 1),
        ne(0),
        empty(_empty),
        empty_key(std::get<0>(empty)),
        key_hash(_key_hash) {
    table = sequence<T>::uninitialized(m);
    clear();
    init_counts();
  }

  void analyze() {
    std::vector<size_t> probe_lengths;
    size_t chain_length = 0;
    for (size_t i = 0; i < m; i++) {
      if (std::get<0>(table[i]) == empty_key) {
        if (chain_length > 0) {
          probe_lengths.emplace_back(chain_length);
        }
        chain_length = 0;
      } else {
        chain_length++;
      }
    }

    std::sort(probe_lengths.begin(), probe_lengths.end());
    std::cout << "Analyzed table, m = " << m << " ne = " << ne
              << " num_probes = " << probe_lengths.size() << std::endl;
    for (long k = probe_lengths.size() - 1; k > 0; k--) {
      std::cout << probe_lengths[k] << " ";
      if (probe_lengths.size() - k > 250) break;
    }
    std::cout << std::endl << "End of table analysis" << std::endl;
  }

  void maybe_resize(size_t n_inc) {
    size_t nt = ne + n_inc;
    if (nt > (0.25 * m)) {
      size_t old_m = m;
      auto old_t = std::move(table);
      m = ((size_t)1 << parlay::log2_up((size_t)(10 * nt)));
      if (m == old_m) {
        return;
      }
      mask = m - 1;
      ne = 0;
      table = sequence<T>::uninitialized(m);
      clear();
      parallel_for(0, old_m, [&](size_t i) {
        if (std::get<0>(old_t[i]) != empty_key) {
          insert(old_t[i]);
        }
      });
      update_nelms();
    }
  }

  iter_kv<K, V> get_iter(K k) {
    size_t h = firstIndex(k);
    return iter_kv<K, V>(k, h, mask, table.begin(), empty_key);
  }

  bool insert(std::tuple<K, V> kv) {
    K& k = std::get<0>(kv);
    V& v = std::get<1>(kv);
    size_t h = firstIndex(k);
    while (1) {
      if (std::get<0>(table[h]) == empty_key &&
          gbbs::atomic_compare_and_swap(&table[h], empty, kv)) {
        size_t wn = worker_id();
        cts[wn * kResizableTableCacheLineSz]++;
        return 1;
      }
      if (std::get<0>(table[h]) == k && std::get<1>(table[h]) == v) {
        return false;
      }
      h = incrementIndex(h, mask);
    }
    return 0;
  }

  bool insert_seq(std::tuple<K, V> kv) {
    K& k = std::get<0>(kv);
    V& v = std::get<1>(kv);
    size_t h = firstIndex(k);
    while (1) {
      if (std::get<0>(table[h]) == empty_key) {
        table[h] = kv;
        return 1;
      }
      h = incrementIndex(h, mask);
    }
    return 0;
  }

  size_t num_appearances(K k) {
    size_t h = firstIndex(k);
    size_t ct = 0;
    while (1) {
      if (std::get<0>(table[h]) == k) {
        ct++;
      } else if (std::get<0>(table[h]) == empty_key) {
        return ct;
      }
      h = incrementIndex(h, mask);
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
      h = incrementIndex(h, mask);
    }
    return 0;
  }

  bool contains(K k, V v) {
    size_t h = firstIndex(k);
    while (1) {
      if (std::get<0>(table[h]) == k && std::get<1>(table[h]) == v) {
        return 1;
      } else if (std::get<0>(table[h]) == empty_key) {
        return 0;
      }
      h = incrementIndex(h, mask);
    }
    return 0;
  }

  template <class F>
  void map(F& f) {
    parallel_for(0, m, [&](size_t i) {
      if (std::get<0>(table[i]) != empty_key) {
        f(table[i]);
      }
    });
  }

  sequence<T> entries() {
    auto pred = [&](T& t) { return std::get<0>(t) != empty_key; };
    auto table_seq = gbbs::make_slice<T>(table, m);
    return parlay::filter(table_seq, pred);
  }

  void clear() {
    parallel_for(0, m, [&](size_t i) { table[i] = empty; });
  }
};

template <class K, class V, class KeyHash>
inline resizable_table<K, V, KeyHash> make_resizable_table(
    size_t m, std::tuple<K, V> empty, KeyHash key_hash) {
  return resizable_table<K, V, KeyHash>(m, empty, key_hash);
}

}  // namespace gbbs
