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
class tomb_table {
 public:
  using T = std::tuple<K, V>;

  size_t m;
  size_t mask;
  T empty;
  K empty_key;
  K tomb_key;
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
  inline size_t decrementIndex(size_t h) {return hashToRange(h-1);}
  inline bool lessIndex(size_t a, size_t b) {return 2 * hashToRange(a - b) > m;}

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

  tomb_table() : m(0) {
    mask = 0;
    alloc = false;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  tomb_table(size_t _m, T _empty, T _tomb_key, KeyHash _key_hash, long inp_space_mult=-1)
      : empty(_empty),
        empty_key(std::get<0>(empty)),
        tomb_key(_tomb_key),
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
  tomb_table(size_t _m, T _empty, T _tomb_key, KeyHash _key_hash, T* _tab, bool clear=true)
      : m(_m),
        mask(m - 1),
        empty(_empty),
        empty_key(std::get<0>(empty)),
        tomb_key(_tomb_key),
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


  // change to tombstone if found
  bool deleteVal(K k) {
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == k) {
        return pbbslib::CAS(&std::get<0>(table[h]), k, tomb_key);
      } else if (std::get<0>(table[h]) == empty_key) {//when delete, only stop when seeing empty
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
   }

  void maybe_resize(size_t n_inc, size_t ne) {
      size_t nt = ne + n_inc;
      if (nt > (0.9 * m)) {
        size_t old_m = m;
        auto old_t = table;
        m = ((size_t)1 << pbbslib::log2_up((size_t)(2 * nt)));
        if (m == old_m) {
          return;
        }
        mask = m - 1;
        ne = 0;
        size_t line_size = 64;
        size_t bytes = ((m * sizeof(T)) / line_size + 1) * line_size;
        table = (T*)pbbs::aligned_alloc(line_size, bytes);
        clearA(table, m, empty);
        parallel_for(0, old_m, [&] (size_t i) {
          if (std::get<0>(old_t[i]) != empty_key && std::get<0>(old_t[i]) != tomb_key) {
            insert(old_t[i]);
          }
        });
        // update_nelms();
        if (alloc) {
          pbbslib::free_array(old_t);
        }
        alloc = true;
      }
    }


  bool insert(std::tuple<K, V> kv) {
    K k = std::get<0>(kv);
    size_t h = firstIndex(k);
    while (true) {
      K prev_key = std::get<0>(table[h]);
      if (prev_key == empty_key || prev_key == tomb_key) {
        if (pbbslib::CAS(&std::get<0>(table[h]), prev_key, k)) {
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
      K prev_key = std::get<0>(table[h]);
      if (prev_key == empty_key || prev_key == tomb_key) {
        if (pbbslib::CAS(&std::get<0>(table[h]), prev_key, k)) {
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
      if (std::get<0>(table[h]) == empty_key || std::get<0>(table[h]) == tomb_key ) {
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

  // bool insert_check(std::tuple<K, V> kv, bool* abort) {
  //   if (*abort) { return false; }
  //   K k = std::get<0>(kv);
  //   size_t h = firstIndex(k);
  //   size_t n_probes = 0;
  //   while (true) {
  //     if (std::get<0>(table[h]) == empty_key) {
  //       if (pbbslib::CAS(&std::get<0>(table[h]), empty_key, k)) {
  //         std::get<1>(table[h]) = std::get<1>(kv);
  //         return true;
  //       }
  //     }
  //     if (std::get<0>(table[h]) == k) {
  //       return false;
  //     }
  //     h = incrementIndex(h);
  //     n_probes++;
  //     if (n_probes > 10000) {
  //       *abort = true;
  //       break;
  //     }
  //   }
  //   return false;
  // }

  // void mark_seq(K k) {
  //   size_t h = firstIndex(k);
  //   while (true) {
  //     if (std::get<0>(table[h]) == empty_key) {
  //       return;
  //     }
  //     if (std::get<0>(table[h]) == k) {
  //       std::get<0>(table[h]) = empty_key - 1;
  //       return;
  //     }
  //     h = incrementIndex(h);
  //   }
  // }

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

  sequence<T> entries() const {
    auto pred = [&](const T& t) { return std::get<0>(t) != empty_key && std::get<0>(t) != tomb_key; };
    auto table_seq = pbbslib::make_sequence<T>(table, m);
    return pbbslib::filter(table_seq, pred);
  }

  void clear() {
    parallel_for(0, m, [&] (size_t i) { table[i] = empty; });
  }
};


}  // namespace pbbslib