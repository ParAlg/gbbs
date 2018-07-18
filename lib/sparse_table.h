#pragma once

#include <tuple>

#include "sequence_ops.h"
#include "utilities.h"

using namespace std;

template <class K, class V, class KeyHash>
class sparse_table {
 public:
  using T = tuple<K, V>;

  size_t m;
  size_t mask;
  T empty;
  K empty_key;
  T* table;
  bool alloc;
  KeyHash& key_hash;

  static void clearA(T* A, long n, T kv) {
    parallel_for(size_t i = 0; i < n; i++) A[i] = kv;
  }

  inline size_t hashToRange(size_t h) { return h & mask; }
  inline size_t firstIndex(K& k) { return hashToRange(key_hash(k)); }
  inline size_t incrementIndex(size_t h) { return hashToRange(h + 1); }

  void del() {
    if (alloc) {
      free(table);
      alloc = false;
    }
  }

  sparse_table() : m(0) {
    mask = 0;
    alloc = false;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_table(size_t _m, T _empty, KeyHash _key_hash)
      : m((size_t)1 << pbbs::log2_up((size_t)(1.1 * _m))),
        mask(m - 1),
        empty(_empty),
        empty_key(get<0>(empty)),
        key_hash(_key_hash) {
    size_t line_size = 64;
    size_t bytes = ((m * sizeof(T)) / line_size + 1) * line_size;
    table = (T*)aligned_alloc(line_size, bytes);
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
        empty_key(get<0>(empty)),
        key_hash(_key_hash) {
    clearA(table, m, empty);
    alloc = false;
  }

  bool insert(tuple<K, V> kv) {
    K k = get<0>(kv);
    size_t h = firstIndex(k);
    while (1) {
      if (get<0>(table[h]) == empty_key) {
        if (pbbs::CAS(&std::get<0>(table[h]), empty_key, k)) {
          std::get<1>(table[h]) = std::get<1>(kv);
          return 1;
        }
      }
      if (std::get<0>(table[h]) == k) {
        return false;
      }
      h = incrementIndex(h);
    }
    return 0;
  }

  bool insert_seq(tuple<K, V> kv) {
    K k = get<0>(kv);
    size_t h = firstIndex(k);
    while (1) {
      if (get<0>(table[h]) == empty_key) {
        table[h] = kv;
        return 1;
      }
      if (std::get<0>(table[h]) == k) {
        return false;
      }
      h = incrementIndex(h);
    }
    return 0;
  }

  void mark_seq(K k) {
    size_t h = firstIndex(k);
    while (1) {
      if (get<0>(table[h]) == empty_key) {
        return;
      }
      if (std::get<0>(table[h]) == k) {
        get<0>(table[h]) = empty_key - 1;
        return;
      }
      h = incrementIndex(h);
    }
    return 0;
  }

  bool contains(K k) {
    size_t h = firstIndex(k);
    while (1) {
      if (get<0>(table[h]) == k) {
        return 1;
      } else if (get<0>(table[h]) == empty_key) {
        return 0;
      }
      h = incrementIndex(h);
    }
    return 0;
  }

  auto entries() {
    T* out = newA(T, m);
    auto pred = [&](T& t) { return get<0>(t) != empty_key; };
    size_t new_m = pbbs::filterf(table, out, m, pred);
    return make_array_imap<T>(out, new_m);
  }

  void clear() {
    parallel_for(size_t i = 0; i < m; i++) { table[i] = empty; }
  }
};

template <class K, class V, class KeyHash>
auto make_sparse_table(size_t m, tuple<K, V> empty, KeyHash key_hash) {
  return sparse_table<K, V, KeyHash>(m, empty, key_hash);
}

template <class K, class V, class KeyHash>
auto make_sparse_table(tuple<K, V>* tab, size_t m, tuple<K, V> empty,
                       KeyHash key_hash) {
  return sparse_table<K, V, KeyHash>(m, empty, key_hash, tab);
}
