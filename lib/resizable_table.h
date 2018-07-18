#pragma once

#include "sequence_ops.h"
#include "utilities.h"
#include <tuple>

// TODO: see if striding by an entire page reduces times more.
#define CACHE_STRIDE 128

static inline size_t hashToRange(const size_t& h, const size_t& mask) {
  return h & mask;
}
static inline size_t incrementIndex(const size_t& h, const size_t& mask) {
  return hashToRange(h + 1, mask);
}

template <class K, class V>
struct iter_kv {
  using T = std::tuple<K, V>;
  K k;
  size_t h;
  size_t mask;
  T* table;
  K empty_key;
  iter_kv(K _k, size_t _h, size_t _mask, T* _table, K _empty_key)
      : k(_k), h(_h), mask(_mask), table(_table), empty_key(_empty_key) {}

  // Finds the location of the first key
  bool init() {
    while (1) {
      if (std::get<0>(table[h]) == empty_key) {
        return false;
      } else if (std::get<0>(table[h]) == k) {
        return true;
      }
      h = incrementIndex(h, mask);
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
  T empty;
  K empty_key;
  T* table;
  bool alloc;
  KeyHash& key_hash;
  size_t* cts;
  size_t ne;

  static void clearA(T* A, long n, T kv) {
    parallel_for(size_t i = 0; i < n; i++) A[i] = kv;
  }

  inline size_t firstIndex(K& k) { return hashToRange(key_hash(k), mask); }

  void del() {
    if (alloc) {
      free(table);
      free(cts);
      alloc = false;
    }
  }

  void init_counts() {
    size_t workers = nworkers();
    cts = newA(size_t, CACHE_STRIDE * workers);
    for (size_t i = 0; i < workers; i++) {
      cts[i * CACHE_STRIDE] = 0;
    }
  }

  size_t update_nelms() {
    size_t workers = nworkers();
    for (size_t i = 0; i < workers; i++) {
      ne += cts[i * CACHE_STRIDE];
      cts[i * CACHE_STRIDE] = 0;
    }
  }

  resizable_table() : m(0), ne(0) {
    mask = 0;
    alloc = false;
    init_counts();
  }

  resizable_table(size_t _m, T _empty, KeyHash _key_hash, T* backing,
                  bool _alloc = false)
      : m(_m),
        mask(m - 1),
        ne(0),
        empty(_empty),
        empty_key(std::get<0>(empty)),
        table(backing),
        key_hash(_key_hash),
        alloc(_alloc) {
    clearA(table, m, empty);
    init_counts();
  }

  resizable_table(size_t _m, T _empty, KeyHash _key_hash)
      : m((size_t)1 << pbbs::log2_up((size_t)(1.1 * _m))),
        mask(m - 1),
        ne(0),
        empty(_empty),
        empty_key(std::get<0>(empty)),
        key_hash(_key_hash) {
    size_t line_size = 64;
    size_t bytes = ((m * sizeof(T)) / line_size + 1) * line_size;
    table = (T*)aligned_alloc(line_size, bytes);
    clearA(table, m, empty);
    init_counts();
    alloc = true;
  }

  void maybe_resize(size_t n_inc) {
    size_t nt = ne + n_inc;
    if (nt > (0.9 * m)) {
      size_t old_m = m;
      auto old_t = table;
      m = ((size_t)1 << pbbs::log2_up((size_t)(2 * nt)));
      if (m == old_m) {
        // should investigate
        return;
      }
      mask = m - 1;
      ne = 0;
      size_t line_size = 64;
      size_t bytes = ((m * sizeof(T)) / line_size + 1) * line_size;
      table = (T*)aligned_alloc(line_size, bytes);
      clearA(table, m, empty);
      parallel_for(size_t i = 0; i < old_m; i++) {
        if (std::get<0>(old_t[i]) != empty_key) {
          insert(old_t[i]);
        }
      }
      update_nelms();
      if (alloc) {
        free(old_t);
      }
      alloc = true;
    }
  }

  auto get_iter(K k) {
    size_t h = firstIndex(k);
    return iter_kv<K, V>(k, h, mask, table, empty_key);
  }

  bool insert(std::tuple<K, V> kv) {
    K& k = std::get<0>(kv);
    V& v = std::get<1>(kv);
    size_t h = firstIndex(k);
    while (1) {
      if (std::get<0>(table[h]) == empty_key && pbbs::CAS(&table[h], empty, kv)) {
        size_t wn = get_worker_num();
        cts[wn * CACHE_STRIDE]++;
        return 1;
      }
      if (std::get<0>(table[h]) == k && std::get<1>(table[h]) == v) {
        return false;
      }
      h = incrementIndex(h, mask);
    }
    return 0;
  }

  // TODO: finish this implementation.
  size_t get_label_set(K k, V* arr) {
    size_t h = firstIndex(k);
    size_t ct = 0;
    while (1) {
      if (std::get<0>(table[h]) == k) {
        arr[ct++] = std::get<1>(table[h]);
      } else if (std::get<0>(table[h]) == empty_key) {
        return ct;
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
    parallel_for(size_t i = 0; i < m; i++) {
      if (std::get<0>(table[i]) != empty_key) {
        f(table[i]);
      }
    }
  }

  auto entries() {
    T* out = newA(T, m);
    auto pred = [&](T& t) { return std::get<0>(t) != empty_key; };
    size_t new_m = pbbs::filterf(table, out, m, pred);
    return make_array_imap<T>(out, new_m);
  }

  void clear() {
    parallel_for(size_t i = 0; i < m; i++) { table[i] = empty; }
  }
};

template <class K, class V, class KeyHash>
auto make_resizable_table(size_t m, std::tuple<K, V> empty, KeyHash key_hash) {
  return resizable_table<K, V, KeyHash>(m, empty, key_hash);
}
