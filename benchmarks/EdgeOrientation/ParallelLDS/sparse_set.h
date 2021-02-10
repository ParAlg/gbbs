#pragma once

#include <tuple>
#include <cassert>

#include "gbbs/bridge.h"

namespace gbbs {

template <class K>
class sparse_set {
 public:

  static constexpr K kEmptyKey = std::numeric_limits<K>::max();
  static constexpr K kTombstone = std::numeric_limits<K>::max() - 1;
  static constexpr double kSpaceMult = 1.25;

  uintE mask;  // table.size() - 1
  uintE elms_in_table;
  uintE tombstones_in_table;
  parlay::sequence<K> table_seq;  // table.size() is a power of two
  K* table;  // = table_seq.begin()

  size_t size() const {
    return table_seq.size();
  }

  size_t num_elms() const {
    return elms_in_table;
  }

  static void clearA(K* A, size_t n, K k) {
    parallel_for(0, n, [&] (size_t i) { A[i] = k; });
  }

  static inline bool valid(K k) {
    return (k != kEmptyKey) && (k != kTombstone);
  }

  inline size_t hashToRange(size_t h) const { return h & mask; }
  inline size_t firstIndex(K& k) const { return hashToRange(parlay::hash32(k)); }
  inline size_t incrementIndex(size_t h) const { return hashToRange(h + 1); }

  // Updates elms_in_table by the number of incoming elements (guaranteed not to be
  // in the table).
  void resize(size_t incoming) {
    size_t total = elms_in_table + incoming;
    // std::cout << "total = " << total << " n_elms = " << elms_in_table << " incoming = " << incoming << std::endl;
    if (total * kSpaceMult + tombstones_in_table >= table_seq.size()) {
      tombstones_in_table = 0;
    //if (total * kSpaceMult >= table_seq.size()) {
      size_t new_size = (1 << parlay::log2_up((size_t)(kSpaceMult * total) + 1));
      new_size = std::max(new_size, (size_t)8);
      // std::cout << "new_size = " << new_size << std::endl;
      auto new_table = parlay::sequence<K>(new_size, kEmptyKey);
      auto old_table = std::move(table_seq);
      table_seq = std::move(new_table);
      table = table_seq.begin();
      mask = table_seq.size() - 1;
      // std::cout << "old_table_size = " << old_table.size() << std::endl;
      auto old_tab = old_table.begin();
      parallel_for(0, old_table.size(), [&] (size_t i) {
        if (valid(old_tab[i])) {
          insert(old_tab[i]);
        }
      });
    }
    elms_in_table += incoming;
  }

  void resize_down(size_t removed) {
    if (removed == 0) return;
    size_t total = elms_in_table - removed;
    tombstones_in_table += removed;

    if (total == 0) {
      auto old_table = std::move(table_seq);
      table_seq = parlay::sequence<K>();
      mask = 0;
      tombstones_in_table = 0;
    } else if (total * kSpaceMult <= table_seq.size() / 2) {
      size_t new_size = (1 << parlay::log2_up((size_t)(kSpaceMult * total)));
      new_size = std::max(new_size, (size_t)8);  // some minimal size
      auto new_table = parlay::sequence<K>(new_size, kEmptyKey);
      auto old_table = std::move(table_seq);
      table_seq = std::move(new_table);
      table = table_seq.begin();
      mask = table_seq.size() - 1;
      auto old_tab = old_table.begin();
      parallel_for(0, old_table.size(), [&] (size_t i) {
        if (valid(old_tab[i])) {
          insert(old_tab[i]);
        }
      });
      tombstones_in_table = 0;
    }
    elms_in_table -= removed;

    if ((tombstones_in_table + elms_in_table) > 0.8*elms_in_table) {
      auto elts = entries();
      clearA(table_seq.begin(), table_seq.size(), kEmptyKey);
      parallel_for(0, elts.size(), [&] (size_t i) {
        insert(elts[i]);
      });
      tombstones_in_table = 0;
    }

  }

  sparse_set() : mask(0), elms_in_table(0), tombstones_in_table(0) {}

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_set(size_t _m, long inp_space_mult=-1) {
    double space_mult = 1.1;
    if (inp_space_mult != -1) space_mult = inp_space_mult;
    auto m = (size_t)1 << parlay::log2_up((size_t)(space_mult * _m) + 1);
    m = std::max(m, (size_t)8);
    mask = m - 1;
    table_seq = parlay::sequence<K>(m, kEmptyKey);
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

  // Can be called concurrently. Must ensure that the caller has called resize
  // before concurrent insertions to update num_elms appropriately.
  bool insert(K k) {
    size_t h = firstIndex(k);
    while (true) {
      auto read = table[h];
      if (read == kEmptyKey || read == kTombstone) {
        if (pbbslib::CAS(&table[h], read, k)) {
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

  bool remove(K k) {
    size_t h = firstIndex(k);
    while (true) {
      if (table[h] == kEmptyKey) {
        return false;
      }
      if (table[h] == k) {
        bool succ = pbbslib::CAS(&table[h], k, kTombstone);
        return succ;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  // Requires that no element in s is currently in the structure (s has been
  // prefiltered using contains()).
  template <class Seq>
  void append(Seq& s) {
    resize(s.size());  // resize updates num_elms.
    parallel_for(0, s.size(), [&] (size_t i) {
      insert(s[i]);
    });
  }

  bool contains(K k) const {
    if (table_seq.size() > 0) {
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
    auto pred = [&](const K& k) { return valid(k); };
    return parlay::filter(parlay::make_slice(table_seq), pred);
  }

  void clear() {
    parallel_for(0, table_seq.size(), [&] (size_t i) { table[i] = kEmptyKey; });
  }
};

}  // namespace gbbs
