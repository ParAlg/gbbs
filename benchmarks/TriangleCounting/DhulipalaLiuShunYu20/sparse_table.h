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
class sparse_table {
 public:
  using T = std::tuple<K, V>;
  using KT = K;

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

  bool not_empty(K k){
    return k != empty_key;
  }
  bool not_empty(T kv){
    return std::get<0>(kv) != empty_key;
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

  sparse_table() : m(0) {
    mask = 0;
    alloc = false;
  }

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
  sparse_table(size_t _m, T _empty, KeyHash _key_hash, long inp_space_mult=-1)
      : empty(_empty),
        empty_key(std::get<0>(empty)),
        tomb_key(std::get<0>(empty)),
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
  sparse_table(size_t _m, T _empty, KeyHash _key_hash, T* _tab, bool clear=true)
      : m(_m),
        mask(m - 1),
        empty(_empty),
        empty_key(std::get<0>(empty)),
        tomb_key(std::get<0>(empty)),
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

  // Pre-condition: k must be present in T.
  // do not support removing and inserting at the same time
  // needs to be more thoroughly tested
  // currently always returns true
  bool deleteVal(K v) {
    size_t i = firstIndex(v);
    int cmp = 1;

    // find first element less than or equal to v in priority order
    size_t j = i;
    T c = table[j];

    if (std::get<0>(c) == empty_key) return true;
      
    // find first location with priority less or equal to v's priority
    // while(c != empty && (cmp = hashStruct.cmp(v,hashStruct.getKey(c))) != 0) {
    while(std::get<0>(c) != empty_key && (cmp = key_hash.cmp(v,std::get<0>(c))) != 0) {
      j = incrementIndex(j);
      c = table[j];
    }
    cmp=(std::get<0>(c) == empty_key)?1:key_hash.cmp(v,std::get<0>(c));
    while (1) {
      // Invariants:
      //   v is the key that needs to be deleted
      //   j is our current index into TA
      //   if v appears in TA, then at least one copy must appear at or before j
      //   c = TA[j] at some previous time (could now be changed)
      //   i = h(v)
      //   cmp = compare v to key of c (1 if greater, 0 equal, -1 less)
      
      if (cmp != 0){//why doesn't the following work as the condition???
	//c==empty || hashStruct.cmp(v,hashStruct.getKey(c)) != 0) {
        // v does not match key of c, need to move down one and exit if
        // moving before h(v)
	if (j == i) return true;
  	j = decrementIndex(j);
  	c = table[j];
  	cmp = (std::get<0>(c) == empty_key) ? 1 : key_hash.cmp(v, std::get<0>(c));
      } else { // found v at location j (at least at some prior time)

  	// Find next available element to fill location j.
        // This is a little tricky since we need to skip over elements for
        // which the hash index is greater than j, and need to account for
        // things being moved downwards by others as we search.
        // Makes use of the fact that values in a cell can only decrease
        // during a delete phase as elements are moved from the right to left.
  	size_t jj = incrementIndex(j);
  	T x = table[jj];
  	while (std::get<0>(x) != empty_key && lessIndex(j, firstIndex(std::get<0>(x)))) {
  	  jj = incrementIndex(jj);
  	  x = table[jj];
  	}
  	size_t jjj = decrementIndex(jj);
  	while (jjj != j) {
  	  T y = table[jjj];
  	  if (std::get<0>(y) == empty_key || !lessIndex(j, firstIndex(std::get<0>(y)))) {
  	    x = y;
  	    jj = jjj;
  	  }
  	  jjj = decrementIndex(jjj);
  	}

  	// try to copy the the replacement element into j
  	// if (pbbslib::CAS(&std::get<0>(table[j]),std::get<0>(c),std::get<0>(x))) {
    //   std::get<1>(table[j]) = std::get<1>(x);
    if (pbbslib::CAS(&table[j],c,x)) {
          // swap was successful
          // if the replacement element was empty, we are done
  	  if (std::get<0>(x) == empty_key) return true;

  	  // Otherwise there are now two copies of the replacement element x
          // delete one copy (probably the original) by starting to look at jj.
          // Note that others can come along in the meantime and delete
          // one or both of them, but that is fine.
  	  v = std::get<0>(x);
  	  j = jj;
  	  i = firstIndex(v);
  	}
  	c = table[j];
  	cmp = (std::get<0>(c) == empty_key) ? 1 : key_hash.cmp(v, std::get<0>(c));
      }
    }
  }

  // Pre-condition: k must be present in T.
  // do not support updating and deleting/insert_f at the same time
  inline void updateSeq(K k, V val) {
    size_t h = idx(k);
    std::get<1>(table[h]) = val;
  }

  void maybe_resize(size_t nt) {
      if (nt > (0.9 * m) || nt < (m/4)) {
        size_t old_m = m;
        auto old_t = table;
        m = ((size_t)1 << pbbslib::log2_up((size_t)(1.2 * nt) + 1));
        if (m == old_m) {
          return;
        }
        mask = m - 1;
        // ne = 0;
        // size_t line_size = 128;
        // size_t bytes = ((m * sizeof(T)) / line_size + 1) * line_size;
        // table = (T*)pbbs::aligned_alloc(line_size, bytes);
        table = pbbslib::new_array_no_init<T>(m);
        clearA(table, m, empty);
        parallel_for(0, old_m, [&] (size_t i) {
          if (std::get<0>(old_t[i]) != empty_key) {
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

    void maybe_resize(size_t n_inc, size_t ne){
      size_t nt = ne + n_inc;
      maybe_resize(nt); 
    }


  bool insert(std::tuple<K, V> kv) {
    K k = std::get<0>(kv);
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == empty_key) {
        if (pbbslib::CAS(&std::get<0>(table[h]), empty_key, k)) {
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
        if (pbbslib::CAS(&std::get<0>(table[h]), empty_key, k)) {
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

  bool insert_check(std::tuple<K, V> kv, bool* abort) {
    if (*abort) { return false; }
    K k = std::get<0>(kv);
    size_t h = firstIndex(k);
    size_t n_probes = 0;
    while (true) {
      if (std::get<0>(table[h]) == empty_key) {
        if (pbbslib::CAS(&std::get<0>(table[h]), empty_key, k)) {
          std::get<1>(table[h]) = std::get<1>(kv);
          return true;
        }
      }
      if (std::get<0>(table[h]) == k) {
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
    auto pred = [&](const T& t) { return std::get<0>(t) != empty_key; };
    auto table_seq = pbbslib::make_sequence<T>(table, m);
    return pbbslib::filter(table_seq, pred);
  }

  struct getKey{
    getKey(){}
    K operator  () (const T& kv) {
      return std::get<0>(kv);
    }
  };
  sequence<K> keys() const {
    auto pred = [&](const K& k) { return k != empty_key; };
    auto table_seq = pbbslib::make_delayed<K>(m, getKey());
    return pbbslib::filter(table_seq, pred);
  }

  size_t entries_out(range<T> seq_out) const {
    auto pred = [&](const T& t) { return std::get<0>(t) != empty_key; };
    auto table_seq = pbbslib::make_sequence<T>(table, m);
    return pbbslib::filter_out(table_seq, seq_out, pred);
  }

  void clear() {
    parallel_for(0, m, [&] (size_t i) { table[i] = empty; });
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
                                                     KeyHash key_hash, bool clear=true) {
  return sparse_table<K, V, KeyHash>(m, empty, key_hash, tab, clear);
}

}  // namespace pbbslib