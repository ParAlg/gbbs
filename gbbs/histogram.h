// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010-2016 Guy Blelloch, Laxman Dhulipala, and the PBBS team
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

#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <assert.h>

#include "pbbslib/counting_sort_no_transpose.h"
#include "bridge.h"
#include "macros.h"
#include "sequential_ht.h"

namespace gbbs {

// Tunable parameters
constexpr const size_t _hist_max_buckets = 1024;
constexpr const size_t _hist_seq_threshold = 4096;

template <typename E, class B>
struct get_bucket {
  std::tuple<E, int>* hash_table;
  size_t table_mask;
  size_t low_mask;
  size_t bucket_mask;
  int num_buckets;
  int k;
  B& I;

  std::tuple<E*, int> heavy_hitters(size_t n, size_t count) {
    E* sample = pbbslib::new_array_no_init<E>(count);
    for (size_t i = 0; i < count; i++) {
      sample[i] = I[pbbslib::hash32(i) % n];
    }
    std::sort(sample, sample + count);

    // only keep those with at least three copies
    int j = 0;
    int c = 0;
    for (size_t i = 1; i < count; i++) {
      if (sample[i] == sample[i - 1]) {
        if (c++ == 1) {
          sample[j++] = sample[i];
        }
      } else {
        c = 0;
      }
    }
    return std::make_tuple(sample, j);
  }

  std::tuple<E, int>* make_hash_table(E* entries, size_t n, size_t table_size,
                                      size_t _table_mask) {
    using ttype = std::tuple<E, int>;
    auto table = pbbslib::new_array_no_init<ttype>(table_size);
    for (size_t i = 0; i < table_size; i++) table[i] = std::make_pair(0, -1);
    size_t n_distinct = 0;
    for (size_t i = 0; i < n; i++) {
      size_t h = pbbslib::hash32(entries[i]) & _table_mask;
      while (std::get<1>(table[h]) != -1) {
        h = (h + 1) & _table_mask;
      }
      table[h] = std::make_pair(entries[i], n_distinct++);
    }
    return table;
  }

  get_bucket(B& A, size_t n, size_t bits) : I(A) {
    num_buckets = 1 << bits;
    bucket_mask = num_buckets - 1;
    low_mask = ~((size_t)15);
    int count = 2 * num_buckets;
    int table_size = 4 * count;  // tune
    table_mask = table_size - 1;

    E* sample;
    std::tie(sample, k) = heavy_hitters(n, count);

    if (k > 0) {
      hash_table = make_hash_table(sample, k, table_size, table_mask);
    }
    pbbslib::free_array(sample);
  }

  ~get_bucket() {
    if (k > 0) {
      pbbslib::free_array(hash_table);
    }
  }

  size_t operator()(size_t i) {
    if (k > 0) {
      size_t h = pbbslib::hash32(I[i]) & table_mask;
      while (true) {
        E elm;
        int ct;
        std::tie(elm, ct) = hash_table[h];
        if (ct == -1) {
          break;
        }
        if (elm == I[i] && ct != -1) {
          return ct + num_buckets;
        }
        h = (h + 1) & table_mask;
      }
    }
    return pbbslib::hash32(I[i] & low_mask) & bucket_mask;
  }
};

template <class K, class V>
struct hist_table {
  using KV = std::tuple<K, V>;
  KV empty;
  KV* table;
  size_t size;
  hist_table(KV _empty, size_t _size) : empty(_empty), size(_size) {
    table = pbbslib::new_array_no_init<KV>(size);
    par_for(0, size, 2048, [&] (size_t i) { table[i] = empty; });
  }
  hist_table() {}

  void resize(size_t req_size) {
    if (req_size > size) {
      size_t rounded_size = (1L << pbbslib::log2_up<size_t>(req_size));
      pbbslib::free_array(table);
      table = pbbslib::new_array_no_init<KV>(rounded_size);
      size = rounded_size;
      par_for(0, size, 2048, [&] (size_t i) { table[i] = empty; });
      debug(std::cout << "resized to: " << size << "\n";);
    }
  }

  void del() {
    if (table) {
      pbbslib::free_array(table);
    }
  }
};



// Parallelizes across buckets, but does not use light/heavy buckets
template <class O, class K, class V, class A, class Apply>
inline std::pair<size_t, O*> histogram_medium(A& get_key, size_t n,
                                              Apply& apply_f,
                                              hist_table<K, V>& ht) {
  using KV = std::tuple<K, V>;
  size_t sqrt = (size_t)ceil(pow(n, 0.5));
  size_t num_buckets = (size_t)(n < 20000000) ? (sqrt / 5) : sqrt;

  num_buckets = std::max(1 << pbbslib::log2_up(num_buckets), 1);
  num_buckets = std::min(num_buckets, _hist_max_buckets);

  // (1) count-sort based on bucket
  size_t low_mask = ~((size_t)15);
  size_t bucket_mask = num_buckets - 1;
  auto gb = [&](uintE i) {
    return pbbslib::hash32(get_key[i] & low_mask) & bucket_mask;
  };

  K* elms;
  size_t* counts;
  size_t num_blocks;
  if (num_buckets <= 256) {
    std::tie(elms, counts, num_blocks) =
        pbbslib::_count_sort<uint8_t, size_t, K>(get_key, gb, n, (uintE)num_buckets);
  } else {
    std::tie(elms, counts, num_blocks) =
        pbbslib::_count_sort<uint16_t, size_t, K>(get_key, gb, n, (uintE)num_buckets);
  }
  size_t block_size = ((n - 1) / num_blocks) + 1;

#define S_STRIDE 64
  size_t* bkt_counts = pbbslib::new_array_no_init<size_t>(num_buckets * S_STRIDE);
  par_for(0, num_buckets, 1, [&] (size_t i) {
    bkt_counts[i * S_STRIDE] = 0;
    if (i == (num_buckets - 1)) {
      size_t ct = 0;
      for (size_t j = 0; j < num_blocks; j++) {
        size_t start = std::min(j * block_size, n);
        size_t end = std::min(start + block_size, n);
        ct += (end - start) - counts[j * num_buckets + i];
      }
      bkt_counts[i * S_STRIDE] = ct;
    } else {
      size_t ct = 0;
      for (size_t j = 0; j < num_blocks; j++) {
        ct += counts[j * num_buckets + i + 1] - counts[j * num_buckets + i];
      }
      bkt_counts[i * S_STRIDE] = ct;
    }
  });

  sequence<size_t> out_offs = sequence<size_t>(num_buckets + 1);
  sequence<size_t> ht_offs = sequence<size_t>(num_buckets + 1);

  // (2) process each bucket, compute the size of each HT and scan (seq)
  ht_offs[0] = 0;
  size_t min_size = std::numeric_limits<size_t>::max();
  size_t max_size = 0;
  for (size_t i = 0; i < num_buckets; i++) {
    size_t size = bkt_counts[i * S_STRIDE];
    size_t ht_size = 0;
    if (size > 0) {
      ht_size = 1 << pbbslib::log2_up((intT)(size + 1));
    }
    ht_offs[i + 1] = ht_offs[i] + ht_size;
    if (size < min_size) {
      min_size = size;
    }
    if (size > max_size) {
      max_size = size;
    }
  }

  ht.resize(ht_offs[num_buckets]);
  KV* table = ht.table;
  auto empty = ht.empty;

  // (3) insert elms into per-bucket hash table (par)
  {
    auto inner_for = [&](size_t i) {
      size_t ht_start = ht_offs[i];
      size_t ht_size = ht_offs[i + 1] - ht_start;

      size_t k = 0;
      if (ht_size > 0) {
        KV* my_ht = &(table[ht_start]);
        sequentialHT<K, V> S(my_ht, ht_size, empty);

        for (size_t j = 0; j < num_blocks; j++) {
          size_t start = std::min(j * block_size, n);
          size_t end = std::min(start + block_size, n);
          size_t ct = 0;
          size_t off = 0;
          if (i == (num_buckets - 1)) {
            off = counts[j * num_buckets + i];
            ct = (end - start) - off;
          } else {
            off = counts[j * num_buckets + i];
            ct = counts[j * num_buckets + i + 1] - off;
          }
          off += start;
          for (size_t l = 0; l < ct; l++) {
            K a = elms[off + l];
            S.insertAdd(a);
          }
        }

        k = S.compactIntoSelf(apply_f);
      }
      out_offs[i] = k;
    };
    par_for(0, num_buckets, 1, [&] (size_t i) { inner_for(i); });
  }

  // (4) scan
  size_t ct = 0;
  for (size_t i = 0; i < num_buckets; i++) {
    size_t s = ct;
    ct += out_offs[i];
    out_offs[i] = s;
  }
  out_offs[num_buckets] = ct;
  uintT num_distinct = ct;

  O* res = pbbslib::new_array_no_init<O>(ct);

  // (5) map compacted hts to output, clear hts
  par_for(0, num_buckets, 1, [&] (size_t i) {
    size_t o = out_offs[i];
    size_t k = out_offs[(i + 1)] - o;

    size_t ht_start = ht_offs[i];
    size_t ht_size = ht_offs[i + 1] - ht_start;

    if (ht_size > 0) {
      KV* my_ht = &(table[ht_start]);

      for (size_t j = 0; j < k; j++) {
        res[o + j] = my_ht[j];
        my_ht[j] = empty;
      }
    }
  });

  pbbslib::free_array(elms);
  pbbslib::free_array(counts);
  pbbslib::free_array(bkt_counts);
  return std::make_pair(num_distinct, res);
}

// Applies light/heavy buckets
template <class O, class K, class V, class A, class Apply>
inline std::pair<size_t, O*> histogram(A& get_key, size_t n, Apply& apply_f,
                                       hist_table<K, V>& ht) {
  using KV = std::tuple<K, V>;
  int nworkers = num_workers();

  // sequential elision
  if (n < _hist_seq_threshold || nworkers == 1) {
    size_t pn = pbbslib::log2_up((intT)(n + 1));
    size_t rs = 1L << pn;
    ht.resize(rs);
    sequentialHT<K, V> S(ht.table, n, 1.0f, ht.empty);
    size_t ct = 0;
    for (size_t i = 0; i < n; i++) {
      K k = get_key[i];
      ct += S.insertAdd(k);
    }
    O* out = pbbslib::new_array_no_init<O>(ct);
    size_t k = S.compactInto(apply_f, out);
    return std::make_pair(k, out);
  }

  if (n < 5000000) {
    return histogram_medium<O>(get_key, n, apply_f, ht);
  }

  size_t sqrt = (size_t)ceil(pow(n, 0.5));
  size_t num_buckets = (size_t)(n < 20000000) ? (sqrt / 5) : sqrt;

  num_buckets = std::max(1 << pbbslib::log2_up(num_buckets), 1);
  num_buckets = std::min(num_buckets, _hist_max_buckets);
  size_t bits = pbbslib::log2_up(num_buckets);

  auto gb = get_bucket<K, A>(get_key, n, bits);
  size_t num_heavy = gb.k;
  if (num_heavy == 0) {
    return histogram_medium<O>(get_key, n, apply_f, ht);
  }

  bool heavy = (num_heavy > 0);
  size_t num_total_buckets = (heavy) ? 2 * num_buckets : num_buckets;
  size_t num_actual_buckets = num_buckets + num_heavy;

  K* elms;
  size_t* counts;
  size_t num_blocks;
  if (num_total_buckets <= 256) {
    std::tie(elms, counts, num_blocks) = pbbslib::_count_sort<uint8_t, size_t, K>(
        get_key, gb, n, (uintE)num_total_buckets);
  } else {
    std::tie(elms, counts, num_blocks) = pbbslib::_count_sort<uint16_t, size_t, K>(
        get_key, gb, n, (uintE)num_total_buckets);
  }

  size_t block_size = ((n - 1) / num_blocks) + 1;

#define S_STRIDE 64
  size_t* bkt_counts = pbbslib::new_array_no_init<size_t>(num_total_buckets * S_STRIDE);
  par_for(0, num_actual_buckets, 1, [&] (size_t i) {
    bkt_counts[i * S_STRIDE] = 0;
    if (i == (num_total_buckets - 1)) {
      size_t ct = 0;
      for (size_t j = 0; j < num_blocks; j++) {
        size_t start = std::min(j * block_size, n);
        size_t end = std::min(start + block_size, n);
        ct += (end - start) - counts[j * num_total_buckets + i];
      }
      bkt_counts[i * S_STRIDE] = ct;
    } else {
      size_t ct = 0;
      for (size_t j = 0; j < num_blocks; j++) {
        ct += counts[j * num_total_buckets + i + 1] -
              counts[j * num_total_buckets + i];
      }
      bkt_counts[i * S_STRIDE] = ct;
    }
  });

  sequence<size_t> out_offs = sequence<size_t>(num_buckets + 1);
  sequence<size_t> ht_offs = sequence<size_t>(num_buckets + 1);

  using MO = std::optional<O>;
  MO heavy_cts_stk[128];
  MO* heavy_cts;
  if (heavy) {
    if (num_heavy > 128) {
      heavy_cts = pbbslib::new_array_no_init<MO>(num_heavy);
    } else {
      heavy_cts = heavy_cts_stk;
    }
    //      for (size_t i=0; i<num_heavy; i++) {
    //        heavy_cts[i] = std::nullopt;
    //        std::cout << "cnt i = " << i << " = " <<
    //        bkt_counts[(num_buckets+i)*S_STRIDE] << "\n";
    //      }
  }

  // (2) process each bucket, compute the size of each HT and scan (seq)
  ht_offs[0] = 0;
  size_t min_size = std::numeric_limits<size_t>::max();
  size_t max_size = 0;
  for (size_t i = 0; i < num_buckets; i++) {
    size_t size = bkt_counts[i * S_STRIDE];
    size_t ht_size = 0;
    if (size > 0) {
      ht_size = 1 << pbbslib::log2_up((intT)(size + 1));
    }
    ht_offs[i + 1] = ht_offs[i] + ht_size;
    if (size < min_size) {
      min_size = size;
    }
    if (size > max_size) {
      max_size = size;
    }
  }

  //    if (n > 100000) {
  //      std::cout << "n = " << n << " min size = " << min_size << " max size =
  //      " << max_size << " avg size = " << avg_size << "\n"; //" k = " << gb.k
  //      <<
  //      "\n";
  //    }

  ht.resize(ht_offs[num_buckets]);
  KV* table = ht.table;
  auto empty = ht.empty;

  // (3) insert elms into per-bucket hash table (par)
  {
    auto for_inner = [&](size_t i) {
      if (i < num_buckets) {
        size_t ht_start = ht_offs[i];
        size_t ht_size = ht_offs[i + 1] - ht_start;

        size_t k = 0;
        if (ht_size > 0) {
          KV* my_ht = &(table[ht_start]);
          sequentialHT<K, V> S(my_ht, ht_size, empty);

          for (size_t j = 0; j < num_blocks; j++) {
            size_t start = std::min(j * block_size, n);
            size_t end = std::min(start + block_size, n);
            size_t ct = 0;
            size_t off = 0;
            if (i == (num_total_buckets - 1)) {
              off = counts[j * num_total_buckets + i];
              ct = (end - start) - off;
            } else {
              off = counts[j * num_total_buckets + i];
              ct = counts[j * num_total_buckets + i + 1] - off;
            }
            off += start;
            for (size_t l = 0; l < ct; l++) {
              K a = elms[off + l];
              S.insertAdd(a);
            }
          }

          k = S.compactIntoSelf(apply_f);
        }
        out_offs[i] = k;
      } else {
        // heavy bucket
        size_t bkt_id = i - num_buckets;
        K key = 0; // initializing to get rid of -Wmaybe-uninitialized
        bool is_set = false;
        size_t total_ct = 0;
        for (size_t j = 0; j < num_blocks; j++) {
          size_t start = std::min(j * block_size, n);
          size_t end = std::min(start + block_size, n);
          size_t ct = 0;
          size_t off = 0;
          if (i == (num_total_buckets - 1)) {
            off = counts[j * num_total_buckets + i];
            ct = (end - start) - off;
          } else {
            off = counts[j * num_total_buckets + i];
            ct = counts[j * num_total_buckets + i + 1] - off;
          }
          off += start;
          if (!is_set && ct) {
            key = elms[off];
            is_set = true;
          }
          total_ct += ct;
        }
        assert(is_set);

        std::optional<O> value = apply_f(std::make_tuple(key, total_ct));
        heavy_cts[bkt_id] = value;
      }
    };
    par_for(0, num_actual_buckets, 1, [&] (size_t i)
                    { for_inner(i); });
  }

  // (4) scan
  size_t ct = 0;
  for (size_t i = 0; i < num_buckets; i++) {
    size_t s = ct;
    ct += out_offs[i];
    out_offs[i] = s;
  }

  out_offs[num_buckets] = ct;
  size_t heavy_start = ct;

  if (heavy) {
    for (size_t i = 0; i < num_heavy; i++) {
      if (heavy_cts[i].has_value()) {
        ct++;
      }
    }
  }
  uintT num_distinct = ct;

  O* res = pbbslib::new_array_no_init<O>(ct);

  // (5) map compacted hts to output, clear hts
  par_for(0, num_buckets, 1, [&] (size_t i) {
    size_t o = out_offs[i];
    size_t k = out_offs[(i + 1)] - o;

    size_t ht_start = ht_offs[i];
    size_t ht_size = ht_offs[i + 1] - ht_start;

    if (ht_size > 0) {
      KV* my_ht = &(table[ht_start]);

      for (size_t j = 0; j < k; j++) {
        res[o + j] = my_ht[j];
        my_ht[j] = empty;
      }
    }
  });

  if (heavy) {
    size_t heavy_off = 0;
    for (size_t i = 0; i < num_heavy; i++) {
      auto& value = heavy_cts[i];
      if (value.has_value()) {
        res[heavy_start + heavy_off++] = *value;
      }
    }
  }

  pbbslib::free_array(elms);
  pbbslib::free_array(counts);
  pbbslib::free_array(bkt_counts);
  if (heavy && num_heavy > 128) {
    pbbslib::free_array(heavy_cts);
  }

  return std::make_pair(num_distinct, res);
}

template <class E, class O, class K, class V, class A, class Reduce,
          class Apply>
inline std::pair<size_t, O*> seq_histogram_reduce(A& get_elm, size_t n,
                                                  Reduce& reduce_f,
                                                  Apply& apply_f,
                                                  hist_table<K, V>& ht) {
  size_t pn = pbbslib::log2_up((intT)(n + 1));
  size_t rs = 1L << pn;
  ht.resize(rs);
  sequentialHT<K, V> S(ht.table, n, 1.0f, ht.empty);
  for (size_t i = 0; i < n; i++) {
    E a = get_elm[i];
    reduce_f(S, a);
  }
  O* out = pbbslib::new_array_no_init<O>(n);
  size_t k = S.compactInto(apply_f, out);
  return std::make_pair(k, out);
}

// Issue: want to make the type that's count-sort'd independent of K,V
// A : <E> elms
// B : inmap<keys> (numeric so that we can hash)
// E : type of intermediate elements (what we count sort)
// Q : reduction (seqHT<uintE, E>&, E&) -> void
// F : std::tuple<K, Elm> -> std::optional<uintE>
template <class E, class O, class K, class V, class A, class B, class Reduce,
          class Apply>
inline std::pair<size_t, O*> histogram_reduce(A& get_elm, B& get_key, size_t n,
                                              Reduce& reduce_f, Apply& apply_f,
                                              hist_table<K, V>& ht) {
  typedef std::tuple<K, V> KV;

  int nworkers = num_workers();

  if (n < _hist_seq_threshold || nworkers == 1) {
    auto r = seq_histogram_reduce<E, O>(get_elm, n, reduce_f, apply_f, ht);
    return r;
  }

  size_t sqrt = (size_t)ceil(pow(n, 0.5));
  size_t num_buckets = (size_t)(n < 20000000) ? (sqrt / 5) : sqrt;

  num_buckets = std::max(1 << pbbslib::log2_up(num_buckets), 1);
  num_buckets = std::min(num_buckets, _hist_max_buckets);

  // (1) count-sort based on bucket
  size_t low_mask = ~((size_t)15);
  size_t bucket_mask = num_buckets - 1;
  auto gb = [&](uintE i) {
    return pbbslib::hash32(get_key[i] & low_mask) & bucket_mask;
  };

  auto p = pbbslib::_count_sort<int16_t, size_t, E>(get_elm, gb, n, (uintE)num_buckets);

  auto elms = std::get<0>(p);  // count-sort'd
  // laid out as num_buckets (row), blocks (col)
  size_t* counts = std::get<1>(p);
  size_t num_blocks = std::get<2>(p);
  size_t block_size = ((n - 1) / num_blocks) + 1;

#define S_STRIDE 64
  size_t* bkt_counts = pbbslib::new_array_no_init<size_t>(num_buckets * S_STRIDE);
  par_for(0, num_buckets, 1, [&] (size_t i) {
    bkt_counts[i * S_STRIDE] = 0;
    if (i == (num_buckets - 1)) {
      size_t ct = 0;
      for (size_t j = 0; j < num_blocks; j++) {
        size_t start = std::min(j * block_size, n);
        size_t end = std::min(start + block_size, n);
        ct += (end - start) - counts[j * num_buckets + i];
      }
      bkt_counts[i * S_STRIDE] = ct;
    } else {
      size_t ct = 0;
      for (size_t j = 0; j < num_blocks; j++) {
        ct += counts[j * num_buckets + i + 1] - counts[j * num_buckets + i];
      }
      bkt_counts[i * S_STRIDE] = ct;
    }
  });

  sequence<size_t> out_offs = sequence<size_t>(num_buckets + 1);
  sequence<size_t> ht_offs = sequence<size_t>(num_buckets + 1);

  // (2) process each bucket, compute the size of each HT and scan (seq)
  ht_offs[0] = 0;
  for (size_t i = 0; i < num_buckets; i++) {
    size_t size = bkt_counts[i * S_STRIDE];
    size_t ht_size = 0;
    if (size > 0) {
      ht_size = 1 << pbbslib::log2_up((intT)(size + 1));
    }
    ht_offs[i + 1] = ht_offs[i] + ht_size;
  }

  ht.resize(ht_offs[num_buckets]);
  O* out = pbbslib::new_array_no_init<O>(ht_offs[num_buckets]);

  // (3) insert elms into per-bucket hash table (par)
  {
    auto for_inner = [&](size_t i) {
      size_t ht_start = ht_offs[i];
      size_t ht_size = ht_offs[i + 1] - ht_start;

      size_t k = 0;
      KV* table = ht.table;
      if (ht_size > 0) {
        KV* my_ht = &(table[ht_start]);
        sequentialHT<K, V> S(my_ht, ht_size, ht.empty);

        O* my_out = &(out[ht_start]);

        for (size_t j = 0; j < num_blocks; j++) {
          size_t start = std::min(j * block_size, n);
          size_t end = std::min(start + block_size, n);
          size_t ct = 0;
          size_t off = 0;
          if (i == (num_buckets - 1)) {
            off = counts[j * num_buckets + i];
            ct = (end - start) - off;
          } else {
            off = counts[j * num_buckets + i];
            ct = counts[j * num_buckets + i + 1] - off;
          }
          off += start;
          for (size_t l = 0; l < ct; l++) {
            E a = elms[off + l];
            reduce_f(S, a);
          }
        }

        k = S.compactInto(apply_f, my_out);
      }
      out_offs[i] = k;
    };
    par_for(0, num_buckets, 1, [&] (size_t i) { for_inner(i); });
  }

  // (4) scan
  size_t ct = 0;
  for (size_t i = 0; i < num_buckets; i++) {
    size_t s = ct;
    ct += out_offs[i];
    out_offs[i] = s;
  }
  out_offs[num_buckets] = ct;
  uintT num_distinct = ct;

  O* res = pbbslib::new_array_no_init<O>(ct);

  // (5) map compacted hts to output, clear hts
  par_for(0, num_buckets, 1, [&] (size_t i) {
    size_t o = out_offs[i];
    size_t k = out_offs[(i + 1)] - o;

    size_t ht_start = ht_offs[i];
    size_t ht_size = ht_offs[i + 1] - ht_start;

    if (ht_size > 0) {
      O* my_out = &(out[ht_start]);

      for (size_t j = 0; j < k; j++) {
        res[o + j] = my_out[j];
      }
    }
  });

  pbbslib::free_array(elms);
  pbbslib::free_array(counts);
  pbbslib::free_array(bkt_counts);
  pbbslib::free_array(out);
  return std::make_pair(num_distinct, res);
}

}  // namespace gbbs
