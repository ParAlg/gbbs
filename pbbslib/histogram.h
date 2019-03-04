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
#include <math.h>
#include <stdio.h>
#include <cstdint>
#include <vector>
#include <algorithm>
#include "utilities.h"
#include "counting_sort.h"

namespace pbbs {

  template <typename s_size_t, typename Seq>
  sequence<s_size_t> seq_histogram(Seq const &A, size_t m) {
    sequence<s_size_t> counts(m);
    for (size_t i = 0; i < m; i++)
      counts[i] = 0;

    for (size_t i = 0; i < A.size(); i++)
      counts[A[i]]++;
    return counts;
  }

  template <typename Seq, typename Iter>
  void _seq_count(Seq const &In, range<Iter> counts) {
    for (size_t i = 0; i < counts.size(); i++) counts[i] = 0;
    for (size_t j = 0; j < In.size(); j++) counts[In[j]]++;
  }

  template <typename s_size_t, typename Seq>
  sequence<s_size_t> _count(Seq const &In, size_t num_buckets) {
    sequence<s_size_t> counts(num_buckets);
    size_t n = In.size();

    if (n < ((size_t) 1 << 14)) {
      _seq_count(In, counts.slice());
      return counts;
    }

    size_t num_threads = num_workers();
    size_t num_blocks = std::min((size_t) (1 + n/(num_buckets*32)),
				 num_threads*4);
    size_t block_size = ((n-1)/num_blocks) + 1;
    size_t m = num_blocks * num_buckets;

    sequence<s_size_t> block_counts(m);

    // count each block
    auto block_f = [&] (size_t i) {
      size_t start = std::min(i * block_size, n);
      size_t end = std::min((i+1) * block_size, n);
      auto bc = block_counts.slice(i*num_buckets,(i+1)*num_buckets);
      _seq_count(In.slice(start,end), bc);
    };
    parallel_for(0, num_blocks, block_f, 1);

    auto bucket_f = [&] (size_t j) {
      size_t sum = 0;
      for (size_t i = 0; i < num_blocks; i++)
	sum += block_counts[i*num_buckets+j];
      counts[j] = sum;
    };
    if (m >= (1 << 14))
      parallel_for(0, num_buckets, bucket_f, 1);
    else
      for (size_t j = 0; j < num_buckets; j++) bucket_f(j);
    return counts;
  }

  template <typename Seq>
  struct get_bucket {
    using E = typename Seq::value_type;
    using HE = std::pair<E,int>;
    std::vector<HE> hash_table;
    size_t table_mask;
    size_t low_mask;
    size_t bucket_mask;
    int num_buckets;
    int k;
    const Seq& I;

    get_bucket() {}

    std::pair<std::vector<E>,int> heavy_hitters(const Seq &A, size_t count) {
      size_t n = A.size();
      std::vector<E> sample(count);
      for (size_t i = 0; i < count; i++)
	sample[i] = A[hash64(i)%n];
      std::sort(sample.begin(), sample.end());

      // only keep those with at least two copies
      int k = 0;
      int c = 0;
      for (size_t i = 1; i < count; i++)
	if (sample[i] == sample[i-1]) {
	  if (c++ == 0) sample[k++] = sample[i];
	} else c = 0;
      return std::make_pair(sample,k);
    }

    void make_hash_table(std::vector<E> const &entries, size_t n,
			 size_t table_size, size_t table_mask) {
      hash_table = std::vector<HE>(table_size);
      for (size_t i=0; i < table_size; i++)
	hash_table[i] = std::make_pair(0,-1);

      for (size_t i = 0; i < n; i++)
	hash_table[hash64(entries[i])&table_mask] = std::make_pair(entries[i],i);
    }

    get_bucket(Seq const &A, size_t bits) : I(A) {
      num_buckets = 1 << bits;
      bucket_mask = num_buckets-1;
      low_mask = ~((size_t) 15);
      int count = 2 * num_buckets;
      int table_size = 4 * count;
      table_mask = table_size-1;

      std::vector<E> sample;
      std::tie(sample,k) = heavy_hitters(A, count);

      if (k > 0)
	make_hash_table(sample, k, table_size, table_mask);
    }

    size_t operator() (size_t i) const {
      if (k > 0) {
        std::pair<E,int> h = hash_table[hash64(I[i])&table_mask];
	if (h.first == I[i] && h.second != -1)
	  return h.second + num_buckets;
      }
      return pbbs::hash64_2(I[i] & low_mask) & bucket_mask;
    }

  };

  
  template <typename s_size_t, typename Seq>
  sequence<s_size_t> histogram(Seq const &A, size_t m) {
    size_t n = A.size();
    size_t bits;
    using T = typename Seq::value_type;
    
    if (n < (1 << 27)) bits = (log2_up(n) - 7)/2;
    // for large n selected so each bucket fits into cache
    else bits = (log2_up(n) - 17);
    size_t num_buckets = (1<<bits);
    if (m < n / num_buckets)
      return  _count<s_size_t>(A, m);
    if (n < (1 << 13))
      return seq_histogram<s_size_t>(A , m);

    //timer t("histogram", false);

    // Returns a map (hash) from key to bucket.
    // Keys with many elements (big) have their own bucket while
    // others share a bucket.
    // Keys that share low 4 bits get same bucket unless big.
    // This is to avoid false sharing.
    get_bucket<Seq> x(A, bits-1);
    auto get_buckets = delayed_seq<size_t>(n, x);
    sequence<T> B(n);

    // first buckets based on hash using a counting sort
    sequence<size_t> bucket_offsets
      = count_sort(A, B.slice(), get_buckets, num_buckets);
    //t.next("send to buckets");
    
    // note that this is cache line alligned
    sequence<s_size_t> counts(m, 0);
    
    // now sequentially process each bucket
    auto bucket_f = [&] (size_t i) {
      size_t start = bucket_offsets[i];
      size_t end = bucket_offsets[i+1];

      // small buckets have indices in bottom half
      if (i < num_buckets/2)
	for (size_t i = start; i < end; i++)
	  counts[B[i]]++;

      // large buckets have indices in top half
      else if (end > start)
	counts[B[i]] = end-start;
    };
    parallel_for(0, num_buckets, bucket_f, 1);
    //t.next("within buckets ");
    return counts;
  }
}
