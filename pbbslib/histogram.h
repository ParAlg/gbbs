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
    parallel_for(0, num_blocks, [&] (size_t i) {
      size_t start = std::min(i * block_size, n);
      size_t end = std::min((i+1) * block_size, n);
      auto bc = block_counts.slice(i*num_buckets,(i+1)*num_buckets);
      _seq_count(In.slice(start,end), bc);
    }, 1);


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

  // The idea is to return a hash function that maps any items
  // that appear many times into their own bucket.
  // Otherwise items can end up in the same bucket.
  // E is the type of element
  // HashEq must contain an hash function E -> size_t
  //    and an equality function E x E -> bool
  template <typename E, typename HashEq>
  struct get_bucket_old {
    using HE = std::pair<E,int>;
    sequence<HE> hash_table;
    size_t table_mask;
    size_t bucket_mask;
    size_t num_buckets;
    bool heavy_hitters;
    const HashEq heq;

    // creates a structure from a sequence of elements
    // bits is the number of bits that will be returned by the hash function
    // items that appear many times will be mapped individually into
    //    the top half [2^{bits-1},2^{bits})
    //    and light items shared into the bottom half [0,2^{bits-1})
    template <typename Seq>
    get_bucket_old(Seq const &A, HashEq const &heq, size_t bits) : heq(heq) {
      size_t n = A.size();
      size_t low_bits = bits - 1;  // for the bottom half
      num_buckets = 1 << low_bits; // in bottom half
      size_t count = 2 * num_buckets;
      size_t table_size = 4 * count;
      table_mask = table_size-1;
      
      hash_table = sequence<HE>(table_size, std::make_pair(E(),-1));
      
      // insert sample into hash table with one less than the
      // count of how many times appears (since it starts with -1)
      for (size_t i = 0; i < count; i++) {
	E s = A[hash64(i)%n];
	size_t idx = heq.hash(s) & table_mask;
	while (1) {
	  if (hash_table[idx].second == -1) {
	    hash_table[idx] = std::make_pair(s,0);
	    break;}
	  else if (heq.eql(hash_table[idx].first, s)) {
	      hash_table[idx].second += 1;
	      break;
	    }
	  else idx = (idx + 1) & table_mask;
	}
      }

      // keep in the hash table if at least three copies and give kept items
      // consecutive numbers.   k will be total kept items.
      size_t k = 0;
      for (size_t i = 0; i < table_size; i++) {
	if (hash_table[i].second > 1) {
	  E key = hash_table[i].first;
	  size_t idx = heq.hash(key) & table_mask;
	  hash_table[idx] = std::make_pair(key, k++);
	}
	else hash_table[i].second = -1;
      }

      heavy_hitters = (k > 0);
      bucket_mask = heavy_hitters ? num_buckets-1 : 2*num_buckets-1;
    }

    // the hash function.
    // uses chosen id if key appears many times (top half)
    // otherwise uses (heq.hash(v) % num_buckets) directly (bottom half)
    size_t operator() (E v) const {
      if (heavy_hitters) {
        auto h = hash_table[heq.hash(v) & table_mask];
	if (h.second != -1 && heq.eql(h.first, v))
	  return h.second + num_buckets; // top half
      }
      return heq.hash(v) & bucket_mask; // bottom half
    }

  };

  // mask low so elements on same cache line end up
  // in same bucket
  template <typename intType>
  struct int_hasheq_mask_low {
    static size_t hash(intType a) {return hash64_2(a & ~((size_t) 31));}
    static bool eql(intType a, intType b) {return a == b;}
  };

  template <typename s_size_t, typename Seq>
  sequence<s_size_t> histogram(Seq const &A, size_t m) {
    size_t n = A.size();
    using T = typename Seq::value_type;

    // #bits is selected so each block fits into L3 cache
    //   assuming an L3 cache of size 1M per thread
    // the counting sort uses 2 x input size due to copy
    size_t cache_per_thread = 1000000;
    size_t bits = std::max<size_t>(log2_up(1 + 2 * (size_t) sizeof(T) * n / cache_per_thread),
				   4);
    //if (bits == 0) 
    //  return seq_histogram<s_size_t>(A, m);
    size_t num_buckets = (1<<bits);
    if (m < n / num_buckets)
      return  _count<s_size_t>(A, m);
    if (n < (1 << 13))
      return seq_histogram<s_size_t>(A , m);

    timer t("histogram", false);
    
    sequence<T> B(n);
    sequence<T> Tmp(n);

    // gb is a map (hash) from key to bucket.
    // Keys with many elements (big) have their own bucket while
    // others share a bucket.
    // Keys that share low 4 bits get same bucket unless big.
    // This is to avoid false sharing.
    get_bucket_old<T,int_hasheq_mask_low<T>> gb(A, int_hasheq_mask_low<T>(), bits);
    t.next("head");

    // first buckets based on hash using a counting sort
    sequence<size_t> bucket_offsets =
      integer_sort_(A.slice(), B.slice(), Tmp.slice(), gb,
		    bits, num_buckets, false);
    t.next("send to buckets");

    // note that this is cache line alligned
    sequence<s_size_t> counts(m, (s_size_t) 0);
    t.next("initialize buckets");

    // now in parallel across the buckets, sequentially process each bucket
    parallel_for(0, num_buckets, [&] (size_t i) {
      size_t start = bucket_offsets[i];
      size_t end = bucket_offsets[i+1];
      size_t cut =  gb.heavy_hitters ? num_buckets/2 : num_buckets;
      // small buckets have indices in bottom half
      if (i < cut)
	for (size_t j = start; j < end; j++)
	  counts[B[j]]++;

      // large buckets have indices in top half
      else if (end > start) {
	counts[B[start]] = end-start;
      }
    }, 1);
    t.next("within buckets ");
    return counts;
  }
}
