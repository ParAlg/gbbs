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
#include <stdio.h>
#include <math.h>
#include "utilities.h"
#include "sequence_ops.h"
#include "transpose.h"
#include "histogram.h"

// Supports functions that take a seq of key-value pairs, collects all the
// keys with the same value, and sums them up with a binary function (monoid)
//
// For the first one the keys must be integers in the range [0,num_buckets).
// It returns a sequence of length num_buckets, one sum per possible key value.
//
//   template <typename Seq, typename M>
//   sequence<typename Seq::value_type::second_type>
//   collect_reduce(Seq const &A, M const &monoid, size_t num_buckets);
//
// For the second one keys can be any integer values
// It returns a sequence of key-value pairs.  If a key appeared at least
//   once, an entry with the sum for that key will appear in the output.
// The output is not necessarily sorted by key value
//
//   template <typename Seq, typename M>
//   sequence<typename Seq::value_type>
//   collect_reduce_sparse(Seq const &A, M const &monoid);

namespace pbbs {

  // the following parameters can be tuned
  constexpr const size_t CR_SEQ_THRESHOLD = 8192;

  template <class Seq, class OutSeq, class M>
  void seq_collect_reduce_few(Seq const &A, OutSeq &&Out, M const &monoid, size_t num_buckets) {
    size_t n = A.size();
    for (size_t i = 0; i < num_buckets; i++) Out[i] = monoid.identity;
    for (size_t j = 0; j < n; j++) {
      size_t k = A[j].first;
      Out[k] = monoid.f(Out[k],A[j].second);
    }
  }

  template <class Seq, class M>
  sequence <typename Seq::value_type::second_type>
  seq_collect_reduce_few(Seq const &A, M const &monoid, size_t num_buckets) {
      using val_type = typename Seq::value_type::second_type;
      sequence<val_type> Out(num_buckets);
      seq_collect_reduce_few(A, Out, monoid, num_buckets);
      return Out;
    }

  // This one is for few buckets (e.g. less than 2^16)
  //  A is a sequence of key-value pairs
  //  monoid has fields m.identity and m.f (a binary associative function)
  //  all keys must be smaller than num_buckets
  template <typename Seq, typename M>
  sequence<typename Seq::value_type::second_type>
  collect_reduce_few(Seq const &A, M const &monoid, size_t num_buckets) {
    using val_type = typename Seq::value_type::second_type;
    size_t n = A.size();
    timer t;
    t.start();
    
    // pad to 16 buckets to avoid false sharing (does not affect results)
    num_buckets = std::max(num_buckets, (size_t) 16);

    //size_t num_blocks = ceil(pow(n/num_buckets,0.5));
    size_t num_threads = num_workers();
    size_t num_blocks = std::min(4*num_threads, n/num_buckets/64);

    num_blocks = 1 << log2_up(num_blocks);

    sequence<val_type> Out(num_buckets);

    // if insufficient parallelism, do sequentially
    if (n < CR_SEQ_THRESHOLD || num_blocks == 1 || num_threads == 1)
      return seq_collect_reduce_few(A, monoid, num_buckets);

    size_t block_size = ((n-1)/num_blocks) + 1;
    size_t m = num_blocks * num_buckets;

    sequence<val_type> OutM(m);

    sliced_for(n, block_size, [&] (size_t i, size_t start, size_t end) {
	seq_collect_reduce_few(A.slice(start,end),
			       OutM.slice(i*num_buckets,(i+1)*num_buckets),
			       monoid, num_buckets);
      });

    parallel_for (0, num_blocks, [&] (size_t i) {
	val_type o_val = monoid.identity;
	for (size_t j = 0; j < num_blocks; j++)
	  o_val = monoid.f(o_val, OutM[i + j*num_buckets]);
	Out[i] = o_val;
      }, 1);

    return Out;
  }

  template <typename E>
  struct pair_hasheq_mask_low {
    static inline size_t hash(E a) {return hash64_2(a.first & ~((size_t) 15));}
    static inline bool eql(E a, E b) {return a.first == b.first;}
  };

  template <typename Seq, typename M>
  sequence<typename Seq::value_type::second_type>
  collect_reduce(Seq const &A, M const &monoid, size_t num_buckets) {
    using T = typename Seq::value_type;
    using val_type = typename T::second_type;
    size_t n = A.size();
    
    // #bits is selected so each block fits into L3 cache
    //   assuming an L3 cache of size 1M per thread
    // the counting sort uses 2 x input size due to copy
    size_t cache_per_thread = 1000000;
    size_t bits = log2_up(2 * (size_t) sizeof(val_type) * n / cache_per_thread);
    size_t num_blocks = (1<<bits);
    
    if (num_buckets <= 4 * num_blocks) 
      return collect_reduce_few(A, monoid, num_buckets);

    // Returns a map (hash) from key to block.
    // Keys with many elements (big) have their own block while
    // others share a block.
    // Keys that share low 4 bits get same block unless big.
    // This is to avoid false sharing.
    //auto get_i = [&] (size_t i) -> size_t {return A[i].first;};
    //auto s = delayed_seq<size_t>(n,get_i);

    using hasheq = pair_hasheq_mask_low<T>;
    get_bucket<T,hasheq> gb(A, hasheq(), bits);
    //get_bucket<decltype(s)> x(s, bits-1);
    //auto get_blocks = delayed_seq<size_t>(n, x);
    sequence<T> B = sequence<T>::no_init(n);
    sequence<T> Tmp = sequence<T>::no_init(n);

    // first partition into blocks based on hash using a counting sort
    sequence<size_t> block_offsets;
    //bool single_block;
    //std::tie(block_offsets, single_block)
      //= count_sort(A, B.slice(), Tmp.slice(), get_blocks, num_blocks);
    block_offsets = integer_sort_(A, B.slice(), Tmp.slice(), gb, 
				  bits, num_blocks, false);
    // note that this is cache line alligned
    sequence<val_type> sums(num_buckets, monoid.identity);
       
    // now process each block in parallel
    parallel_for(0, num_blocks, [&] (size_t i) {
	size_t start = block_offsets[i];
	size_t end = block_offsets[i+1];
	size_t cut =  gb.heavy_hitters ? num_buckets/2 : num_buckets;
	
	// small blocks have indices in bottom half
	if (i < cut)
	  for (size_t i = start; i < end; i++) {
	    size_t j = B[i].first;
	    sums[j] = monoid.f(sums[j], B[i].second);
	  }

	// large blocks have indices in top half
	else if (end > start) {
	  auto x = [&] (size_t i) {return B[i].second;};
	  auto vals = delayed_seq<size_t>(n, x);
	  sums[B[i].first] = reduce(vals, monoid);
	}
      }, 1);
    return sums;
  }

  // this one is for more buckets than the length of A (i.e. sparse)
  //  A is a sequence of key-value pairs
  //  monoid has fields m.identity and m.f (a binary associative function)
  template <typename Seq, typename HashEq, typename M>
  sequence<typename Seq::value_type>
  collect_reduce_sparse(Seq const &A, HashEq hasheq, M const &monoid) {
    using T = typename Seq::value_type;
    using val_type = typename T::second_type;
    
    timer t("collect_reduce_sparse", false);
    size_t n = A.size();

    if (n < 1000) {
      auto cmp = [] (T a, T b) {return a.first < b.first;};
      sequence<T> B = sample_sort(A, cmp);
      size_t j = 0;
      for (size_t i = 1; i < n; i++) {
	if (B[i].first == B[j].first)
	  B[j].second = monoid.f(B[j].second, B[i].second);
	else B[j++] = B[i];
      };
      return sequence<T>(j, [&] (size_t i) {return B[i];});
    }

    // #bits is selected so each block fits into L3 cache
    //   assuming an L3 cache of size 1M per thread
    // the counting sort uses 2 x input size due to copy
    size_t cache_per_thread = 1000000;
    size_t bits = log2_up((size_t) ((1.2 * 2 * sizeof(T) * n) / (float) cache_per_thread));
    size_t num_buckets = (1<<bits);

    // Returns a map (hash) from key to bucket.
    // Keys with many elements (big) have their own bucket while
    // others share a bucket.
    // Keys that share low 4 bits get same bucket unless big.
    // This is to avoid false sharing.
    sequence<T> B = sequence<T>::no_init(n);
    sequence<T> Tmp = sequence<T>::no_init(n);
    
    // first buckets based on hash using a counting sort
    get_bucket<T,HashEq> gb(A, hasheq, bits);
    sequence<size_t> bucket_offsets =
      integer_sort_(A.slice(), B.slice(), Tmp.slice(), gb,
		    bits, num_buckets, false);
    t.next("sort to blocks");
    
    // note that this is cache line alligned
    size_t num_tables = gb.heavy_hitters ? num_buckets/2 : num_buckets;
    size_t bucket_size = (n - 1) / num_tables + 1;
    float factor = 1.2;
    if (bucket_size < 128000)
      factor += (17 - log2_up(bucket_size))*.15;
    size_t table_size = (factor * bucket_size);
    size_t total_table_size = table_size * num_tables;
    sequence<T> table = sequence<T>::no_init(total_table_size);
    sequence<size_t> sizes(num_tables + 1);
	
    // now in parallel process each bucket sequentially
    parallel_for(0, num_tables, [&] (size_t i) {
      T* my_table = table.begin() + i * table_size;

      sequence<bool> flags(table_size, false);
      // clear tables
      //for (size_t i = 0; i < table_size; i++)
      // assign_uninitialized(my_table[i], T(empty, identity));

      // insert small bucket (ones with multiple different items)
      size_t start = bucket_offsets[i];
      size_t end = bucket_offsets[i+1];
      if ((end-start) > table_size)
	cout << "error in collect_reduce: " << (end-start) << ", " << table_size
	     << ", " << total_table_size << ", " << n << endl;
      for (size_t j = start; j < end; j++) {
	size_t idx = B[j].first;
	size_t k = ((uint) hasheq.hash(B[j])) % table_size;
	while (flags[k] && my_table[k].first != idx)
	  k = (k + 1 == table_size) ? 0 : k + 1;
	if (flags[k])
	  my_table[k] = T(idx, monoid.f(my_table[k].second, B[j].second));
	else {
	  flags[k] = true;
	  assign_uninitialized(my_table[k], T(idx, B[j].second));
	}
      }

      // now if there are any "heavy hitters" (buckets with a single item)
      // insert them
      if (gb.heavy_hitters) {
	size_t start_l = bucket_offsets[num_tables + i];
	size_t len = bucket_offsets[num_tables + i + 1] - start_l;
	if (len > 0) {
	  auto f = [&] (size_t i) -> val_type {return B[i+start_l].second;};
	  auto s = delayed_seq<val_type>(len, f);
	  val_type x = reduce(s, monoid);
	  size_t j = 0;
	  while (flags[j])
	    j = (j + 1 == table_size) ? 0 : j + 1;
	  assign_uninitialized(my_table[j], T(B[start_l].first, x));
	}
      }

      // pack tables down to bottom
      size_t j=0;
      for (size_t i = 0; i < table_size; i++)
	if (flags[i]) 
	  move_uninitialized(my_table[j++], my_table[i]);
      sizes[i] = j;

    }, 0);
    t.next("hash into tables");

    sizes[num_tables] = 0;
    size_t total = scan_inplace(sizes.slice(), addm<size_t>());

    // copy packed tables into contiguous result
    sequence<T> result = sequence<T>::no_init(total);
    auto copy_f = [&] (size_t i) {
      size_t d_offset = sizes[i];
      size_t s_offset = i * table_size;
      size_t len = sizes[i+1] - sizes[i];
      for (size_t j = 0; j < len; j++)
	move_uninitialized(result[d_offset+j], table[s_offset+j]);
    };
    parallel_for(0, num_tables, copy_f, 1);
    t.next("copy subresults");
    return result;
  }
}
