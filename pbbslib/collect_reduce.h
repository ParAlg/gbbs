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

// TODO
//  make few buckets special case of many buckets
//  collect_reduce mutates inputs

namespace pbbs {

  // the following parameters can be tuned
  constexpr const size_t CR_SEQ_THRESHOLD = 8192;

  template<class OutVal, class InSeq, class KeySeq, class M>
  void _seq_collect_reduce(InSeq In, KeySeq Keys, OutVal* Out, size_t num_buckets, M monoid) {
    size_t n = In.size();
    for (size_t i = 0; i < num_buckets; i++) Out[i] = monoid.identity;
    for (size_t j = 0; j < n; j++) {
      size_t k = Keys[j];
      Out[k] = monoid.f(Out[k],In[j]);
    }
  }

  template<class OutVal, class InSeq, class KeySeq, class M>
  sequence<OutVal> seq_collect_reduce(InSeq In, KeySeq Keys,
				      size_t num_buckets, M monoid) {
    OutVal* Out = new_array<OutVal>(num_buckets);
    _seq_collect_reduce(In, Keys, Out, num_buckets, monoid);
    return sequence<OutVal>(Out,num_buckets);
  }

  // this one is for few buckets (e.g. less than 2^16)
  //  In is the sequence of values to be reduced
  //  Keys is the sequence of bucket numbers
  //  num_buckets is the number of buckets (all keys need to be less)
  //  monoid has fields m.identity and m.f (a binary associative function)
  template<class OutVal, class InSeq, class KeySeq, class M>
  sequence<OutVal> collect_reduce(InSeq const &In, KeySeq const &Keys,
				  size_t num_buckets, M monoid) {
    size_t n = In.size();
    timer t;
    t.start();
    
    // pad to 16 buckets to avoid false sharing (does not affect results)
    num_buckets = std::max(num_buckets, (size_t) 16);

    //size_t num_blocks = ceil(pow(n/num_buckets,0.5));
    size_t num_threads = num_workers();
    size_t num_blocks = std::min(4*num_threads, n/num_buckets/64);

    num_blocks = 1 << log2_up(num_blocks);

    sequence<OutVal> Out(num_buckets);

    // if insufficient parallelism, sort sequentially
    if (n < CR_SEQ_THRESHOLD || num_blocks == 1 || num_threads == 1)
      return seq_collect_reduce<OutVal>(In,Keys,num_buckets, monoid);

    size_t block_size = ((n-1)/num_blocks) + 1;
    size_t m = num_blocks * num_buckets;

    OutVal *OutM = new_array<OutVal>(m);
    
    auto block_f = [&] (size_t i) {
      size_t start = std::min(i * block_size, n);
      size_t end =  std::min(start + block_size, n);
      _seq_collect_reduce<OutVal>(In.slice(start,end), Keys.slice(start,end),
				  OutM + i*num_buckets, num_buckets, monoid);
    };
    parallel_for (0, num_blocks, block_f, 1);
    
    auto sum_buckets = [&] (size_t i) {
      OutVal o_val = monoid.identity;
      for (size_t j = 0; j < num_blocks; j++)
	o_val = monoid.f(o_val, OutM[i + j*num_buckets]);
      Out[i] = o_val;
    };
    parallel_for(0, num_buckets, sum_buckets, 1);
    delete_array(OutM, m);
    return Out;
  }

  // this one is for many buckets but no more than len(A) buckets
  //      (e.g. more than 2^16)
  //  A is the input sequence
  //  m is the number of buckets
  //  get_index is a function that gets the bucket from an element of A
  //  get_val is a function that gets the value from an element of A
  //  monoid has fields m.identity and m.f (a binary associative function)
  template <typename OT, typename Seq, typename F, typename G, typename M>
  sequence<OT> collect_reduce(Seq const &A, size_t m,
			      F get_index, G get_val, M monoid) {
    using T = typename Seq::value_type;
    size_t n = A.size();
    size_t bits;

    if (n < (1 << 27)) bits = (log2_up(n) - 7)/2;
    // for large n selected so each bucket fits into cache
    else bits = (log2_up(n) - 17);
    size_t num_buckets = (1<<bits);

    // Returns a map (hash) from key to bucket.
    // Keys with many elements (big) have their own bucket while
    // others share a bucket.
    // Keys that share low 4 bits get same bucket unless big.
    // This is to avoid false sharing.
    auto get_i = [&] (size_t i) -> size_t {return get_index(A[i]);};
    auto s = delayed_seq<size_t>(n,get_i);
    get_bucket<decltype(s)> x(s, bits-1);
    auto get_buckets = delayed_seq<size_t>(n, x);
    //sequence<size_t> get_buckets(n, [&] (size_t i) {return 0;});
    sequence<T> B(n);

    // first buckets based on hash using a counting sort
    sequence<size_t> bucket_offsets 
      = count_sort(A, B.slice(), get_buckets, num_buckets);

    // note that this is cache line alligned
    sequence<OT> sums(m, monoid.identity);
       
    // now process each bucket in parallel
    parallel_for(0, num_buckets, [&] (size_t i) {
	size_t start = bucket_offsets[i];
	size_t end = bucket_offsets[i+1];

	// small buckets have indices in bottom half
	if (i < num_buckets/2)
	  for (size_t i = start; i < end; i++) {
	    size_t j = get_index(B[i]);
	    sums[j] = monoid.f(sums[j], get_val(B[i]));
	  }

	// large buckets have indices in top half
	else if (end > start) {
	  auto x = [&] (size_t i) {return get_val(B[i]);};
	  auto vals = delayed_seq<size_t>(n, x);
	  sums[get_index(B[i])] = reduce(vals, monoid);
	}
      }, 1);
    return sums;
  }

  // this one is for more buckets than the length of A (i.e. sparse)
  //  A is a sequence of key-value pairs
  //  monoid has fields m.identity and m.f (a binary associative function)
  template <typename Seq, typename M>
  sequence<typename Seq::value_type>
  collect_reduce_pair(Seq const &A, M const &monoid) {
    using T = typename Seq::value_type;
    using key_type = typename T::first_type;
    using val_type = typename T::second_type;
    key_type empty = (key_type) -1;
    val_type identity = monoid.identity;
    
    timer t;
    t.start();
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
      
    size_t bits;
    if (n < (1 << 19)) bits = (log2_up(n))/2;
    else if (n < (1 << 23)) bits = (log2_up(n) - 3)/2;
    else if (n < (1 << 27)) bits = (log2_up(n) - 5)/2;
    // for large n selected so each bucket fits into cache
    else bits = (log2_up(n) - 17);
    size_t num_buckets = (1<<bits);

    // Returns a map (hash) from key to bucket.
    // Keys with many elements (big) have their own bucket while
    // others share a bucket.
    // Keys that share low 4 bits get same bucket unless big.
    // This is to avoid false sharing.
    auto get_i = [&] (size_t i) -> size_t {return A[i].first;};
    auto s = delayed_seq<size_t>(n,get_i);
    get_bucket<decltype(s)> x(s, bits-1);
    auto get_buckets = delayed_seq<size_t>(n, x);
    sequence<T> B = sequence<T>::no_init(n);
    
    // first buckets based on hash using a counting sort
    sequence<size_t> bucket_offsets 
      = count_sort(A, B.slice(), get_buckets, num_buckets);
      //= integer_sort_2(A, B.slice(), x, bits);
    //t.next("sort");
    
    // note that this is cache line alligned
    size_t num_tables = num_buckets/2;
    size_t sz = n / num_tables;
    float factor = 1.2;
    if (sz < 128000) factor += (17 - log2_up(sz))*.15;
    uint table_size = (factor * sz);
    size_t total_table_size = table_size * num_tables;
    sequence<T> table = sequence<T>::no_init(total_table_size);
    sequence<size_t> sizes(num_tables + 1);
    //t.next("alloc");
	
    // now in parallel process each bucket sequentially
    auto hash_f = [&] (size_t i) {
      T* my_table = table.begin() + i * table_size;

      // clear tables
      for (size_t i = 0; i < table_size; i++)
	assign_uninitialized(my_table[i], T(empty, identity));

      // insert small bucket (ones with multiple different items)
      size_t start = bucket_offsets[i];
      size_t end = bucket_offsets[i+1];
      if ((end-start) > table_size)
	cout << "error in collect_reduce: " << (end-start) << ", " << table_size
	     << ", " << total_table_size << ", " << n << endl;
      for (size_t i = start; i < end; i++) {
	size_t idx = B[i].first;
	size_t j = ((uint) hash64_2(idx)) % table_size;
	while (my_table[j].first != empty &&
	       my_table[j].first != idx)
	  j = (j + 1 == table_size) ? 0 : j + 1;
	val_type v = my_table[j].second;
	my_table[j] = T(idx, monoid.f(v, B[i].second));
      }

      // now insert large bucket (ones with single item)
      size_t start_l = bucket_offsets[num_tables + i];
      size_t len = bucket_offsets[num_tables + i + 1] - start_l;
      if (len > 0) {
	auto f = [&] (size_t i) -> val_type {return B[i+start_l].second;};
	auto s = delayed_seq<val_type>(len, f);
	val_type x = reduce(s, monoid);
	size_t k = 0;
	while (my_table[k].first != empty) k++;
	my_table[k] = T(B[start_l].first, x);
      }

      // pack tables down to bottom
      size_t j=0;
      for (size_t i = 0; i < table_size; i++)
	if (my_table[i].first != empty) 
	  my_table[j++] = my_table[i];
      sizes[i] = j;

    };
    parallel_for(0, num_tables, hash_f, 1);
    //t.next("hash");

    sizes[num_tables] = 0;
    size_t total = scan_inplace(sizes.slice(), addm<size_t>());
    //t.next("scan");

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
    //t.next("copy");
    return result;
  }

  template <typename K, typename V, typename M>
  sequence<V> collect_reduce(sequence<std::pair<K,V>> const &A, size_t m, M monoid) {
    using P = std::pair<K,V>;
    auto get_index = [] (P v) {return v.first;};
    auto get_val = [] (P v) {return v.second;};
    return collect_reduce<V>(A, m, get_index, get_val, monoid);
  }
}
