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

namespace pbbs {

  // the following parameters can be tuned
  constexpr const size_t CR_SEQ_THRESHOLD = 8192;

  template<class OutVal, class InSeq, class KeySeq, class Identity, class F>
  void _seq_collect_reduce(InSeq In, KeySeq Keys, OutVal* Out, size_t num_buckets, Identity I, F f) {
    size_t n = In.size();
    for (size_t i = 0; i < num_buckets; i++) Out[i] = I;
    for (size_t j = 0; j < n; j++) {
      size_t k = Keys[j];
      Out[k] = f(Out[k],In[j]);
    }
  }

  template<class OutVal, class InSeq, class KeySeq, class Identity, class F>
  sequence<OutVal> seq_collect_reduce(InSeq In, KeySeq Keys,
				      size_t num_buckets, Identity I, F f) {
    OutVal* Out = new_array<OutVal>(num_buckets);
    _seq_collect_reduce(In, Keys, Out, num_buckets, I, f);
    return sequence<OutVal>(Out,num_buckets);
  }

  // this one is for few buckets (e.g. less than 2^16)
  //  In is the sequence of values to be reduced
  //  Keys is the sequence of bucket numbers
  //  num_buckets is the number of buckets (all keys need to be less)
  //  I is the identity for add
  //  f is a function that adds two values
  template<class OutVal, class InSeq, class KeySeq, class Identity, class F>
  sequence<OutVal> collect_reduce(InSeq In, KeySeq Keys,
				  size_t num_buckets, Identity I, F f) {
    size_t n = In.size();
    timer t;
    t.start();
    
    // pad to 16 buckets to avoid false sharing (does not affect results)
    num_buckets = std::max(num_buckets, (size_t) 16);

    //size_t num_blocks = ceil(pow(n/num_buckets,0.5));
    size_t num_threads = num_workers();
    size_t num_blocks = std::min(4*num_threads, n/num_buckets/64);

    num_blocks = 1 << log2_up(num_blocks);

    OutVal* Out = new_array<OutVal>(num_buckets);

    // if insufficient parallelism, sort sequentially
    if (n < CR_SEQ_THRESHOLD || num_blocks == 1 || num_threads == 1)
      return seq_collect_reduce<OutVal>(In,Keys,num_buckets, I, f);

    size_t block_size = ((n-1)/num_blocks) + 1;
    size_t m = num_blocks * num_buckets;

    OutVal *OutM = new_array<OutVal>(m);
    
    //parallel_for (size_t i = 0; i < num_blocks; ++i) {
    auto block_f = [&] (size_t i) {
      size_t start = std::min(i * block_size, n);
      size_t end =  std::min(start + block_size, n);
      _seq_collect_reduce<OutVal>(In.slice(start,end), Keys.slice(start,end),
				  OutM + i*num_buckets, num_buckets, I, f);
    };
    par_for (0, num_blocks, 1, block_f);
    
    auto sum_buckets = [&] (size_t i) {
      OutVal O = I;
      for (size_t j = 0; j < num_blocks; j++)
	O = f(O, OutM[i + j*num_buckets]);
      Out[i] = O;
    };
    par_for(0, num_buckets, 1, sum_buckets);
    delete_array(OutM, m);
    return sequence<OutVal>(Out,num_buckets);
  }

  // this one is for many buckets but no more than len(A) buckets
  //      (e.g. more than 2^16)
  //  A is the input sequence
  //  m is the number of buckets
  //  get_index is a function that gets the bucket from an element of A
  //  get_val is a function that gets the value from an element of A
  //  identity is the identity for add
  //  add is a function that adds two values
  template <typename OT, typename Seq, typename F, typename G,
    typename I, typename H>
  sequence<OT> collect_reduce(Seq A, size_t m,
			      F get_index, G get_val, I identity, H add) {
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
    auto s = make_sequence<size_t>(n,get_i);
    get_bucket<decltype(s)> x(s, bits-1);
    auto get_buckets = make_sequence<size_t>(n, x);

    // first buckets based on hash using a counting sort
    sequence<size_t> bucket_offsets 
      = count_sort(A, A, get_buckets, num_buckets);

    // note that this is cache line alligned
    sequence<OT> sums(m, identity);
       
    // now process each bucket in parallel
    auto bucket_f = [&] (size_t i) {
      size_t start = bucket_offsets[i];
      size_t end = bucket_offsets[i+1];

      // small buckets have indices in bottom half
      if (i < num_buckets/2)
	for (size_t i = start; i < end; i++) {
	  size_t j = get_index(A[i]);
	  sums[j] = add(sums[j], get_val(A[i]));
	}

      // large buckets have indices in top half
      else if (end > start) {
	auto x = [&] (size_t i) {return get_val(A[i]);};
	auto vals = make_sequence<size_t>(n, x);
	sums[get_index(A[i])] = reduce(vals, add);
      }
    };
    par_for(0, num_buckets, 1, bucket_f);
    return sums;
  }

  // this one is for more buckets than the length of A
  //  A is the input sequence
  //  get_index is a function that gets the bucket from an element of A
  //  get_val is a function that gets the value from an element of A
  //  The template parameter R is a reducer that must supply:
  //      identity and add
  //  The template parameter P is a pair type:  pair<key_time,val_type>
  //  The other template parameters should be inferred
  template <typename R, typename P, typename Seq, typename F, typename G>
  sequence<P> collect_reduce_pair(Seq A, F get_index, G get_val) {
    using key_type = typename P::first_type;
    using val_type = typename P::second_type;
    key_type empty = (key_type) -1;
    val_type identity = R::identity();
    
    timer t;
    t.start();
    size_t n = A.size();

    if (n < 1000) {
      sequence<P> x(n, [&] (size_t i) {return P(get_index(A[i]), get_val(A[i]));});
      auto cmp = [] (P a, P b) {return a.first < b.first;};
      sequence<P> y = sample_sort(x, cmp);
      size_t j = 0;
      for (size_t i = 1; i < n; i++) {
	if (A[i].first == A[j].first)
	  A[j].second = R::add(A[j].second, A[i].second);
	else A[j++] = A[i];
      };
      return sequence<P>(j, [&] (size_t i) {return A[i];});
    }
      
    size_t bits;
    if (n < (1 << 19)) bits = (log2_up(n))/2;
    else if (n < (1 << 23)) bits = (log2_up(n) - 3)/2;
    else if (n < (1 << 27)) bits = (log2_up(n) - 5)/2;
    // for large n selected so each bucket fits into cache
    else bits = (log2_up(n) - 16);
    size_t num_buckets = (1<<bits);

    // Returns a map (hash) from key to bucket.
    // Keys with many elements (big) have their own bucket while
    // others share a bucket.
    // Keys that share low 4 bits get same bucket unless big.
    // This is to avoid false sharing.
    auto get_i = [&] (size_t i) -> size_t {return get_index(A[i]);};
    auto s = make_sequence<size_t>(n,get_i);
    get_bucket<decltype(s)> x(s, bits-1);
    auto get_buckets = make_sequence<size_t>(n, x);

    // first buckets based on hash using a counting sort
    sequence<size_t> bucket_offsets 
      = count_sort(A, A, get_buckets, num_buckets);
    //t.next("sort");
    
    // note that this is cache line alligned
    size_t num_tables = num_buckets/2;
    //size_t table_size = ((size_t) 1) << table_bits;
    size_t sz = n / num_tables;
    float factor = 1.2;
    if (sz < 128000) factor += (17 - log2_up(sz))*.15;
    uint table_size = (factor * sz);
    //size_t mask = table_size-1;
    size_t total_table_size = table_size * num_tables;
    sequence<P> table
      = sequence<P>::alloc_no_init(total_table_size);
    sequence<size_t> sizes(num_tables + 1);
    //t.next("alloc");
	
    // now in parallel process each bucket sequentially
    auto hash_f = [&] (size_t i) {
      P* my_table = table.as_array() + i * table_size;

      // clear tables
      for (size_t i = 0; i < table_size; i++)
	assign_uninitialized(my_table[i], P(empty, identity));

      // insert small bucket (ones with multiple different items)
      size_t start = bucket_offsets[i];
      size_t end = bucket_offsets[i+1];
      if ((end-start) > table_size)
	cout << "yikes: " << (end-start) << ", " << table_size << ", " << total_table_size << ", " << n << endl;
      for (size_t i = start; i < end; i++) {
	size_t idx = get_index(A[i]);
	size_t j = ((uint) hash64_2(idx)) % table_size;
	while (my_table[j].first != empty &&
	       my_table[j].first != idx)
	  j = (j + 1 == table_size) ? 0 : j + 1;
	val_type v = my_table[j].second;
	my_table[j] = P(idx, R::add(v, get_val(A[i])));
      }

      // now insert large bucket (ones with single item)
      size_t start_l = bucket_offsets[num_tables + i];
      size_t len = bucket_offsets[num_tables + i + 1] - start_l;
      if (len > 0) {
	auto f = [&] (size_t i) -> val_type {return get_val(A[i+start_l]);};
	val_type x = reduce(make_sequence<val_type>(len, f), R::add);
	size_t k = 0;
	while (my_table[k].first != empty) k++;
	my_table[k] = P(get_index(A[start_l]), x);
      }

      // pack tables down to bottom
      size_t j=0;
      for (size_t i = 0; i < table_size; i++)
	if (my_table[i].first != empty) 
	  my_table[j++] = my_table[i];
      sizes[i] = j;

    };
    par_for(0, num_tables, 1, hash_f);
    //t.next("hash");

    sizes[num_tables] = 0;
    size_t total = scan_add(sizes, sizes);
    //t.next("scan");

    // copy packed tables into contiguous result
    sequence<P> result = sequence<P>::alloc_no_init(total);
    auto copy_f = [&] (size_t i) {
      size_t d_offset = sizes[i];
      size_t s_offset = i * table_size;
      size_t len = sizes[i+1] - sizes[i];
      for (size_t j = 0; j < len; j++)
	move_uninitialized(result[d_offset+j], table[s_offset+j]);
    };
    par_for(0, num_tables, 1, copy_f);
    //t.next("copy");
    return result;
  }

  template <typename R, typename Seq>
  Seq collect_reduce_pair(Seq A) {
    using P = typename Seq::T;
    auto get_index = [] (P v) {return v.first;};
    auto get_val = [] (P v) {return v.second;};
    return collect_reduce_pair<R,P>(A, get_index, get_val);
  }
    
  template <typename K, typename V, typename I, typename H>
  sequence<V> collect_reduce(sequence<std::pair<K,V>> A, size_t m, I identity, H add) {
    using P = std::pair<K,V>;
    auto get_index = [] (P v) {return v.first;};
    auto get_val = [] (P v) {return v.second;};
    return collect_reduce<V>(A, m, get_index, get_val, identity, add);
  }
}
