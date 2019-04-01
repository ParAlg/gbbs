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

// TODO
// Make sure works for inplace or not with regards to move_uninitialized

namespace pbbs {

  // the following parameters can be tuned
  constexpr const size_t SEQ_THRESHOLD = 8192;
  constexpr const size_t BUCKET_FACTOR = 32;
  constexpr const size_t LOW_BUCKET_FACTOR = 16;

  // Sequential internal version
  template <typename b_size_t, typename s_size_t, typename InSeq, typename KeySeq>
  void _seq_count_sort(InSeq In, typename InSeq::T* Out, KeySeq Keys,
		       s_size_t* counts, size_t num_buckets) {
    size_t n = In.size();
    b_size_t tmp[n];

    //sequence<b_size_t> tmp(n);

    for (size_t i = 0; i < num_buckets; i++)
      counts[i] = 0;
    for (size_t j = 0; j < n; j++) {
      size_t k = tmp[j] = Keys[j];
      if (k >= num_buckets) abort();
      counts[k]++;
    }
    
    // generate offsets
    size_t s = 0;
    for (size_t i = 0; i < num_buckets; i++) {
      size_t t = counts[i];
      counts[i] = s;
      s += t;
    }
    for (size_t j = 0; j < n; j++) {
      size_t k = counts[tmp[j]]++;
      // needed for types with self defined assignment or initialization
      // otherwise equivalent to: Out[k+start] = In[j+start];
      move_uninitialized(Out[k], In[j]);
    }
    // convert offsets back to counts
    s = 0;
    for (size_t i = 0; i < num_buckets; i++) {
      size_t t = counts[i];
      counts[i] = t - s;
      s = t;
    }
  }

   // Sequential internal version
  template <typename s_size_t, typename InSeq, typename KeySeq>
  void _seq_count(InSeq In, KeySeq Keys,
		  s_size_t* counts, size_t num_buckets) {
    size_t n = In.size();

    for (size_t i = 0; i < num_buckets; i++)
      counts[i] = 0;
    for (size_t j = 0; j < n; j++) {
      size_t k = Keys[j];
      if (k >= num_buckets) abort();
      counts[k]++;
    }
  }

  template <typename s_size_t, typename InSeq, typename KeySeq>
  void _seq_write(InSeq In, typename InSeq::T* Out, KeySeq Keys,
		  s_size_t* offsets, size_t num_buckets) {
  
    for (size_t j = 0; j < In.size(); j++) {
      size_t k = offsets[Keys[j]]++;
      // needed for types with self defined assignment or initialization
      // otherwise equivalent to: Out[k+start] = In[j+start];
      move_uninitialized(Out[k], In[j]);
    }
  }

  // Sequential version
  template <typename b_size_t, typename InS, typename OutS, typename KeyS>
  sequence<size_t> seq_count_sort(InS& In, OutS& Out, KeyS& Keys, size_t num_buckets) {
    using T = typename InS::T;
    size_t n = In.size();

    size_t* counts = new_array_no_init<size_t>(num_buckets+1);
    T* B = new_array_no_init<T>(n);
    _seq_count_sort<b_size_t,size_t>(In, B, Keys, counts, num_buckets);
    for (size_t i = 0; i < n ; i++)
      Out[i] = B[i]; 
    free_array(B);
    size_t c = 0;
    for (size_t i=0; i < num_buckets; i++) {
      size_t x = counts[i];
      counts[i] = c;
      c = c + x;
    }
    counts[num_buckets] = n;
    return sequence<size_t>(counts,num_buckets+1);
  }

  // Parallel internal version
  template <typename b_size_t, typename s_size_t, 
    typename InS, typename OutS, typename KeyS>
  sequence<size_t> _count_sort(InS& In, OutS& Out, KeyS& Keys,
			       size_t num_buckets) {
    timer t;
    t.start();
    using T = typename InS::T;
    size_t n = In.size();
    size_t num_threads = num_workers();

    // pad to 16 buckets to avoid false sharing (does not affect results)
    num_buckets = std::max(num_buckets, (size_t) 16);

    // if not given, then use heuristic to choose num_blocks
    size_t sqrt = (size_t) ceil(pow(n,0.5));
    size_t num_blocks = 
      (size_t) (n < (1<<24)) ? (sqrt/16) : ((n < (1<<28)) ? sqrt/2 : sqrt);
    if (2*num_blocks < num_threads) num_blocks *= 2;

    num_blocks = 1 << log2_up(num_blocks);

    // if insufficient parallelism, sort sequentially
    if (n < SEQ_THRESHOLD || num_blocks == 1 || num_threads == 1) {
      return seq_count_sort<b_size_t>(In,Out,Keys,num_buckets);}

    size_t block_size = ((n-1)/num_blocks) + 1;
    size_t m = num_blocks * num_buckets;
    
    T *B = new_array_no_init<T>(n);
    s_size_t *counts = new_array_no_init<s_size_t>(m,1);
    //t.next("head");

    // sort each block
    auto block_f = [&] (size_t i) {
      s_size_t start = std::min(i * block_size, n);
      s_size_t end =  std::min(start + block_size, n);
      _seq_count_sort<b_size_t,s_size_t>(In.slice(start,end), B+start, 
					 Keys.slice(start,end),
					 counts + i*num_buckets, num_buckets);
    };
    parallel_for(0, num_blocks, block_f, 1);
    //t.next("count");
    T* C = Out.begin();
    size_t* bucket_offsets = transpose_buckets(B, C,
					       counts, n, block_size,
					       num_blocks, num_buckets);
    //t.next("transpose");
    free_array(B);
    return sequence<size_t>(bucket_offsets,num_buckets+1);
  }

  // Sequential version
  template <typename InS, typename OutS, typename KeyS>
  sequence<size_t> seq_count_sort2(InS& In, OutS& Out, KeyS& Keys,
				   size_t num_buckets) {
    sequence<size_t> counts(num_buckets+1);

    // count size of each bucket
    _seq_count(In, Keys, counts.begin(), num_buckets);

    // generate offsets for buckets
    size_t s = 0;
    for (size_t i = 0; i < num_buckets; i++) {
      size_t t = counts[i];
      counts[i] = s;
      s += t;
    }
    counts[num_buckets] = s;

    // send to destination
    _seq_write(In, Out.begin(), Keys, counts.begin(), num_buckets);
    return counts;
  }

  // Parallel internal version
  template <typename b_size_t, typename s_size_t, 
	    typename InS, typename OutS, typename KeyS>
  sequence<size_t> _count_sort2(InS& In, OutS& Out, KeyS& Keys,
				size_t num_buckets) {
    timer t;
    t.start();
    using T = typename InS::T;
    size_t n = In.size();
    size_t num_threads = num_workers();
    cout << num_threads << endl;

    // pad to 16 buckets to avoid false sharing (does not affect results)
    num_buckets = std::max(num_buckets, (size_t) 16);

    // if not given, then use heuristic to choose num_blocks
    size_t sqrt = (size_t) ceil(pow(n,0.5));
    size_t num_blocks = 
      (size_t) (n < (1<<24)) ? (sqrt/16) : ((n < (1<<28)) ? sqrt/2 : sqrt);
    if (2*num_blocks < num_threads) num_blocks *= 2;
    if (sizeof(T) <= 4) num_blocks = num_blocks/2;

    //num_blocks = 1 << log2_up(num_blocks);
    //cout << num_blocks << endl;
    //num_blocks = 1024;

    // if insufficient parallelism, sort sequentially
    if (n < SEQ_THRESHOLD || num_blocks == 1 || num_threads == 1) {
      return seq_count_sort2(In,Out,Keys,num_buckets);}

    size_t block_size = ((n-1)/num_blocks) + 1;
    size_t m = num_blocks * num_buckets;
    
    //T *B = new_array_no_init<T>(n);
    s_size_t *counts = new_array_no_init<s_size_t>(m,1);
    if (n > 1000000000) t.next("head");

    // sort each block
    parallel_for(0, num_blocks,  [&] (size_t i) {
	s_size_t start = std::min(i * block_size, n);
	s_size_t end =  std::min(start + block_size, n);
	_seq_count(In.slice(start,end), Keys.slice(start,end),
		   counts + i*num_buckets, num_buckets);
      },1);

    if (n > 1000000000) t.next("count");

    sequence<size_t> bucket_offsets = sequence<size_t>::no_init(num_buckets+1);
    parallel_for(0, num_buckets, [&] (size_t i) {
	size_t v = 0;
	for (size_t j= 0; j < num_blocks; j++) 
	  v += counts[j*num_buckets + i];
	bucket_offsets[i] = v;
      }, 1 + 1024/num_blocks);
    bucket_offsets[num_buckets] = 0;
    size_t total = scan_inplace(bucket_offsets, addm<size_t>());
    if (total != n) abort();
    sequence<s_size_t> dest_offsets = sequence<s_size_t>::no_init(num_blocks*num_buckets);
    parallel_for(0, num_buckets, [&] (size_t i) {
	size_t v = bucket_offsets[i];
	size_t start = i * num_blocks;
	for (size_t j= 0; j < num_blocks; j++) {
	  dest_offsets[start+j] = v;
	  v += counts[j*num_buckets + i];
	}
      }, 1 + 1024/num_blocks);

    parallel_for(0, num_blocks, [&] (size_t i) {
	size_t start = i * num_buckets;
	for (size_t j= 0; j < num_buckets; j++)
	  counts[start+j] = dest_offsets[j*num_blocks + i];
      }, 1 + 1024/num_buckets);

    // transpose<s_size_t>(counts, dest_offsets.begin()).trans(num_blocks,
    // 							    num_buckets);
    // size_t sum = scan_inplace(dest_offsets, addm<s_size_t>());
    // if (sum != n) abort();
    // transpose<s_size_t>(dest_offsets.begin(), counts).trans(num_buckets,
    // 							    num_blocks);
        if (n > 1000000000) t.next("scan");

    parallel_for(0, num_blocks,  [&] (size_t i) {
	s_size_t start = std::min(i * block_size, n);
	s_size_t end =  std::min(start + block_size, n);
	_seq_write(In.slice(start,end), Out.begin(),
		   Keys.slice(start,end),
		   counts + i*num_buckets, num_buckets);
      },1);

    if (n > 1000000000) t.next("move");
    
    // for (s_size_t i=0; i < num_buckets; i++) {
    //   bucket_offsets[i] = dest_offsets[i*num_blocks];
    //   //cout << i << ", " << bucket_offsets[i] << endl;
    // }
    // // last element is the total size n
    // bucket_offsets[num_buckets] = n;

    //t.next("transpose");
    //free_array(B);
    if (n > 1000000000) {
      //for (size_t i=0; i < 10; i++) cout << bucket_offsets[i] << endl;
    }
    return bucket_offsets;
  }

  template <typename s_size_t, typename InS, typename OutS, typename KeyS>
  sequence<size_t> _count_sort_size(InS& In, OutS& Out, KeyS& Keys, size_t num_buckets) {
    if (num_buckets <= 256)
      if (true)
	return _count_sort2<uint8_t,s_size_t>(In, Out, Keys, num_buckets);
      else return _count_sort<uint8_t,s_size_t>(In, Out, Keys, num_buckets);
    else if  (num_buckets <= (1 << 16))
      return _count_sort2<uint16_t,s_size_t>(In, Out, Keys, num_buckets);
    else
      return _count_sort<s_size_t,s_size_t>(In, Out, Keys, num_buckets);
  }

  // Parallel version
  template <typename InS, typename OutS, typename KeyS>
  sequence<size_t> count_sort(InS& In, OutS& Out, KeyS& Keys, size_t num_buckets) {
    size_t n = In.size();
    size_t max32 = ((size_t) 1) << 32;
    if (n < max32 && num_buckets < max32)
      // use 4-byte counters when larger ones not needed
      return _count_sort_size<uint32_t>(In, Out, Keys, num_buckets);
    return _count_sort_size<size_t>(In, Out, Keys, num_buckets);
  }
}
