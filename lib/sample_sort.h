// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010-2016 Guy Blelloch and Harsha Vardhan Simhadri and the PBBS team
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

// This file is basically the cache-oblivious sorting algorithm from:
//
// Low depth cache-oblivious algorithms.
// Guy E. Blelloch, Phillip B. Gibbons and  Harsha Vardhan Simhadri.
// Proc. ACM symposium on Parallelism in algorithms and architectures (SPAA), 2010

#pragma once
#include <math.h>
#include <stdio.h>
#include <string.h> 
#include "utilities.h"
#include "sequence_ops.h"
#include "quicksort.h"
//#include "merge_sort.h"
#include "transpose.h"

namespace pbbs {

  // the following parameters can be tuned
  constexpr const size_t QUICKSORT_THRESHOLD = 16384;
  constexpr const size_t OVER_SAMPLE = 8;

  // generates counts in Sc for the number of keys in Sa between consecutive
  // values of Sb
  // Sa and Sb must be sorted
  template<typename E, typename BinPred, typename s_size_t>
  void merge_seq (E* sA, E* sB, s_size_t* sC,
		  size_t lA, size_t lB, BinPred f) {
    if (lA==0 || lB==0) return;
    E *eA = sA+lA;
    E *eB = sB+lB;
    for (size_t i=0; i <= lB; i++) sC[i] = 0;
    while(1) {
      while (f(*sA, *sB)) {(*sC)++; if (++sA == eA) return;}
      sB++; sC++;
      if (sB == eB) break;
      if (!(f(*(sB-1),*sB))) {
	while (!f(*sB, *sA)) {(*sC)++; if (++sA == eA) return;}
	sB++; sC++;
	if (sB == eB) break;
      }
    } 
    *sC = eA-sA;
  }

  template<typename s_size_t = size_t, class Seq, typename BinPred>
  auto sample_sort_ (Seq A, const BinPred& f, bool inplace = false, bool stable = false)
    -> sequence<typename Seq::T> {
    using E = typename Seq::T;
    size_t n = A.size();

    //if (n == 0) return A;
    if (n < QUICKSORT_THRESHOLD) {
      sequence<E> B;
      auto y = [&] (size_t i) -> E {return A[i];};
      if (inplace) B = A.as_sequence();
      else B = sequence<E>(n, y);
      quicksort(B.begin(), n, f);
      return B;
    } else {
      timer t;
      size_t bucket_quotient = 5;
      size_t block_quotient = 5;
      if (is_pointer(A[0])) {
	bucket_quotient = 2;
	block_quotient = 3;
      } else if (sizeof(E) > 8) {
	bucket_quotient = block_quotient = 4;
      }
      size_t sqrt = (size_t) ceil(pow(n,0.5));
      size_t num_blocks = 1 << log2_up((sqrt/block_quotient) + 1);
      size_t block_size = ((n-1)/num_blocks) + 1;
      size_t num_buckets = (sqrt/bucket_quotient) + 1;
      size_t sample_set_size = num_buckets * OVER_SAMPLE;
      size_t m = num_blocks*num_buckets;
    
      E* sample_set = new_array<E>(sample_set_size);

      // generate "random" samples with oversampling
      auto hash_f = [&] (size_t j) {
	sample_set[j] = A[hash64(j)%n];};
      parallel_for(0, sample_set_size, hash_f, 1000);
      
      // sort the samples
      quicksort(sample_set, sample_set_size, f);

      // subselect samples at even stride
      E* pivots = new_array<E>(num_buckets-1);
      auto pivot_f = [&] (size_t k) {
	pivots[k] = sample_set[OVER_SAMPLE*k];};
      parallel_for(0, num_buckets-1, pivot_f, 1000);

      delete_array(sample_set,sample_set_size);

      sequence<E> Bs;
      if (inplace) Bs = A.as_sequence();
      else Bs = sequence<E>(new_array_no_init<E>(n,1), n);
      E* B = Bs.begin();
      E *C = new_array_no_init<E>(n,1);
      
      // sort each block and merge with samples to get counts for each bucket
      s_size_t *counts = new_array_no_init<s_size_t>(m+1,1);
      counts[m] = 0;
      parallel_for (0, num_blocks,  [&] (size_t i) {
	  size_t start = i * block_size;
	  size_t l = (i < num_blocks - 1) ? block_size : n - start;
	  if (!inplace)
	    for (size_t j = start;  j < start + l; j++) 
	      assign_uninitialized(B[j], A[j]);
	  if (stable) // if stable then use mergesort
	    {} //merge_sort_(sequence<E>(B+start, l), sequence<E>(C+start,l), f, 1);
	  else {
#if defined(OPENMP)
	    quicksort_serial(B+start, l, f);
#else
	    quicksort(B+start, l, f);
#endif
	  }
	  merge_seq(B + start, pivots, counts + i*num_buckets,
		    l, num_buckets-1, f);
	}, 1);
      
      // move data from blocks to buckets
      size_t* bucket_offsets = transpose_buckets(B, C, counts, n, block_size,
						 num_blocks, num_buckets);
      delete_array(counts, m);
      
      // sort within each bucket
      parallel_for (0, num_buckets, [&] (size_t i) {
	size_t start = bucket_offsets[i];
	size_t end = bucket_offsets[i+1];
	size_t l = end-start;

	// middle buckets need not be sorted if two consecutive pivots
	// are equal
	if (i == 0 || i == num_buckets - 1 || f(pivots[i-1],pivots[i])) {
	  if (stable)
	    {} //merge_sort_(sequence<E>(C+start,l), sequence<E>(B+start,l), f, 1);
	  else {
#if defined(OPENMP)
	    quicksort_serial(C+start, l, f);
#else
	    quicksort(C+start, l, f);
#endif
	  }
	}
	if (inplace) // move back to B
	  parallel_for (start, end, [&] (size_t k) {
	   move_uninitialized(B[k], C[k]);}, 1000);
	}, 1);
      
      delete_array(pivots,num_buckets-1);
      delete_array(bucket_offsets,num_buckets+1 );
      if (inplace) {free_array(C); return Bs;}
      else {free_array(B); return sequence<E>(C,n,true);}
    }
  }

  template<class Seq, typename BinPred>
  auto sample_sort (Seq &A, const BinPred& f, bool inplace = false, bool stable = false)
    -> sequence<typename Seq::T> {
    if (A.size() < ((size_t) 1) << 32)
      return sample_sort_<unsigned int>(A.slice(), f, inplace, stable);
    else return sample_sort_<size_t>(A.slice(), f, inplace, stable);
  }
    
  template<typename E, typename BinPred, typename s_size_t>
  void sample_sort (E* A, s_size_t n, const BinPred& f, bool stable = false) {
    sequence<E> B(A,A+n);
    sample_sort_<s_size_t>(B, f, true, stable);
  }
}
