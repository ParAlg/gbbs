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
#include "merge_sort.h"
#include "transpose.h"
#include "bucket_sort.h"
#include "get_time.h"

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

  template<class Iter, class SeqB, typename BinPred>
  void sort_small_(range<Iter> A, SeqB B, const BinPred& f,
		   bool inplace = false, bool stable = false) {
    using T = typename SeqB::value_type;
    size_t n = A.size();
    if (stable) {
      if (inplace) {
	sequence<T> C(n);
	merge_sort_(A.slice(), C.slice(), f, true);
      } else ; merge_sort_(A.slice(), B.slice(), f, false);
    } else {
      if (!inplace) 
	parallel_for(0, n, [&] (size_t i) { B[i] = A[i];});
      quicksort(B.begin(), n, f);
    }
  }

  template<class SeqA, class SeqB, typename BinPred>
  void sort_small_(SeqA A, SeqB B, const BinPred& f,
		   bool inplace = false, bool stable = false) {
    if (inplace) std::cout << "bad inplace arg in sort" << std::endl;
    using T = typename SeqA::value_type;
    size_t n = A.size();
    if (stable) {
      sequence<T> C(A);
      merge_sort_(C.slice(), B.slice(), f, false);
    } else {
      parallel_for(0, n, [&] (size_t i) { B[i] = A[i];});
      quicksort(B.begin(), n, f);
    }
  }

  template<typename E, typename BinPred, typename s_size_t>
  void sample_sort (E* A, s_size_t n, const BinPred& f, bool stable = false);

  template<typename E, typename BinPred>
  void seq_sort_inplace(range<E*> A, BinPred f, bool stable) {
#if defined(OPENMP)
    quicksort_serial(A.begin(), A.size(), f);
#else
    if (((sizeof(E) > 8) || is_pointer(A[0])) && !stable) 
      quicksort(A.begin(), A.size(), f);
    else bucket_sort(A, f, stable);
#endif
  }
  
  template<typename s_size_t = size_t, class SeqA, class SeqB, typename BinPred>
  void sample_sort_ (SeqA A, SeqB B, const BinPred& f,
		     bool inplace = false, bool stable = false) {
    using T = typename SeqA::value_type;
    size_t n = A.size();

    //if (n == 0) return A;
    if (n < QUICKSORT_THRESHOLD) {
      sort_small_(A,B, f, inplace, stable);
    } else {
      timer t("sample sort", false);
      size_t bucket_quotient = 4;
      size_t block_quotient = 4;
      if (is_pointer(A[0])) {
	bucket_quotient = 2;
	block_quotient = 3;
      } else if (sizeof(T) > 8) {
	bucket_quotient = 3;
	block_quotient = 3;
      }
      size_t sqrt = (size_t) ceil(pow(n,0.5));
      size_t num_blocks = 1 << log2_up((sqrt/block_quotient) + 1);
      size_t block_size = ((n-1)/num_blocks) + 1;
      size_t num_buckets = (sqrt/bucket_quotient) + 1;
      size_t sample_set_size = num_buckets * OVER_SAMPLE;
      size_t m = num_blocks*num_buckets;
      
      // generate "random" samples with oversampling
      sequence<T> sample_set(sample_set_size, [&] (size_t i) {
	  return A[hash64(i)%n];});
      
      // sort the samples
      quicksort(sample_set.begin(), sample_set_size, f);

      // subselect samples at even stride
      sequence<T> pivots(num_buckets-1, [&] (size_t i) {
	  return sample_set[OVER_SAMPLE*i];});

      sequence<T> C = sequence<T>::no_init(n);
      t.next("head");
      
      // sort each block and merge with samples to get counts for each bucket
      s_size_t *counts = new_array_no_init<s_size_t>(m+1,1);
      counts[m] = 0;
      parallel_for (0, num_blocks,  [&] (size_t i) {
	  size_t start = std::min(n, i * block_size);
	  size_t end = std::min(n, start + block_size);
	  size_t l = end-start;
	  if (inplace)
	    for (size_t j = start;  j < start + l; j++) 
	      move_uninitialized(C[j], A[j]);
	  else
	    for (size_t j = start;  j < start + l; j++) 
	      assign_uninitialized(C[j], A[j]);
	  seq_sort_inplace(C.slice(start,end), f, stable);
	  merge_seq(C.begin() + start, pivots.begin(), counts + i*num_buckets,
		    l, num_buckets-1, f);
	}, 1);
      t.next("first sort");

      // move data from blocks to buckets
      size_t* bucket_offsets = transpose_buckets(C.begin(), B.begin(),
						 counts, n, block_size,
						 num_blocks, num_buckets);
      t.next("transpose");
      
      // sort within each bucket
      parallel_for (0, num_buckets, [&] (size_t i) {
	  size_t start = bucket_offsets[i];
	  size_t end = bucket_offsets[i+1];

	  // buckets need not be sorted if two consecutive pivots are equal
	  if (i == 0 || i == num_buckets - 1 || f(pivots[i-1],pivots[i])) {
	    seq_sort_inplace(B.slice(start,end), f, stable);
	  }
	},1);
      t.next("second sort");
      delete_array(bucket_offsets,num_buckets+1 );
    }
  }

  template<class Seq, typename BinPred>
  auto sample_sort (Seq const &A, const BinPred& f, bool stable = false)
    -> sequence<typename Seq::value_type> {
    using T = typename Seq::value_type;
    sequence<T> R = sequence<T>::no_init(A.size());
    if (A.size() < ((size_t) 1) << 32)
      sample_sort_<unsigned int>(A.slice(), R.slice(), f, false, stable);
    else sample_sort_<size_t>(A.slice(), R.slice(), f, false, stable);
    return R;
  }

  template<class Iter, typename BinPred>
  void sample_sort_inplace (range<Iter> A, const BinPred& f, bool stable = false) {
    if (A.size() < ((size_t) 1) << 32)
      sample_sort_<unsigned int>(A.slice(), A.slice(), f, true, stable);
    else sample_sort_<size_t>(A.slice(), A.slice(), f, true, stable);
  }
    
  template<typename E, typename BinPred, typename s_size_t>
  void sample_sort (E* A, s_size_t n, const BinPred& f, bool stable) {
    range<E*> B(A,A+n);
    sample_sort_inplace(B, f, stable);
  }
}
