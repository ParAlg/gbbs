// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010-2016 Guy Blelloch and Harsha Vardhan Simhadri and the PBBS
// team
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
// Proc. ACM symposium on Parallelism in algorithms and architectures (SPAA),
// 2010

#pragma once
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "bucket_sort.h"
#include "get_time.h"
#include "merge_sort.h"
#include "quicksort.h"
#include "sequence_ops.h"
#include "transpose.h"
#include "utilities.h"

namespace pbbs {

// the following parameters can be tuned
constexpr const size_t QUICKSORT_THRESHOLD = 16384;
constexpr const size_t OVER_SAMPLE = 8;

// generates counts in Sc for the number of keys in Sa between consecutive
// values of Sb
// Sa and Sb must be sorted
template <typename E, typename Compare, typename s_size_t>
void merge_seq(E* sA, E* sB, s_size_t* sC, size_t lA, size_t lB, Compare f) {
  if (lA == 0 || lB == 0) return;
  E* eA = sA + lA;
  E* eB = sB + lB;
  for (size_t i = 0; i <= lB; i++) sC[i] = 0;
  while (1) {
    while (f(*sA, *sB)) {
      (*sC)++;
      if (++sA == eA) return;
    }
    sB++;
    sC++;
    if (sB == eB) break;
    if (!(f(*(sB - 1), *sB))) {
      while (!f(*sB, *sA)) {
        (*sC)++;
        if (++sA == eA) return;
      }
      sB++;
      sC++;
      if (sB == eB) break;
    }
  }
  *sC = eA - sA;
}

template <typename T, typename Compare>
void seq_sort_inplace(range<T*> A, const Compare& less, bool stable) {
#if defined(OPENMP)
  quicksort_serial(A.begin(), A.size(), less);
#else
  if (((sizeof(T) > 8) || is_pointer(A[0])) && !stable)
    quicksort(A.begin(), A.size(), less);
  else
    bucket_sort(A, less, stable);
#endif
}

template <class Seq, class Iter, typename Compare>
void seq_sort_(Seq const& In, range<Iter> Out, const Compare& less,
               bool inplace = false, bool stable = false) {
  size_t l = In.size();
  if (inplace)
    for (size_t j = 0; j < l; j++) move_uninitialized(Out[j], In[j]);
  else
    for (size_t j = 0; j < l; j++) assign_uninitialized(Out[j], In[j]);
  seq_sort_inplace(Out, less, stable);
}

template <class Seq, class Iter, typename Compare>
void small_sort_(Seq const& In, range<Iter> Out, const Compare& less,
                 bool inplace = false, bool stable = false) {
  if (inplace)
    std::cout << "bad inplace arg in sort" << std::endl;
  else
    seq_sort_(In, Out, less, false, stable);
}

template <class Iter1, class Iter2, typename Compare>
void small_sort_(range<Iter1> In, range<Iter2> Out, const Compare& less,
                 bool inplace = false, bool stable = false) {
  if (inplace)
    seq_sort_inplace(In, less, stable);
  else
    seq_sort_(In, Out, less, false, stable);
}

// if inplace, then In and Out must be the same, i.e. it copies back to itsefl
// if inplace the copy constructor or assignment is never called on the elements
// if not inplace, then the copy constructor is called once per element
template <typename s_size_t = size_t, class SeqIn, class SeqOut,
          typename Compare>
void sample_sort_(SeqIn In, SeqOut Out, const Compare& less,
                  bool inplace = false, bool stable = false) {
  using T = typename SeqIn::value_type;
  size_t n = In.size();

  if (n < QUICKSORT_THRESHOLD) {
    small_sort_(In, Out, less, inplace, stable);
  } else {
    timer t("sample sort", false);
    // The larger these are, the more comparisons are done but less
    // overhead for the transpose.
    size_t bucket_quotient = 4;
    size_t block_quotient = 4;
    if (is_pointer(In[0])) {
      bucket_quotient = 2;
      block_quotient = 3;
    } else if (sizeof(T) > 8) {
      bucket_quotient = 3;
      block_quotient = 3;
    }
    size_t sqrt = (size_t)ceil(pow(n, 0.5));
    size_t num_blocks = 1 << log2_up((sqrt / block_quotient) + 1);
    size_t block_size = ((n - 1) / num_blocks) + 1;
    size_t num_buckets = (sqrt / bucket_quotient) + 1;
    size_t sample_set_size = num_buckets * OVER_SAMPLE;
    size_t m = num_blocks * num_buckets;

    // generate "random" samples with oversampling
    sequence<T> sample_set(sample_set_size,
                           [&](size_t i) { return In[hash64(i) % n]; });

    // sort the samples
    quicksort(sample_set.begin(), sample_set_size, less);

    // subselect samples at even stride
    sequence<T> pivots(num_buckets - 1,
                       [&](size_t i) { return sample_set[OVER_SAMPLE * i]; });

    sequence<T> Tmp = sequence<T>::no_init(n);
    t.next("head");

    // sort each block and merge with samples to get counts for each bucket
    s_size_t* counts = new_array_no_init<s_size_t>(m + 1, 1);
    counts[m] = 0;
    sliced_for(n, block_size, [&](size_t i, size_t start, size_t end) {
      seq_sort_(In.slice(start, end), Tmp.slice(start, end), less, inplace,
                stable);
      merge_seq(Tmp.begin() + start, pivots.begin(), counts + i * num_buckets,
                end - start, num_buckets - 1, less);
    });
    t.next("first sort");

    // move data from blocks to buckets
    size_t* bucket_offsets =
        transpose_buckets(Tmp.begin(), Out.begin(), counts, n, block_size,
                          num_blocks, num_buckets);
    t.next("transpose");
    // sort within each bucket
    parallel_for(0, num_buckets,
                 [&](size_t i) {
                   size_t start = bucket_offsets[i];
                   size_t end = bucket_offsets[i + 1];

                   // buckets need not be sorted if two consecutive pivots are
                   // equal
                   if (i == 0 || i == num_buckets - 1 ||
                       less(pivots[i - 1], pivots[i])) {
                     seq_sort_inplace(Out.slice(start, end), less, stable);
                   }
                 },
                 1);
    t.next("second sort");
    delete_array(bucket_offsets, num_buckets + 1);
  }
}

template <class Seq, typename Compare>
auto sample_sort(Seq const& A, const Compare& less, bool stable = false)
    -> sequence<typename Seq::value_type> {
  using T = typename Seq::value_type;
  sequence<T> R = sequence<T>::no_init(A.size());
  if (A.size() < ((size_t)1) << 32)
    sample_sort_<unsigned int>(A.slice(), R.slice(), less, false, stable);
  else
    sample_sort_<size_t>(A.slice(), R.slice(), less, false, stable);
  return R;
}

template <class Iter, typename Compare>
void sample_sort_inplace(range<Iter> A, const Compare& less,
                         bool stable = false) {
  if (A.size() < ((size_t)1) << 32)
    sample_sort_<unsigned int>(A.slice(), A.slice(), less, true, stable);
  else
    sample_sort_<size_t>(A.slice(), A.slice(), less, true, stable);
}

template <typename E, typename Compare, typename s_size_t>
void sample_sort(E* A, s_size_t n, const Compare& less, bool stable) {
  range<E*> B(A, A + n);
  sample_sort_inplace(B, less, stable);
}

}  // namespace pbbs
