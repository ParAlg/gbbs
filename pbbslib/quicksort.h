// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and the PBBS team
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

#include <algorithm>
#include "sequence_ops.h"
#include "utilities.h"

namespace pbbs {

// use different parameters for pointer and non-pointer types
// and depending on size
template <typename E>
bool is_pointer(E x) {
  return 0;
}
template <typename E>
bool is_pointer(E* x) {
  return 1;
}
// template<typename E, typename V> bool is_pointer(std::pair<E*,V> x) {return
// 1;}

template <class T>
bool base_case(T* x, size_t n) {
  bool large = is_pointer(x[0]) || (sizeof(x) > 8);
  return large ? (n < 16) : (n < 24);
}

template <class E, class BinPred>
void insertion_sort(E* A, size_t n, const BinPred& f) {
  for (size_t i = 0; i < n; i++) {
    E v = A[i];
    E* B = A + i;
    while (--B >= A && f(v, *B)) *(B + 1) = *B;
    *(B + 1) = v;
  }
}

// sorts 5 elements taken at even stride and puts them at the front
template <class E, class BinPred>
void sort5(E* A, size_t n, const BinPred& f) {
  size_t size = 5;
  size_t m = n / (size + 1);
  for (size_t l = 0; l < size; l++) std::swap(A[l], A[m * (l + 1)]);
  insertion_sort(A, size, f);
}

// splits based on two pivots (p1 and p2) into 3 parts:
//    less than p1, greater than p2, and the rest in the middle
// If the pivots are the same, returns a true flag to indicate middle need not
// be sorted
template <class E, class BinPred>
std::tuple<E*, E*, bool> split3(E* A, size_t n, const BinPred& f) {
  sort5(A, n, f);
  E p1 = A[1];
  E p2 = A[3];
  if (!f(A[0], A[1])) p1 = p2;  // if few elements less than p1, then set to p2
  if (!f(A[3], A[4]))
    p2 = p1;  // if few elements greater than p2, then set to p1
  E* L = A;
  E* R = A + n - 1;
  // set up initial invariants
  while (f(*L, p1)) L++;
  while (f(p2, *R)) R--;
  E* M = L;
  // invariants:
  //  below L is less than p1,
  //  above R is greatert than p2
  //  between L and M are between p1 and p2 inclusive
  //  between M and R are unprocessed
  while (M <= R) {
    // E mid = *M;
    if (f(*M, p1)) {
      std::swap(*M, *L);
      L++;
    } else if (f(p2, *M)) {
      std::swap(*M, *R);
      if (f(*M, p1)) {
        std::swap(*L, *M);
        L++;
      }
      // if (f(*R,p1)) {*M = *L; *L = *R; L++;}
      // else {*M = *R;}
      //*R = mid;
      R--;
      while (f(p2, *R)) R--;
    }
    M++;
  }
  return std::make_tuple(L, M, !f(p1, p2));
}

template <class E, class BinPred>
void quicksort_serial(E* A, size_t n, const BinPred& f) {
  while (!base_case(A, n)) {
    E* L;
    E* M;
    bool mid_eq;
    std::tie(L, M, mid_eq) = split3(A, n, f);
    if (!mid_eq) quicksort_serial(L, M - L, f);
    quicksort_serial(M, A + n - M, f);
    n = L - A;
  }
  insertion_sort(A, n, f);
}

template <class E, class BinPred>
void quicksort(E* A, size_t n, const BinPred& f) {
  if (n < (1 << 10))
    quicksort_serial(A, n, f);
  else {
    E* L;
    E* M;
    bool mid_eq;
    std::tie(L, M, mid_eq) = split3(A, n, f);
    auto left = [&]() { quicksort(A, L - A, f); };
    auto mid = [&]() { quicksort(L, M - L, f); };
    auto right = [&]() { quicksort(M, A + n - M, f); };

    if (!mid_eq)
      par_do3(left, mid, right);
    else
      par_do(left, right);
  }
}

template <class Range, class BinPred>
void quicksort(Range A, const BinPred& f) {
  quicksort(A.begin(), A.size(), f);
}

//// Fully Parallel version below here

template <class SeqA, class BinPred>
std::tuple<size_t, size_t, bool> p_split3(SeqA const& A,
                                          range<typename SeqA::value_type*> B,
                                          const BinPred& f) {
  using E = typename SeqA::value_type;
  size_t n = A.size();
  sort5(A.begin(), n, f);
  E p1 = A[1];
  E p2 = A[3];
  if (!f(A[0], A[1])) p1 = p2;  // if few elements less than p1, then set to p2
  if (!f(A[3], A[4]))
    p2 = p1;  // if few elements greater than p2, then set to p1
  auto flag = [&](size_t i) { return f(A[i], p1) ? 0 : f(p2, A[i]) ? 2 : 1; };
  auto r = split_three(A, B.slice(), delayed_seq<unsigned char>(n, flag),
                       fl_conservative);
  return std::make_tuple(r.first, r.first + r.second, !f(p1, p2));
  // sequence<size_t> r = count_sort(A, B.slice(),
  //				    delayed_seq<unsigned char>(n, flag), 3,
  //true);
  // return std::make_tuple(r[0],r[0]+r[1], !f(p1,p2));
}

// The fully parallel version copies back and forth between two arrays
// inplace: if true then result is put back in In
//     and Out is just used as temp space
//     otherwise result is in Out
//     In and Out cannot be the same (Out is needed as temp space)
// cut_size: is when to revert to  quicksort.
//    If -1 then it uses a default based on number of threads
template <class Iter, class F>
void p_quicksort_(range<Iter> In, range<Iter> Out, const F& f,
                  bool inplace = false, long cut_size = -1) {
  size_t n = In.size();
  if (cut_size == -1)
    cut_size = std::max<long>((3 * n) / num_workers(), (1 << 14));
  if (n < (size_t)cut_size) {
    quicksort(In.begin(), n, f);
    auto copy_out = [&](size_t i) { Out[i] = In[i]; };
    if (!inplace) parallel_for(0, n, copy_out, 2000);
  } else {
    size_t l, m;
    bool mid_eq;
    std::tie(l, m, mid_eq) = p_split3(In, Out, f);
    par_do3(
        [&]() {
          p_quicksort_(Out.slice(0, l), In.slice(0, l), f, !inplace, cut_size);
        },
        [&]() {
          auto copy_in = [&](size_t i) { In[i] = Out[i]; };
          if (!mid_eq)
            p_quicksort_(Out.slice(l, m), In.slice(l, m), f, !inplace,
                         cut_size);
          else if (inplace)
            parallel_for(l, m, copy_in, 2000);
        },
        [&]() {
          p_quicksort_(Out.slice(m, n), In.slice(m, n), f, !inplace, cut_size);
        });
  }
}

template <class SeqA, class F>
sequence<typename SeqA::value_type> p_quicksort(SeqA const& In, const F& f) {
  using T = typename SeqA::value_type;
  sequence<T> Out(In.size());
  p_quicksort_(In.slice(), Out.slice(), f);
  return Out;
}

template <class T, class F>
void p_quicksort_inplace(range<T*> In, const F& f) {
  sequence<T> Tmp = sequence<T>::no_init(In.size());
  p_quicksort_(In, Tmp.slice(), f, true);
}

}  // namespace pbbs
