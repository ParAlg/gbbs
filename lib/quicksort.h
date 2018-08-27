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
#include "seq.h"
#include "utilities.h"

namespace pbbs {

template <class E, class BinPred>
void insertion_sort(E* A, size_t n, const BinPred& f) {
  for (size_t i = 0; i < n; i++) {
    E v = A[i];
    E* B = A + i;
    while (--B >= A && f(v, *B)) *(B + 1) = *B;
    *(B + 1) = v;
  }
}

// sorts 5 elements and puts them at the front
template <class E, class BinPred>
void sort5(E* A, size_t n, const BinPred& f) {
  // E B[5];
  size_t size = 5;
  size_t m = n / (size + 1);
  for (size_t l = 0; l < size; l++)
    // B[l] = A[m*(l+1)];
    std::swap(A[l], A[m * (l + 1)]);
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
    p2 = p1;              // if few elements greater than p2, then set to p1
  E* L = A;               // below L are less than p1
  E* R = A + n - 1;       // above R are greater than p2
  while (f(*L, p1)) L++;  // set up initial invariant
  while (f(p2, *R)) R--;  // set up initial invariant
  E* M = L;               // between L and M are between p1 and p2 inclusive
  while (M <= R) {
    E mid = *M;
    if (f(mid, p1)) {
      std::swap(*M, *L);
      L++;
    } else if (f(p2, mid)) {
      if (f(*R, p1)) {
        *M = *L;
        *L = *R;
        L++;
      } else {
        *M = *R;
      }
      *R = mid;
      R--;
      while (f(p2, *R)) R--;
    }
    M++;
  }
  return std::make_tuple(L, M, !f(p1, p2));
}

template <class E, class BinPred>
void quicksort_serial(E* A, size_t n, const BinPred& f) {
  while (n > 24) {
    std::tuple<E*, E*, bool> X = split3(A, n, f);
    E* L = std::get<0>(X);  // C++ is so incredibly ugly
    E* M = std::get<1>(X);
    bool mid_eq = std::get<2>(X);
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
    std::tuple<E*, E*, bool> X = split3(A, n, f);
    E* L = std::get<0>(X);
    E* M = std::get<1>(X);
    bool mid_eq = std::get<2>(X);
    cilk_spawn quicksort(A, L - A, f);
    if (!mid_eq) cilk_spawn quicksort(L, M - L, f);
    quicksort(M, A + n - M, f);
    cilk_sync;
  }
}

//// Fully Parallel version below here

template <class SeqA, class BinPred>
std::tuple<size_t, size_t, bool> p_split3(SeqA A, SeqA B, const BinPred& f) {
  using E = typename SeqA::T;
  size_t n = A.size();
  sort5(A.get_array(), n, f);
  E p1 = A[1];
  E p2 = A[3];
  //  cout << p1 << ", " << p2 << ", " << endl;
  if (!f(A[0], A[1])) p1 = p2;  // if few elements less than p1, then set to p2
  if (!f(A[3], A[4]))
    p2 = p1;  // if few elements greater than p2, then set to p1
  auto flag = [&](size_t i) { return (A[i] < p1) ? 0 : (A[i] > p2) ? 2 : 1; };
  auto r = split_three(A, B, make_sequence<unsigned char>(n, flag));
  return std::make_tuple(r.first, r.first + r.second, !f(p1, p2));
}

template <class SeqA, class F>
void p_quicksort(SeqA In, SeqA Out, const F& f, bool swap = 0) {
  size_t n = In.size();
  if (n < (1 << 20)) {
    quicksort(In.get_array(), n, f);
    if (!swap)
      for (size_t i = 0; i < n; i++) Out[i] = In[i];
  } else {
    std::tuple<size_t, size_t, bool> X = p_split3(In, Out, f);
    size_t l = std::get<0>(X);
    size_t m = std::get<1>(X);
    bool mid_eq = std::get<2>(X);
    par_do3(
        true, [&]() { p_quicksort(Out.slice(0, l), In.slice(0, l), f, !swap); },
        [&]() {
          if (!mid_eq)
            p_quicksort(Out.slice(l, m), In.slice(l, m), f, !swap);
          else if (swap)
            parallel_for_bc(i, l, m, ((m - l) > pbbs::kSequentialForThreshold),
                            { In[i] = Out[i]; });
        },
        [&]() { p_quicksort(Out.slice(m, n), In.slice(m, n), f, !swap); });
  }
}

}  // namespace pbbs
