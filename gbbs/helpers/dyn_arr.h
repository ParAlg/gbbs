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

// A simple dynamic array implementation.
//
// Note: Have only tested this for simple datatypes, in particular, it uses newA
// (i.e. malloc) to resize, which will cause cryptic segfaults if E requires a
// constructor call to properly initialize memory. see gbbs::new_array(..)
#pragma once

#include "gbbs/bridge.h"

namespace gbbs {

constexpr size_t kDynArrMinBktSize = 2000;

template <class E>
struct dyn_arr {
  E* A;
  size_t size;
  size_t capacity;
  bool alloc;

  dyn_arr() : A(NULL), size(0), capacity(0), alloc(false) {}
  dyn_arr(size_t s) : size(0), capacity(s), alloc(true) {
    A = gbbs::new_array_no_init<E>(s);
  }
  dyn_arr(E* _A, long _size, long _capacity, bool _alloc)
      : A(_A), size(_size), capacity(_capacity), alloc(_alloc) {}

  void del() {
    if (alloc) {
      gbbs::free_array(A, capacity);
      alloc = false;
    }
  }

  sequence<E> to_seq() {
    assert(A);
    // auto ret = sequence<E>(A, size);
    auto ret = sequence<E>(size);
    gbbs::parallel_for(0, size, [&](size_t i) { ret[i] = A[i]; });
    size = 0;
    A = nullptr;
    return ret;
  }

  void clear() { size = 0; }

  inline void resize(size_t n) {
    if (n + size > capacity) {
      size_t new_capacity = std::max(2 * (n + size), (size_t)kDynArrMinBktSize);
      E* nA = gbbs::new_array_no_init<E>(new_capacity);
      parallel_for(0, size, [&](size_t i) { nA[i] = A[i]; });
      if (alloc) {
        gbbs::free_array(A, capacity);
      }
      A = nA;
      capacity = new_capacity;
      alloc = true;
    }
  }

  inline void insert(E val, size_t pos) { A[size + pos] = val; }

  inline void push_back(E val) {
    A[size] = val;
    size++;
  }

  template <class F>
  void map(F f) {
    parallel_for(0, size, [&](size_t i) { f(A[i]); });
  }

  template <class F>
  inline void copyIn(F& f, size_t n) {
    resize(n);
    parallel_for(0, n, [&](size_t i) { A[size + i] = f[i]; });
    size += n;
  }

  template <class F>
  inline void copyInF(F f, size_t n) {
    resize(n);
    parallel_for(0, n, [&](size_t i) { A[size + i] = f(i); });
    size += n;
  }
};
};  // namespace gbbs
