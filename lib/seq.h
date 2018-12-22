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

#include "utilities.h"

template <typename E>
struct sequence {
 public:
  using T = E;

  sequence() : allocated(false) {}

  // copy constructor
  sequence(const sequence& a) : s(a.s), e(a.e), allocated(false) {}

  // move constructor
  sequence(sequence&& b) : s(b.s), e(b.e), allocated(b.allocated) {
    b.s = b.e = NULL;
    b.allocated = false;
  }

  // copy assignment
  sequence& operator=(const sequence& b) {
    if (this != &b) {
      clear();
      s = b.s;
      e = b.e;
      allocated = false;
    }
    return *this;
  }

  // move assignment
  sequence& operator=(sequence&& b) {
    if (this != &b) {
      clear();
      s = b.s;
      e = b.e;
      allocated = b.allocated;
      b.s = b.e = NULL;
      b.allocated = false;
    }
    return *this;
  }

  sequence(const size_t n) : s(pbbs::new_array<E>(n)), allocated(true) {
    e = s + n;
  };

  sequence(const size_t n, T v)
      : s(pbbs::new_array_no_init<E>(n, 1)), allocated(true) {
    e = s + n;
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                    { new ((void*)(s + i)) T(v); });
  };

  template <typename Func>
  sequence(const size_t n, Func fun)
      : s(pbbs::new_array_no_init<E>(n)), allocated(true) {
    e = s + n;
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      T x = fun(i);
      new ((void*)(s + i)) T(x);
    });
  }

  sequence(E* s, const size_t n, bool allocated = false)
      : s(s), e(s + n), allocated(allocated){};

  sequence(E* s, E* e, bool allocated = false)
      : s(s), e(e), allocated(allocated){};

  ~sequence() { clear(); }

  template <typename X, typename F>
  static sequence<X> tabulate(size_t n, F f) {
    X* r = pbbs::new_array_no_init<X>(n);
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                    { new ((void*)(r + i)) X(f(i)); });
    sequence<X> y(r, n);
    y.allocated = true;
    return y;
  }

  sequence copy(sequence& a) {
    return tabulate(e - s, [&](size_t i) { return a[i]; });
  }

  E& operator[](const size_t i) const { return s[i]; }
  E& operator()(const size_t i) const { return s[i]; }

  sequence slice(size_t ss, size_t ee) { return sequence(s + ss, s + ee); }

  void update(const size_t i, T& v) { s[i] = v; }

  size_t size() { return e - s; }

  sequence as_sequence() { return sequence(s, e); }

  // Changed to release memory: careful.
  T* get_array() {
    allocated = false;
    return s;
  }
  T* start() { return s; }
  T* end() { return e; }

  sequence<E> copy() {
    if (allocated) {
      size_t n = e - s;
      auto A = pbbs::new_array<E>(n);
      parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                      { A[i] = s[i]; });
      return sequence(A, n, true);
    } else {
      return sequence(s, e);
    }
  }

  void del() { clear(); }
  void clear() {
    if (allocated) {
      allocated = false;
      pbbs::delete_array<E>(s, e - s);
    }
    s = e = NULL;
  }
  E* s;  // = NULL;
  E* e;  // = NULL;
  bool allocated = false;
};

template <typename E, typename F>
struct func_sequence {
  using T = E;
  func_sequence(size_t n, F& _f) : f(&_f), s(0), e(n){};
  func_sequence(size_t s, size_t e, F& _f) : f(&_f), s(s), e(e){};
  T operator()(const size_t i) { return (*f)(i + s); }
  T operator[](const size_t i) { return (*f)(i + s); }
  func_sequence<T, F> slice(size_t ss, size_t ee) {
    return func_sequence<T, F>(s + ss, s + ee, *f);
  }
  size_t size() { return e - s; }
  sequence<T> as_sequence() {
    return sequence<T>::template tabulate<T>(
        e - s, [&](size_t i) { return (*f)(i + s); });
  }

 private:
  F* f;
  size_t s, e;
};

// used so second template argument can be inferred
template <class E, class F>
inline func_sequence<E, F> make_sequence(size_t n, F f) {
  return func_sequence<E, F>(n, f);
}

template <class E>
inline sequence<E> make_sequence(E* A, size_t n) {
  return sequence<E>(A, n);
}
