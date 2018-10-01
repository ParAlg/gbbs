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

// An index map is an alternative to iterators that can be better
// in certain applications, especially parallel code.
//
// Index maps have an associated lenght n (A.size()),
// and supply the functions:
//     A[i] that  "reads" the i'th index, and
//     A.update(i,v) that updates the i'th index with v.
// Both are only valid for i = [0,n).
// Input index maps only supply the first, output index maps only the
// second, and input/output index maps supply both.

// Index maps are useful when not just reading or writing contiguously in
// an arrary.    They can save time and memory by avoiding copying.
// Some applications include
//
// accessing every k'th element of an array A:
//     make_in_imap<T>(n, [&] (size_t i) {return A[i*k];})
// accessing fields of a structure:
//     make_in_imap<T>(n, [&] (size_t i) {return A[i].my_field;})
// applying a function to each element of an array:
//     make_in_imap<T>(n, [&] (size_t i) {return f(A[i]);})
// defining elements to be some function of i, e.g. its square:
//     make_in_imap<T>(n, [] (size_t i) {return i*i;})
// an arbitrary function can be applied on writing:
//     make_out_imap<T>(n, [&] (size_t i, T v) {A[i] = f(v);})
// or every other element can be written:
//     make_out_imap<T>(n, [&] (size_t i, T v) {A[i*n] = f(v);})
//
// Random access iterators can be converted into an index map:
//     make_iter_map(start, end)
// This creates an input/output map
//
// And the simplest usage is converting an array into an index map:
//     make_array_map(A, n)
// This also creates an input/output map
#pragma once

#include "utilities.h"

template <typename E, typename F>
struct in_imap {
  using T = E;
  F f;
  size_t s, e;
  in_imap(size_t n, F _f) : f(_f), s(0), e(n){};
  in_imap(size_t s, size_t e, F _f) : f(_f), s(s), e(e){};
  T operator[](const size_t i) { return f(i + s); }
  T operator()(const size_t i) { return f(i + s); }
  in_imap<T, F> cut(size_t ss, size_t ee) {
    return in_imap<T, F>(s + ss, s + ee, f);
  }
  in_imap<T, F> slice(size_t ss, size_t ee) {
    return in_imap<T, F>(s + ss, s + ee, f);
  }
  size_t size() { return e - s; }
};

// used so second template argument can be inferred
template <class E, class F>
inline in_imap<E, F> make_in_imap(size_t n, F f) {
  return in_imap<E, F>(n, f);
}

template <typename E, typename F>
struct out_imap {
  using T = E;
  F f;
  size_t s, e;
  out_imap(size_t n, F f) : f(f), s(0), e(n){};
  out_imap(size_t s, size_t e, F f) : f(f), s(s), e(e){};
  out_imap<T, F> cut(size_t ss, size_t ee) {
    return out_imap<T, F>(s + ss, s + ee, f);
  }
  out_imap<T, F> slice(size_t ss, size_t ee) {
    return out_imap<T, F>(s + ss, s + ee, f);
  }
  void update(size_t i, const T& val) { f(i + s, val); }
  size_t size() { return e - s; }
};

// used so second template argument can be inferred
template <class E, class F>
inline out_imap<E, F> make_out_imap(size_t n, F f) {
  return out_imap<E, F>(n, f);
}

template <typename Iterator>
struct iter_imap {
  using T = typename std::iterator_traits<Iterator>::value_type;
  Iterator s;
  const Iterator e;
  iter_imap(const iter_imap& b) : s(b.s), e(b.e) {}
  iter_imap() {}
  iter_imap(Iterator s, Iterator e) : s(s), e(e){};
  T& operator[](const size_t i) const { return s[i]; }
  iter_imap<Iterator> cut(size_t ss, size_t ee) {
    return iter_imap<Iterator>(s + ss, s + ee);
  }
  iter_imap<Iterator> slice(size_t ss, size_t ee) {
    return iter_imap<Iterator>(s + ss, s + ee);
  }
  void update(size_t i, const T& val) { s[i] = val; }
  const size_t size() { return e - s; }
};

// used so template argument can be inferred
template <class Iterator>
inline iter_imap<Iterator> make_iter_imap(Iterator s, Iterator e) {
  return iter_imap<Iterator>(s, e);
}

template <typename E>
struct array_imap {
 public:
  using T = E;

  array_imap() : allocated(false) {}
  // copy constructor
  array_imap(const array_imap& b) : s(b.s), e(b.e), allocated(false) {}
  array_imap(array_imap&& b) : s(b.s), e(b.e), allocated(b.allocated) {
    b.allocated = false;
  }
  array_imap(E* s, size_t n, bool _allocated = false)
      : s(s), e(s + n), allocated(_allocated){};
  array_imap(const size_t n)
      : s(pbbs::new_array<E>(n)), e(s + n), allocated(true){};
  template <class F>
  array_imap(const size_t n, F f)
      : s(pbbs::new_array_no_init<E>(n)), e(s + n), allocated(true) {
    parallel_for_bc(i, 0, n, (n > 2000), { s[i] = f(i); });
  };
  array_imap(const size_t n, const E e)
      : s(pbbs::new_array_no_init<E>(n)), e(s + n), allocated(true) {
    parallel_for_bc(i, 0, n, (n > 2000), { s[i] = e; });
  };
  ~array_imap() { free_mem(); }

  // copy assignment
  array_imap& operator=(const array_imap& m) {
    if (this != &m) {
      free_mem();
      s = m.s;
      e = m.e;
      allocated = false;
    }
    return *this;
  }
  // move assignment
  array_imap& operator=(array_imap&& m) {
    if (this != &m) {
      s = m.s;
      e = m.e;
      allocated = m.allocated;
      m.allocated = false;
    }
    return *this;
  }
  void del() {
    if (allocated) {
      free_mem();
    }
  }
  inline E& operator[](const size_t i) const { return s[i]; }
  inline E& operator()(const size_t i) const { return s[i]; }
  array_imap<T> cut(size_t ss, size_t ee) {
    return array_imap<T>(s + ss, ee - ss);
  }
  array_imap<T> slice(size_t ss, size_t ee) {
    return array_imap<T>(s + ss, ee - ss);
  }
  array_imap<T> copy() {
    if (allocated) {
      size_t n = e - s;
      auto A = pbbs::new_array_no_init<E>(n);
      parallel_for_bc(i, 0, n, (n > 200), { A[i] = s[i]; });
      auto ret = array_imap<E>(A, n);
      ret.allocated = true;
      return ret;
    } else {
      return array_imap<E>(s, (e - s));
    }
  }
  void update(size_t i, const E& val) { s[i] = val; }
  const size_t size() { return e - s; }
  T* get_array() {
    allocated = false;
    return s;
  }
  T* start() { return s; }
  T* end() { return e; }

  void free_mem() {
    if (allocated) {
      allocated = false;
      pbbs::delete_array<E>(s, e - s);
    }
  }

  E* s;
  E* e;
  bool allocated;
};

// used so template argument can be inferred
template <class E>
inline array_imap<E> make_array_imap(E* A, size_t n) {
  return array_imap<E>(A, n);
}
