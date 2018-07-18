// A simple dynamic array implementation.
//
// Note: Have only tested this for simple datatypes, in particular, it uses newA
// (i.e. malloc) to resize, which will cause cryptic segfaults if E requires a
// constructor call to properly initialize memory. see pbbs::new_array(..)
#pragma once

#include "utilities.h"

#define MIN_BKT_SIZE 2048

template <class E>
struct dyn_arr {
  E* A;
  size_t size;
  size_t capacity;
  bool alloc;

  dyn_arr() : A(NULL), size(0), capacity(0), alloc(false) {}
  dyn_arr(size_t s) : size(0), capacity(s), alloc(true) { A = newA(E, s); }
  dyn_arr(E* _A, long _size, long _capacity, bool _alloc) :
    A(_A), size(_size), capacity(_capacity), alloc(_alloc) {}

  void del() {
    if (alloc) {
      free(A);
      alloc = false;
    }
  }

  void clear() { size = 0; }

  inline void resize(size_t n) {
    if (n + size > capacity) {
      size_t new_capacity = std::max(2*(n + size), (size_t)MIN_BKT_SIZE);
      E* nA = newA(E, new_capacity);
      granular_for(i, 0, size, 2000, nA[i] = A[i];);
      if (alloc) {
        free(A);
      }
      A = nA;
      capacity = new_capacity;
      alloc = true;
    }
  }

  inline void insert(E val, size_t pos) {
    A[size + pos] = val;
  }

  inline void push_back(E val) {
    A[size] = val;
    size++;
  }

  template <class F>
  void map(F f) {
    granular_for (i, 0, size, (size > 2000) , {
      f(A[i]);
    });
  }

  template <class F>
  inline void copyIn(F f, size_t n) {
    resize(n);
    granular_for(i, 0, n, 2000, A[size + i] = f[i];);
    size += n;
  }

  template <class F>
  inline void copyInF(F f, size_t n) {
    resize(n);
    granular_for(i, 0, n, 2000, A[size + i] = f(i););
    size += n;
  }

};


