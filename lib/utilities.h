// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2016 Guy Blelloch and the PBBS team
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

#include <ctype.h>
#include <stdlib.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>
#include <type_traits>

#include "macros.h"

#if defined(CILK)
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#define parallel_for cilk_for
#define parallel_for_1 cilk_for
#define parallel_for_16 _Pragma("cilk_grainsize = 16") cilk_for
#define parallel_for_256 _Pragma("cilk_grainsize = 256") cilk_for
size_t nworkers() { return __cilkrts_get_nworkers(); }
static int getWorkers() { return __cilkrts_get_nworkers(); }
static void setWorkers(int n) {
  __cilkrts_end_cilk();
  std::stringstream ss;
  ss << n;
  if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
    std::cerr << "failed to set worker count!" << std::endl;
    std::abort();
  }
}
static inline int get_worker_num() { return __cilkrts_get_worker_number(); }

template <typename Lf, typename Rf>
static void par_do(bool do_parallel, Lf left, Rf right) {
  if (do_parallel) {
    cilk_spawn right();
    left();
    cilk_sync;
  } else {
    left();
    right();
  }
}

template <typename Lf, typename Mf, typename Rf>
static void par_do3(bool do_parallel, Lf left, Mf mid, Rf right) {
  if (do_parallel) {
    cilk_spawn mid();
    cilk_spawn right();
    left();
    cilk_sync;
  } else {
    left();
    mid();
    right();
  }
}

template <typename F>
static void par_for(size_t start, size_t end, size_t granularity, F f) {
  if ((end - start) <= granularity)
    for (size_t i = start; i < end; i++) f(i);
  else {
    size_t mid = (end + start) / 2;
    cilk_spawn par_for(start, mid, granularity, f);
    par_for(mid, end, granularity, f);
    cilk_sync;
  }
}
#else
#define cilk_spawn
#define cilk_sync
#define parallel_for for
#define parallel_for_1 for
#define parallel_for_256 for
#define cilk_for for
static int getWorkers() { return 1; }
static void setWorkers(int n) {}

template <typename Lf, typename Rf>
static void par_do(bool do_parallel, Lf left, Rf right) {
  left();
  right();
}

template <typename Lf, typename Mf, typename Rf>
static void par_do3(bool do_parallel, Lf left, Mf mid, Rf right) {
  left();
  mid();
  right();
}

template <typename F>
static void par_for(size_t start, size_t end, size_t granularity, F f) {
  for (size_t i = start; i < end; i++) f(i);
}
#endif

#include <malloc.h>
struct malloc_init {
  static int i;
  static int j;
};
int malloc_init::i = mallopt(M_MMAP_MAX, 0);
int malloc_init::j = mallopt(M_TRIM_THRESHOLD, -1);

#include <sys/resource.h>
void increase_stack_size() {
  size_t ssize = 1024 * 1024 * 1024;
  // ssize *= 2;
  const rlim_t kStackSize = ssize;
  struct rlimit rl;
  int result;

  result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0) {
    if (rl.rlim_cur < kStackSize) {
      rl.rlim_cur = kStackSize;
      result = setrlimit(RLIMIT_STACK, &rl);
      if (result != 0) {
        fprintf(stderr, "setrlimit returned result = %d\n", result);
      }
    }
  }
}

#define granular_for(_i, _start, _end, _cond, _body) { \
  if (_cond) { \
    {parallel_for(size_t _i=_start; _i < _end; _i++) { \
      _body \
    }} \
  } else { \
    {for (size_t _i=_start; _i < _end; _i++) { \
      _body \
    }} \
  } \
  }

namespace pbbs {

struct empty {};

typedef uint32_t flags;
const flags no_flag = 0;
const flags fl_sequential = 1;
const flags fl_debug = 2;
const flags fl_time = 4;

template <typename T>
inline void assign_uninitialized(T& a, const T& b) {
  new (static_cast<void*>(std::addressof(a))) T(b);
}

template <typename T>
inline void move_uninitialized(T& a, const T& b) {
  new (static_cast<void*>(std::addressof(a))) T(std::move(b));
}

// a 32-bit hash function
uint32_t hash32(uint32_t a) {
  a = (a + 0x7ed55d16) + (a << 12);
  a = (a ^ 0xc761c23c) ^ (a >> 19);
  a = (a + 0x165667b1) + (a << 5);
  a = (a + 0xd3a2646c) ^ (a << 9);
  a = (a + 0xfd7046c5) + (a << 3);
  a = (a ^ 0xb55a4f09) ^ (a >> 16);
  return a;
}

// from numerical recipes
uint64_t hash64(uint64_t u) {
  uint64_t v = u * 3935559000370003845 + 2691343689449507681;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >> 4;
  v *= 4768777513237032717;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v << 5;
  return v;
}

// Does not initialize the array
template <typename E>
E* new_array_no_init(size_t n, bool touch_pages = false) {
  // pads in case user wants to allign with cache lines
  size_t line_size = 64;
  size_t bytes = ((n * sizeof(E)) / line_size + 1) * line_size;
  E* r = (E*)aligned_alloc(line_size, bytes);
  if (r == NULL) {
    fprintf(stderr, "Cannot allocate space");
    exit(1);
  }
  // a hack to make sure tlb is full for huge pages
  if (touch_pages)
    parallel_for(size_t i = 0; i < bytes; i = i + (1 << 21))((bool*)r)[i] = 0;
  return r;
}

// Initializes in parallel
template <typename E>
E* new_array(size_t n) {
  E* r = new_array_no_init<E>(n);
  if (!std::is_trivially_default_constructible<E>::value) {
    if (n > 2048)
      parallel_for(size_t i = 0; i < n; i++) new ((void*)(r + i)) E;
    else
      for (size_t i = 0; i < n; i++) new ((void*)(r + i)) E;
  }
  return r;
}

// Destructs in parallel
template <typename E>
void delete_array(E* A, size_t n) {
  // C++14 -- suppored by gnu C++11
  if (!std::is_trivially_destructible<E>::value) {
    if (n > 2048)
      parallel_for(size_t i = 0; i < n; i++) A[i].~E();
    else
      for (size_t i = 0; i < n; i++) A[i].~E();
  }
  free(A);
}

template <typename ET>
inline bool CAS_GCC(ET* ptr, ET oldv, ET newv) {
  return __sync_bool_compare_and_swap(ptr, oldv, newv);
}

template <class ET>
inline bool CAS(ET *ptr, ET oldv, ET newv) {
  if (sizeof(ET) == 1) {
    return __sync_bool_compare_and_swap((bool*)ptr, *((bool*)&oldv), *((bool*)&newv));
  } else if (sizeof(ET) == 4) {
    return __sync_bool_compare_and_swap((int*)ptr, *((int*)&oldv), *((int*)&newv));
  } else if (sizeof(ET) == 8) {
    return __sync_bool_compare_and_swap((long*)ptr, *((long*)&oldv), *((long*)&newv));
  }
  else {
    std::cout << "CAS bad length : " << sizeof(ET) << std::endl;
    abort();
  }
}

template <class ET>
inline bool CAS128(ET* a, ET b, ET c) {
  return __sync_bool_compare_and_swap_16((__int128*)a, *((__int128*)&b),
                                         *((__int128*)&c));
}

inline long xaddl(long* variable, long value) {
  asm volatile("lock; xaddl %%eax, %2;"
               : "=a"(value)                 // Output
               : "a"(value), "m"(*variable)  // Input
               : "memory");
  return value;
}

inline int xaddi(int* variable, int value) {
  asm volatile("lock; xadd %%eax, %2;"
               : "=a"(value)                 // Output
               : "a"(value), "m"(*variable)  // Input
               : "memory");
  return value;
}

// The conditional should be removed by the compiler
// this should work with pointer types, or pairs of integers
template <class ET>
inline ET xadd(ET* variable, ET value) {
  if (sizeof(ET) == 8) {
    return xaddl((long*)variable, (long)value);
  } else if (sizeof(ET) == 4) {
    return xaddi((int*)variable, (int)value);
  } else {
    std::cout << "xadd bad length" << std::endl;
    abort();
  }
}

template <typename E, typename EV>
inline E fetch_and_add(E* a, EV b) {
  volatile E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!CAS_GCC(a, oldV, newV));
  return oldV;
}

template <typename E, typename EV>
inline void write_add(E* a, EV b) {
  volatile E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!CAS_GCC(a, oldV, newV));
}

template <typename E, typename EV>
inline void write_xor(E* a, EV b) {
  volatile E newV, oldV;
  do {
    oldV = *a;
    newV = oldV ^ b;
  } while (!CAS_GCC(a, oldV, newV));
}

template <typename ET, typename F>
inline bool write_min(ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = *a;
  while (less(b, c) && !(r = CAS_GCC(a, c, b)));
  return r;
}

// returns the log base 2 rounded up (works on ints or longs or unsigned
// versions)
template <class T>
static int log2_up(T i) {
  int a = 0;
  T b = i - 1;
  while (b > 0) {
    b = b >> 1;
    a++;
  }
  return a;
}
};  // namespace pbbs
