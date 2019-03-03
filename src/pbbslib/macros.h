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

#include <cmath>
#include <iostream>
#include <stddef.h>
#include "lib/parallel.h"

// scan/filter macros; used by sequence implementations
#define _SCAN_LOG_BSIZE 10
#define _SCAN_BSIZE (1 << _SCAN_LOG_BSIZE)
#define _F_BSIZE (2 * _SCAN_BSIZE)


template <typename F>
static void par_for(size_t start, size_t end, size_t granularity, F f, bool parallel=true) {
  if (!parallel) {
    for (size_t i=start; i<end; i++) {
      f(i);
    }
  } else {
    parallel_for(start, end, f, granularity);
  }
}

template <typename F>
static void par_for(size_t start, size_t end, F f, bool parallel=true) {
  size_t n = end - start;
  size_t granularity = (n > 100) ? ceil(sqrt(n)) : 100;
  par_for<F>(start, end, granularity, f, parallel);
}

namespace pbbs {

  constexpr const size_t kSequentialForThreshold = 4000;

  template <class ET>
  inline bool CAS(ET* ptr, ET oldv, ET newv) {
    if (sizeof(ET) == 1) {
      return __sync_bool_compare_and_swap((bool*)ptr, *((bool*)&oldv),
                                          *((bool*)&newv));
    } else if (sizeof(ET) == 4) {
      return __sync_bool_compare_and_swap((int*)ptr, *((int*)&oldv),
                                          *((int*)&newv));
    } else if (sizeof(ET) == 8) {
      return __sync_bool_compare_and_swap((long*)ptr, *((long*)&oldv),
                                          *((long*)&newv));
    } else {
      std::cout << "CAS bad length : " << sizeof(ET) << "\n";
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
      std::cout << "xadd bad length"
                << "\n";
      abort();
    }
  }
}  // namespace pbbs
