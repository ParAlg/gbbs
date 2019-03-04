// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2016 Guy Blelloch, Daniel Ferizovic, and the PBBS team
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

// A concurrent allocator for any fixed type T
// Keeps a local pool per processor
// Grabs list_size elements from a global pool if empty, and
// Returns list_size elements to the global pool when local pool=2*list_size
// Keeps track of number of allocated elements.
// Probably more efficient than a general purpose allocator

#pragma once

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "block_allocator.h"

template <typename T>
struct allocator {

  static const size_t log_max_size = 10;
  static const size_t log_min_size = 4;
  static const size_t header_size = 8;
  struct block_allocator allocators[log_max_size-log_min_size+1];

  allocator() {
    for (size_t i = 0; i <= log_max_size; i++) {
      allocators[i] = block_allocator(1 << (i+log_min_size));
    }
  }

  void* alloc(size_t n) {
    size_t np = n + header_size;  // add header
    
    size_t log_size = pbbs::log2_up((size_t) np);
    if (log_size > log_max_size) abort();
    int bucket = log_size - log_min_size;
    uchar* ptr = allocators[bucket].alloc();
    *((size_t*) ptr) = bucket;
    return ptr+header_size;
  }

  void free(void* ptr) {
    uchar* head = ((uchar*) ptr) - header_size;
    int bucket = *((size_t*) head);
    allocators[bucket].free(head);
  }
};
