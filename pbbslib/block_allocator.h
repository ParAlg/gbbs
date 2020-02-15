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

#include <stdlib.h>

#include "concurrent_stack.h"

struct block_allocator {
 private:

  static const size_t default_alloc_size = 1000000;
  static const size_t default_list_size = 1 << 16;
  static const size_t pad_size = 256;

  struct block {
    block* next;
  };

  using block_p = block*;

  struct alignas(64) thread_list {
    size_t sz;
    block_p head;
    block_p mid;
    char cache_line[pad_size];
    thread_list() : sz(0), head(NULL) {};
  };


  block_p initialize_list(block_p);
  block_p get_list();
  concurrent_stack<unsigned char*> pool_roots;
  concurrent_stack<block_p> global_stack;
  thread_list* local_lists;

  size_t list_length;
  size_t max_blocks;
  size_t block_size_;
  std::atomic<size_t> blocks_allocated;
  unsigned char* allocate_blocks(size_t num_blocks);

 public:
  static int thread_count;
  void* alloc();
  void free(void*);
  void reserve(size_t n = default_alloc_size);
  size_t block_size () {return block_size_;}
  size_t num_allocated_blocks() {return blocks_allocated;}
  size_t num_used_blocks();

  ~block_allocator();
  block_allocator(size_t block_size,
		 size_t blocks_count = default_alloc_size,
		 size_t list_length_ = default_list_size,
		 size_t max_blocks_ = 0);

};
