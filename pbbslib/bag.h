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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <atomic>
#include "list_allocator.h"

namespace pbbs {

// a bag of elements of a given type
// supports
//  - creating a singleton bag (constant work)
//  - appending two bags (constant work)
//  - outputting the contents of the bag as a sequence (linear work)
// implemented as a tree, currently only balanced by the append order
template <typename ET>
struct bag {
  using T = ET;
  void* root;
  static const size_t flag = ((size_t)1) << 60;

  struct node {
    size_t size;
    bag left;
    bag right;
  };

  using node_alloc = list_allocator<node>;
  using leaf_alloc = list_allocator<T>;

  static void init() {
    node_alloc::init();
    leaf_alloc::init();
  };

  bag() : root(NULL) {}

  bag(T a) {
    T* x = leaf_alloc::alloc();
    *x = a;
    root = (void*)x;
  }

  static void reserve(size_t n) {
    node_alloc::reserve(n);
    leaf_alloc::reserve(n);
  }

  static bag append(bag a, bag b) {
    if (a.size() == 0) return b;
    if (b.size() == 0) return a;
    node* o = node_alloc::alloc();
    o->left = a;
    o->right = b;
    o->size = a.size() + b.size();
    return bag(o);
  }

  size_t size() {
    if (root == NULL) return 0;
    if (is_node()) return get_node_ptr()->size;
    return 1;
  }

  sequence<T> flatten() {
    sequence<T> out(size());
    flatten_rec(out.begin());
    return out;
  }

 private:
  size_t is_node() { return ((size_t)root) & flag; }
  node* get_node_ptr() { return (node*)(((size_t)root) & ~flag); }
  bag(node* a) { root = (void*)(((size_t)a) | flag); }

  void flatten_rec(T* start) {
    size_t n = size();
    if (n == 1) {
      T* leaf = ((T*)root);
      start[0] = std::move(*leaf);
      leaf_alloc::free(leaf);
    } else if (n > 1) {
      node* x = get_node_ptr();
      size_t nl = x->left.size();
      par_do_if(n > 100, [&]() { (x->left).flatten_rec(start); },
                [&]() { (x->right).flatten_rec(start + nl); });
      node_alloc::free(x);
    }
  }
};

}  // namespace pbbs
