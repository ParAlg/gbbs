// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
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

#include "sequence.h"

namespace pbbs {

// The following supports both "union" that is only safe sequentially
// and "link" that is safe in parallel.  Find is always safe in parallel.
// See:  "Internally deterministic parallel algorithms can be fast"
// Blelloch, Fineman, Gibbons, and Shun
// for a discussion of link/find.
template <class vertexId>
struct unionFind {
  pbbs::sequence<vertexId> parents;

  bool is_root(vertexId u) { return parents[u] < 0; }

  // initialize n elements all as roots
  unionFind(size_t n) { parents = pbbs::sequence<vertexId>(n, -1); }

  vertexId find(vertexId i) {
    if (is_root(i)) return i;
    vertexId p = parents[i];
    if (is_root(p)) return p;

    // find root, shortcutting along the way
    do {
      vertexId gp = parents[p];
      parents[i] = gp;
      i = p;
      p = gp;
    } while (!is_root(p));
    return p;
  }

  // If using "union" then "parents" are used both as
  // parent pointer and for rank (a bit of a hack).
  // When a vertex is a root (negative) then the magnitude
  // of the negative number is its rank.
  // Otherwise it is the vertexId of its parent.
  // cannot be called union since reserved in C
  void union_roots(vertexId u, vertexId v) {
    if (parents[v] < parents[u]) std::swap(u, v);
    // now u has higher rank (higher negative number)
    parents[u] += parents[v];  // update rank of root
    parents[v] = u;            // update parent of other tree
  }

  // Version of union that is safe for parallelism
  // when no cycles are created (e.g. only link from larger
  // to smaller vertexId).
  // Does not use ranks.
  void link(vertexId u, vertexId v) { parents[u] = v; }

  // returns true if successful
  bool tryLink(vertexId u, vertexId v) {
    return (parents[u] == -1 &&
            pbbs::atomic_compare_and_swap(&parents[u], -1, v));
  }
};

}  // namespace pbbs
