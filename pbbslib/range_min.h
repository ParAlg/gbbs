// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2019 Guy Blelloch and the PBBS team
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

#include "parallel.h"
#include "sequence.h"
#include "utilities.h"

namespace pbbs {

// builds a static range minima query structure:
//  range_min(a, less) builds the structure on the sequence a
//    based on the comparison less
//    also takes an optional block_size argument (default = 32)
//  query(i,j) finds the minimum in the range from i to j inclusive of both
//    and returns its index
// Assuming less takes constant time:
//   Build takes O(n log n / block_size) time
//   Query takes O(block_size) time
template <class Seq, class Compare, class Uint = uint>
class range_min {
 public:
  range_min(Seq &_a, Compare _less, long work_block_size = 32)
      : a(_a), less(_less), n(a.size()), block_size(work_block_size) {
    m = 1 + (n - 1) / block_size;
    precomputeQueries();
  }

  Uint query(Uint i, Uint j) {
    // same or adjacent blocks
    if (j - i < block_size) {
      Uint r = i;
      for (long k = i + 1; k <= j; k++) r = min_index(r, k);
      return r;
    }
    long block_i = i / block_size;
    long block_j = j / block_size;
    Uint minl = i;

    // min suffix of first block
    for (long k = minl + 1; k < (block_i + 1) * block_size; k++)
      minl = min_index(minl, k);

    // min prefix of last block
    long minr = block_j * block_size;
    for (long k = minr + 1; k <= j; k++) minr = min_index(minr, k);

    // if adjacent, then done
    if (block_j == block_i + 1) return min_index(minl, minr);

    Uint outOfBlockMin;
    block_i++;
    block_j--;
    if (block_j == block_i)
      outOfBlockMin = table[0][block_i];
    else if (block_j == block_i + 1)
      outOfBlockMin = table[1][block_i];
    else {
      long k = pbbs::log2_up(block_j - block_i + 1) - 1;
      long p = 1 << k;  // 2^k
      outOfBlockMin = min_index(table[k][block_i], table[k][block_j + 1 - p]);
    }
    return min_index(minl, min_index(outOfBlockMin, minr));
  }

 private:
  Seq &a;
  sequence<sequence<Uint>> table;
  Compare less;
  long n, m, depth, block_size;

  Uint min_index(Uint i, Uint j) { return less(a[j], a[i]) ? j : i; }

  void precomputeQueries() {
    depth = log2_up(m + 1);
    table = sequence<sequence<Uint>>(
        depth, [&](size_t i) { return sequence<Uint>(m); });

    // minimums within each block
    sliced_for(n, block_size, [&](size_t i, size_t start, size_t end) {
      long k = start;
      for (size_t j = start + 1; j < end; j++) k = min_index(j, k);
      table[0][i] = k;
    });

    // minimum across layers
    long dist = 1;
    for (long j = 1; j < depth; j++) {
      parallel_for(0, m - dist, [&](size_t i) {
        table[j][i] = min_index(table[j - 1][i], table[j - 1][i + dist]);
      });
      parallel_for(m - dist, m,
                   [&](size_t i) { table[j][i] = table[j - 1][i]; });
      dist *= 2;
    }
  }
};

template <class Seq, class Compare, class Uint = uint>
range_min<Seq, Compare, Uint> make_range_min(Seq &a, Compare less,
                                             long block_size = 32) {
  return range_min<Seq, Compare, Uint>(a, less, block_size);
}
}
