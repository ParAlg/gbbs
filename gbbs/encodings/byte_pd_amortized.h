// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <type_traits>

#include "gbbs/bridge.h"
#include "gbbs/macros.h"

namespace gbbs {
namespace bytepd_amortized {

inline size_t get_virtual_degree(uintE d, uchar* ngh_arr) {
  if (d > 0) {
    return *((uintE*)ngh_arr);
  }
  return 0;
}

__attribute__((always_inline)) inline uintE eatFirstEdge(uchar*& start,
                                                         const uintE source) {
  uchar fb = *start++;
  uintE edgeRead = (fb & 0x3f);
  if (LAST_BIT_SET(fb)) {
    int shiftAmount = 6;
    while (1) {
      uchar b = *start;
      edgeRead |= ((b & 0x7f) << shiftAmount);
      start++;
      if (LAST_BIT_SET(b))
        shiftAmount += EDGE_SIZE_PER_BYTE;
      else
        break;
    }
  }
  return (fb & 0x40) ? source - edgeRead : source + edgeRead;
}

/*
  Reads any edge of an out-edge list after the first edge.
*/
__attribute__((always_inline)) inline uintE eatEdge(uchar*& start) {
  uintE edgeRead = 0;
  int shiftAmount = 0;

  while (1) {
    uchar b = *start;
    edgeRead += ((b & 0x7f) << shiftAmount);
    start++;
    if (LAST_BIT_SET(b))
      shiftAmount += EDGE_SIZE_PER_BYTE;
    else
      break;
  }
  return edgeRead;
}

// Read default weight (expects gbbs::empty)
template <class W, typename std::enable_if<std::is_same<W, gbbs::empty>::value,
                                           int>::type = 0>
__attribute__((always_inline)) inline W eatWeight(uchar*& start) {
  return (W)gbbs::empty();
}

// Read integer weight
template <class W,
          typename std::enable_if<std::is_same<W, intE>::value, int>::type = 0>
__attribute__((always_inline)) inline W eatWeight(uchar*& start) {
  uchar fb = *start++;
  intE edgeRead = (fb & 0x3f);
  if (LAST_BIT_SET(fb)) {
    int shiftAmount = 6;
    while (1) {
      uchar b = *start;
      edgeRead |= ((b & 0x7f) << shiftAmount);
      start++;
      if (LAST_BIT_SET(b))
        shiftAmount += EDGE_SIZE_PER_BYTE;
      else
        break;
    }
  }
  return (fb & 0x40) ? -edgeRead : edgeRead;
}

// Read unsigned int weight
template <class W, typename std::enable_if<std::is_same<W, uint32_t>::value,
                                           int>::type = 0>
__attribute__((always_inline)) inline W eatWeight(uchar*& start) {
  return eatEdge(start);
}

// Read float weight
template <class W,
          typename std::enable_if<std::is_same<W, float>::value, int>::type = 0>
__attribute__((always_inline)) inline W eatWeight(uchar*& start) {
  float wgh = *((float*)start);
  start += sizeof(float);
  return wgh;
}

// Read double weight
template <class W, typename std::enable_if<std::is_same<W, double>::value,
                                           int>::type = 0>
__attribute__((always_inline)) inline W eatWeight(uchar*& start) {
  double wgh = *((double*)start);
  start += sizeof(double);
  return wgh;
}

template <class W,
          typename std::enable_if<std::is_same<W, intE>::value, int>::type = 0>
__attribute__((always_inline)) inline void print_weight(W& wgh) {
  std::cout << wgh << "\n";
}

template <class W,
          typename std::enable_if<!std::is_same<W, intE>::value, int>::type = 0>
inline void print_weight(W& wgh) {}

/*
  Compresses the first edge, writing target-source and a sign bit.
*/
long compressFirstEdge(uchar* start, long offset, long source, long target);

template <class W,
          typename std::enable_if<!std::is_same<W, intE>::value, int>::type = 0>
inline long compressWeight(uchar* start, long offset, W weight) {
  return offset;
}

template <class W,
          typename std::enable_if<std::is_same<W, intE>::value, int>::type = 0>
inline long compressWeight(uchar* start, long offset, W weight) {
  return compressFirstEdge(start, offset, 0, weight);
}

/*
  Should provide the difference between this edge and the previous edge
*/
long compressEdge(uchar* start, long curOffset, uintE e);

template <class W>
struct iter {
  uchar* base;
  uchar* finger;
  uintE src;
  uintT degree;

  uintE num_blocks;
  uintE cur_chunk;
  uintE cur_chunk_degree;

  std::tuple<uintE, W> last_edge;
  uintE read_in_block;
  uintE read_total;

  iter() {}

  iter(uchar* _base, uintT _degree, uintE _src)
      : base(_base),
        src(_src),
        degree(_degree),
        cur_chunk(0),
        cur_chunk_degree(0) {
    if (degree == 0) return;
    uintE virtual_degree = *((uintE*)base);
    num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)(base + sizeof(uintE));

    finger = base + (num_blocks - 1) * sizeof(uintE) + sizeof(uintE);

    uintE start_offset = *((uintE*)finger);
    uintE end_offset = (0 == (num_blocks - 1))
                           ? degree
                           : (*((uintE*)(base + block_offsets[0])));
    cur_chunk_degree = end_offset - start_offset;
    finger += sizeof(uintE);
    if (start_offset < end_offset) {
      std::get<0>(last_edge) = eatFirstEdge(finger, src);
      std::get<1>(last_edge) = eatWeight<W>(finger);
    } else {
      for (cur_chunk = 1; cur_chunk < num_blocks; cur_chunk++) {
        finger = base + block_offsets[cur_chunk - 1];
        start_offset = *((uintE*)finger);
        end_offset = (cur_chunk == (num_blocks - 1))
                         ? degree
                         : (*((uintE*)(base + block_offsets[cur_chunk])));
        cur_chunk_degree = end_offset - start_offset;
        finger += sizeof(uintE);

        if (start_offset < end_offset) {
          std::get<0>(last_edge) = eatFirstEdge(finger, src);
          std::get<1>(last_edge) = eatWeight<W>(finger);
          break;
        }
      }
    }
    read_total = 1;
    read_in_block = 1;
  }

  __attribute__((always_inline)) inline std::tuple<uintE, W> cur() {
    return last_edge;
  }

  __attribute__((always_inline)) inline std::tuple<uintE, W> next() {
    if (read_in_block == cur_chunk_degree) {
      cur_chunk_degree = 0;
      uintE* block_offsets = (uintE*)(base + sizeof(uintE));
      while (cur_chunk_degree == 0) {
        cur_chunk++;
        finger = base + block_offsets[cur_chunk - 1];
        uintE start_offset = *((uintE*)finger);
        uintE end_offset = (cur_chunk == (num_blocks - 1))
                               ? degree
                               : (*((uintE*)(base + block_offsets[cur_chunk])));
        finger += sizeof(uintE);
        cur_chunk_degree = end_offset - start_offset;
      }

      std::get<0>(last_edge) = eatFirstEdge(finger, src);
      std::get<1>(last_edge) = eatWeight<W>(finger);
      read_in_block = 1;
    } else {
      std::get<0>(last_edge) += eatEdge(finger);
      std::get<1>(last_edge) = eatWeight<W>(finger);
      read_in_block++;
    }
    read_total++;
    return last_edge;
  }

  __attribute__((always_inline)) inline bool has_next() {
    return read_total < degree;
  }
};

template <class W>
struct simple_iter {
  uchar* base;
  uchar* finger;
  uintE src;
  uintT degree;

  uintE cur_chunk;

  std::tuple<uintE, W> last_edge;
  uintE proc;

  simple_iter(uchar* _base, uintT _degree, uintE _src)
      : base(_base), src(_src), degree(_degree), cur_chunk(0) {
    if (degree == 0) return;
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    finger = base + (num_blocks - 1) * sizeof(uintE) + 2 * sizeof(uintE);

    std::get<0>(last_edge) = eatFirstEdge(finger, src);
    std::get<1>(last_edge) = eatWeight<W>(finger);
    proc = 1;
  }

  __attribute__((always_inline)) inline std::tuple<uintE, W> cur() {
    return last_edge;
  }

  __attribute__((always_inline)) inline std::tuple<uintE, W> next() {
    if (proc == PARALLEL_DEGREE) {
      finger += sizeof(uintE);  // skip block start
      std::get<0>(last_edge) = eatFirstEdge(finger, src);
      std::get<1>(last_edge) = eatWeight<W>(finger);
      proc = 1;
      cur_chunk++;
    } else {
      std::get<0>(last_edge) += eatEdge(finger);
      std::get<1>(last_edge) = eatWeight<W>(finger);
      proc++;
    }
    return last_edge;
  }

  __attribute__((always_inline)) inline bool has_next() {
    return (cur_chunk * PARALLEL_DEGREE + proc) < degree;
  }
};

// Decode unweighted edges
template <
    class W, class T,
    typename std::enable_if<std::is_same<W, gbbs::empty>::value, int>::type = 0>
void decode(T& t, uchar* edge_start, const uintE& source, const uintT& degree,
            const bool parallel = true) {
  if (degree > 0) {
    uintE virtual_degree = *((uintE*)edge_start);
    size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));
    uchar* nghs_start = edge_start + (num_blocks - 1) * sizeof(uintE) +
                        sizeof(uintE);  // block offs + virtual_degree

    auto wgh = gbbs::empty();
    {  // do first chunk
      uchar* finger = nghs_start;
      uintE start_offset = *((uintE*)finger);
      uintE end_offset = (0 == (num_blocks - 1))
                             ? degree
                             : (*((uintE*)(edge_start + block_offsets[0])));
      finger += sizeof(uintE);

      if (start_offset < end_offset) {  // at least one edge in this block
        uintE ngh = eatFirstEdge(finger, source);
        if (!t(source, ngh, wgh, start_offset)) return;
        for (size_t edgeID = start_offset + 1; edgeID < end_offset; edgeID++) {
          ngh += eatEdge(finger);
          if (!t(source, ngh, wgh, edgeID)) return;
        }
      }
    }
    if ((num_blocks > 2) && parallel) {
      parallel_for(1, num_blocks, 1, [&](size_t i) {
        uchar* finger =
            (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
        uintE start_offset = *((uintE*)finger);
        uintE end_offset = (i == (num_blocks - 1))
                               ? degree
                               : (*((uintE*)(edge_start + block_offsets[i])));
        finger += sizeof(uintE);

        if (start_offset < end_offset) {  // at least one edge in this block
          uintE ngh = eatFirstEdge(finger, source);
          if (!t(source, ngh, wgh, start_offset)) end_offset = 0;
          for (size_t edgeID = start_offset + 1; edgeID < end_offset;
               edgeID++) {
            ngh += eatEdge(finger);
            if (!t(source, ngh, wgh, edgeID)) break;
          }
        }
      });
    } else {
      for (size_t i = 1; i < num_blocks; i++) {
        uchar* finger =
            (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
        uintE start_offset = *((uintE*)finger);
        uintE end_offset = (i == (num_blocks - 1))
                               ? degree
                               : (*((uintE*)(edge_start + block_offsets[i])));
        finger += sizeof(uintE);

        if (start_offset < end_offset) {  // at least one edge in this block
          uintE ngh = eatFirstEdge(finger, source);
          if (!t(source, ngh, wgh, start_offset)) end_offset = 0;
          for (size_t edgeID = start_offset + 1; edgeID < end_offset;
               edgeID++) {
            ngh += eatEdge(finger);
            if (!t(source, ngh, wgh, edgeID)) break;
          }
        }
      }
    }
  }
}

// Decode weighted edges
template <class W, class T,
          typename std::enable_if<!std::is_same<W, gbbs::empty>::value,
                                  int>::type = 0>
inline void decode(T& t, uchar* edge_start, const uintE& source,
                   const uintT& degree, const bool par = true) {
  if (degree > 0) {
    uintE virtual_degree = *((uintE*)edge_start);
    size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));
    uchar* nghs_start = edge_start + (num_blocks - 1) * sizeof(uintE) +
                        sizeof(uintE);  // block offs + virtual_degree

    // TODO: put back par
    for (size_t i = 0; i < num_blocks; i++) {
      uchar* finger =
          (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
      uintE start_offset = *((uintE*)finger);
      uintE end_offset = (i == (num_blocks - 1))
                             ? degree
                             : (*((uintE*)(edge_start + block_offsets[i])));
      finger += sizeof(uintE);

      if (start_offset < end_offset) {  // at least one edge in this block
        uintE ngh = eatFirstEdge(finger, source);
        W wgh = eatWeight<W>(finger);
        if (!t(source, ngh, wgh, start_offset)) return;
        for (size_t edgeID = start_offset + 1; edgeID < end_offset; edgeID++) {
          ngh += eatEdge(finger);
          wgh = eatWeight<W>(finger);
          if (!t(source, ngh, wgh, edgeID)) return;
        }
      }
    }
  }
}

template <class W, class T>
inline void decode_block(T t, uchar* edge_start, const uintE& source,
                         const uintT& degree, uintE block_num) {
  if (degree > 0) {
    uintE virtual_degree = *((uintE*)edge_start);
    size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));
    uchar* nghs_start =
        edge_start + (num_blocks - 1) * sizeof(uintE) + sizeof(uintE);

    uchar* finger = (block_num == 0)
                        ? nghs_start
                        : (edge_start + block_offsets[block_num - 1]);
    uintE start_offset = *((uintE*)finger);
    uintE end_offset =
        (block_num == (num_blocks - 1))
            ? degree
            : (*((uintE*)(edge_start + block_offsets[block_num])));
    finger += sizeof(uintE);

    if (start_offset < end_offset) {  // at least one edge in this block
      uintE ngh = eatFirstEdge(finger, source);
      W wgh = eatWeight<W>(finger);
      t(ngh, wgh, start_offset);
      for (size_t edgeID = start_offset + 1; edgeID < end_offset; edgeID++) {
        ngh += eatEdge(finger);
        wgh = eatWeight<W>(finger);
        t(ngh, wgh, edgeID);
      }
    }
  }
}

template <class W, class T>
inline void decode_block_cond(T t, uchar* edge_start, const uintE& source,
                              const uintT& degree, uintE block_num) {
  if (degree > 0) {
    uintE virtual_degree = *((uintE*)edge_start);
    size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));
    uchar* nghs_start =
        edge_start + (num_blocks - 1) * sizeof(uintE) + sizeof(uintE);

    uchar* finger = (block_num == 0)
                        ? nghs_start
                        : (edge_start + block_offsets[block_num - 1]);
    uintE start_offset = *((uintE*)finger);
    uintE end_offset =
        (block_num == (num_blocks - 1))
            ? degree
            : (*((uintE*)(edge_start + block_offsets[block_num])));
    finger += sizeof(uintE);

    if (start_offset < end_offset) {  // at least one edge in this block
      uintE ngh = eatFirstEdge(finger, source);
      W wgh = eatWeight<W>(finger);
      if (!t(ngh, wgh, start_offset)) return;
      for (size_t edgeID = start_offset + 1; edgeID < end_offset; edgeID++) {
        ngh += eatEdge(finger);
        wgh = eatWeight<W>(finger);
        if (!t(ngh, wgh, edgeID)) break;
      }
    }
  }
}

// r: E -> E -> E
template <class W, class M, class Monoid>
inline typename Monoid::T map_reduce(uchar* edge_start, const uintE& source,
                                     const uintT& degree, M& m, Monoid& reduce,
                                     const bool par = true) {
  using E = typename Monoid::T;
  if (degree > 0) {
    uintE virtual_degree = *((uintE*)edge_start);
    size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));
    uchar* nghs_start = edge_start + (num_blocks - 1) * sizeof(uintE) +
                        sizeof(uintE);  // block offs + virtual_degree

    E stk[100];
    E* block_outputs;
    parlay::sequence<E> alloc;
    if (num_blocks > 100) {
      alloc = parlay::sequence<E>::uninitialized(num_blocks);
      block_outputs = alloc.begin();
    } else {
      block_outputs = (E*)stk;
    }

    parallel_for(0, num_blocks, 1, [&](size_t i) {
      auto cur = reduce.identity;
      uchar* finger =
          (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
      uintE start_offset = *((uintE*)finger);
      uintE end_offset = (i == (num_blocks - 1))
                             ? degree
                             : (*((uintE*)(edge_start + block_offsets[i])));
      finger += sizeof(uintE);

      if (start_offset < end_offset) {
        uintE ngh = eatFirstEdge(finger, source);
        {  // Eat first edge, which is compressed specially
          W wgh = eatWeight<W>(finger);
          cur = reduce.f(cur, m(source, ngh, wgh));
        }
        for (size_t j = start_offset + 1; j < end_offset; j++) {
          ngh += eatEdge(finger);
          W wgh = eatWeight<W>(finger);
          cur = reduce.f(cur, m(source, ngh, wgh));
        }
      }
      block_outputs[i] = cur;
    });

    auto im = gbbs::make_slice(block_outputs, num_blocks);
    E res = parlay::reduce(im, reduce);
    return res;
  } else {
    return reduce.identity;
  }
}

// Merge:
// (constant, constant) -> merge sequentially
// WLOG, we're intersecting (small, large)
// One option is to decompress fully and do the binary-search merge, but can
// we do it directly on the compressed format?
//
// Split the small side on block boundaries (i.e. the first vtx in a block)
//    - this ensures that one side is always recursing with blocks.
//   The issue is that the binary search might turn up in the middle of some
//   block. I think we can get around this by only checking block-starts. So
//   the algorithm is: pick a pivot from the smaller side. Binary search over
//   the array of block starts, and find the first block whose start element
//   is gt than ours. This wastes at most |block-size| amount of work in the
//   intersection.
template <class W>
inline size_t intersect(uchar* l1, uchar* l2, uintE l1_size, uintE l2_size,
                        uintE l1_src, uintE l2_src) {
  if (l1_size == 0 || l2_size == 0) return 0;
  auto it_1 = simple_iter<W>(l1, l1_size, l1_src);
  auto it_2 = simple_iter<W>(l2, l2_size, l2_src);
  size_t i = 0, j = 0, ct = 0;
  while (i < l1_size && j < l2_size) {
    uintE e1 = std::get<0>(it_1.cur());
    uintE e2 = std::get<0>(it_2.cur());
    if (e1 == e2) {
      i++, j++, it_1.next(), it_2.next(), ct++;
    } else if (e1 < e2) {
      i++, it_1.next();
    } else {
      j++, it_2.next();
    }
  }
  return ct;
}

template <class W, class F>
size_t intersect_f(uchar* l1, uchar* l2, uintE l1_size, uintE l2_size,
                   uintE l1_src, uintE l2_src, const F& f) {
  if (l1_size == 0 || l2_size == 0) return 0;
  auto it_1 = simple_iter<W>(l1, l1_size, l1_src);
  auto it_2 = simple_iter<W>(l2, l2_size, l2_src);
  size_t i = 0, j = 0, ct = 0;
  while (i < l1_size && j < l2_size) {
    uintE e1 = std::get<0>(it_1.cur());
    uintE e2 = std::get<0>(it_2.cur());
    if (e1 == e2) {
      f(l1_src, l2_src, e1);
      i++, j++, it_1.next(), it_2.next(), ct++;
    } else if (e1 < e2) {
      i++, it_1.next();
    } else {
      j++, it_2.next();
    }
  }
  return ct;
}

template <class W>
inline std::tuple<uintE, W> get_ith_neighbor(uchar* edge_start, uintE source,
                                             uintE degree, size_t i) {
  uintE virtual_degree = *((uintE*)edge_start);
  size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
  uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));
  uchar* nghs_start =
      edge_start + (num_blocks - 1) * sizeof(uintE) + sizeof(uintE);
  auto blocks_f = [&](size_t j) {
    uintE end = (j == (num_blocks - 1))
                    ? degree
                    : (*((uintE*)(edge_start + block_offsets[j])));
    return end;
  };
  auto blocks_imap = parlay::delayed_seq<size_t>(num_blocks, blocks_f);
  // This is essentially searching a plus_scan'd, incl arr.
  auto lte = [&](const size_t& l, const size_t& r) { return l <= r; };
  size_t block = parlay::binary_search(blocks_imap, i, lte);
  assert(block >= 0);
  assert(block < num_blocks);

  uchar* finger =
      (block > 0) ? (edge_start + block_offsets[block - 1]) : nghs_start;
  uintE start = *((uintE*)finger);
  finger += sizeof(uintE);
  uintE ngh = eatFirstEdge(finger, source);
  W wgh = eatWeight<W>(finger);
  if (i == start) {
    return std::make_tuple(ngh, wgh);
  }
  for (size_t edgeID = start + 1; edgeID < i; edgeID++) {
    ngh += eatEdge(finger);
    wgh = eatWeight<W>(finger);
  }
  return std::make_tuple(ngh, wgh);
}

uintE get_num_blocks(uchar* edge_start, uintE degree);

uintE get_block_degree(uchar* edge_start, uintE degree, uintE block_num);

template <class W>
inline void repack_sequential(const uintE& source, const uintE& degree,
                              uchar* edge_start) {
  uintE virtual_degree = *((uintE*)edge_start);
  size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
  uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));
  uchar* nghs_start =
      edge_start + (num_blocks - 1) * sizeof(uintE) + sizeof(uintE);

  size_t new_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
  uintE offs_stack[100];
  uintE* offs = offs_stack;

  size_t num_in_cur_block = 0;
  size_t cur_block_size = 0;
  size_t cur_block_id = 0;
  uchar tmp[16];
  uintE last_ngh = 0;
  auto compress_next_edge = [&](const uintE& ngh, const W& wgh) {
    if (num_in_cur_block == 0) {
      uintE off = compressFirstEdge(tmp, 0, source, ngh);
      cur_block_size += off;
      off = compressWeight<W>(tmp, 0, wgh);
      cur_block_size += off;
    } else {
      uintE difference = ngh - last_ngh;
      uintE off = compressEdge(tmp, 0, difference);
      cur_block_size += off;
      off = compressWeight<W>(tmp, 0, wgh);
      cur_block_size += off;
    }
    if (num_in_cur_block == PARALLEL_DEGREE) {
      offs[cur_block_id++] = cur_block_size;
      cur_block_size = 0;
      num_in_cur_block = 0;
    }
    last_ngh = ngh;
  };

  // 1. Compute the size of each block
  for (size_t i = 0; i < num_blocks; i++) {
    uchar* finger = (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
    uintE* block_deg_ptr = (uintE*)finger;
    uintE start_offset = *block_deg_ptr;
    uintE end_offset = (i == (num_blocks - 1))
                           ? degree
                           : (*((uintE*)(edge_start + block_offsets[i])));
    uintE block_deg = end_offset - start_offset;
    finger += sizeof(uintE);

    if (block_deg > 0) {
      uintE ngh = eatFirstEdge(finger, source);
      W wgh = eatWeight<W>(finger);
      compress_next_edge(ngh, wgh);
      for (size_t edgeID = start_offset + 1; edgeID < end_offset; edgeID++) {
        // Eat the next 'edge', which is a difference, and reconstruct edge.
        ngh += eatEdge(finger);
        wgh = eatWeight<W>(finger);
        compress_next_edge(ngh, wgh);
      }
    }
  }

  // 2. Scan to compute block offsets
  auto bytes_imap = gbbs::make_slice(offs, new_blocks + 1);
  parlay::scan_inplace(bytes_imap);

  // 3. Compress each block
  nghs_start = edge_start + (new_blocks - 1) * sizeof(uintE) + sizeof(uintE);
  size_t current_offset =
      sizeof(uintE) +
      (new_blocks - 1) * sizeof(uintE);  // virtual_deg + block_offs
  num_in_cur_block = 0;
  cur_block_id = 0;  // reset counters
  size_t num_finished = 0;
  auto write_next_edge = [&](const uintE& ngh, const W& wgh) {
    if (num_in_cur_block == 0) {  // start of a new block
      if (cur_block_id > 0) {
        // write the block_offset from offs
        block_offsets[cur_block_id - 1] = offs[cur_block_id];
      }
      uintE* block_deg = (uintE*)(edge_start + current_offset);
      if (cur_block_id < new_blocks) {
        *block_deg = PARALLEL_DEGREE;
      } else {
        size_t num_left = degree - num_finished;
        *block_deg = num_left;
      }
      current_offset += sizeof(uintE);  // block_deg
      current_offset =
          compressFirstEdge(edge_start, current_offset, source, ngh);
      current_offset = compressWeight<W>(edge_start, current_offset, wgh);
    } else {
      uintE difference = ngh - last_ngh;
      current_offset = compressEdge(edge_start, current_offset, difference);
      current_offset = compressWeight<W>(edge_start, current_offset, wgh);
    }
    if (num_in_cur_block == PARALLEL_DEGREE) {
      cur_block_id++;
      num_in_cur_block = 0;
      num_finished += PARALLEL_DEGREE;
    }
    last_ngh = ngh;
  };

  for (size_t i = 0; i < num_blocks; i++) {
    uchar* finger = (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
    uintE* block_deg_ptr = (uintE*)finger;
    uintE start_offset = *block_deg_ptr;
    uintE end_offset = (i == (num_blocks - 1))
                           ? degree
                           : (*((uintE*)(edge_start + block_offsets[i])));
    uintE block_deg = end_offset - start_offset;
    finger += sizeof(uintE);

    if (block_deg > 0) {
      uintE ngh = eatFirstEdge(finger, source);
      W wgh = eatWeight<W>(finger);
      write_next_edge(ngh, wgh);
      for (size_t edgeID = start_offset + 1; edgeID < end_offset; edgeID++) {
        // Eat the next 'edge', which is a difference, and reconstruct edge.
        ngh += eatEdge(finger);
        wgh = eatWeight<W>(finger);
        write_next_edge(ngh, wgh);
      }
    }
  }
}

template <class W>
inline void repack(const uintE& source, const uintE& degree, uchar* edge_start,
                   std::tuple<uintE, W>* tmp_space, bool par = true) {
  // No need to repack if degree == 0; all other methods abort when the vertex
  // degree is 0.
  if (degree > 0) {
    uintE virtual_degree = *((uintE*)edge_start);
    size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));
    uchar* nghs_start =
        edge_start + (num_blocks - 1) * sizeof(uintE) + sizeof(uintE);

    // 1. Copy all live edges into U
    using uintEW = std::tuple<uintE, W>;
    uintEW tmp_stack[100];
    uintEW* U = tmp_stack;
    parlay::sequence<uintEW> alloc;
    if (degree > 100) {
      alloc = parlay::sequence<uintEW>::uninitialized(degree);
      U = alloc.begin();
    }
    parallel_for(0, num_blocks, 2, [&](size_t i) {
      uchar* finger =
          (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
      uintE start_offset = *((uintE*)finger);
      uintE end_offset = (i == (num_blocks - 1))
                             ? degree
                             : (*((uintE*)(edge_start + block_offsets[i])));
      finger += sizeof(uintE);

      if (start_offset < end_offset) {
        uintE ngh = eatFirstEdge(finger, source);
        W wgh = eatWeight<W>(finger);
        U[start_offset] = std::make_tuple(ngh, wgh);
        for (size_t edgeID = start_offset + 1; edgeID < end_offset; edgeID++) {
          // Eat the next 'edge', which is a difference, and reconstruct edge.
          ngh += eatEdge(finger);
          wgh = eatWeight<W>(finger);
          U[edgeID] = std::make_tuple(ngh, wgh);
        }
      }
    });

    // 2. Repack from edge_start
    size_t new_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE offs_stack[100];
    uintE* offs = offs_stack;
    parlay::sequence<uintE> offs_alloc;
    if ((new_blocks + 1) > 100) {
      offs_alloc = parlay::sequence<uintE>::uninitialized(new_blocks + 1);
      offs = offs_alloc.begin();
    }

    // 3. Compute #bytes per new block
    parallel_for(0, new_blocks, 2, [&](size_t i) {
      size_t start = i * PARALLEL_DEGREE;
      size_t end = start + std::min<size_t>(PARALLEL_DEGREE, degree - start);
      uintE bytes = 0;
      bytes += sizeof(uintE);  // block_deg
      uchar tmp[16];
      auto nw = U[start];
      uintE off = compressFirstEdge(tmp, 0, source, std::get<0>(nw));
      bytes += off;
      off = compressWeight<W>(tmp, 0, std::get<1>(nw));
      bytes += off;
      for (size_t edgeI = start + 1; edgeI < end; edgeI++) {
        uintE difference = std::get<0>(U[edgeI]) - std::get<0>(U[edgeI - 1]);
        off = compressEdge(tmp, 0, difference);
        bytes += off;
        off = compressWeight<W>(tmp, 0, std::get<1>(U[edgeI]));
        bytes += off;
      }
      offs[i] = bytes;
    });

    // 4. Scan to compute offset for each block
    offs[new_blocks] = 0;
    auto bytes_imap = make_slice(offs, offs + new_blocks + 1);
    parlay::scan_inplace(bytes_imap);

    // 5. Repack each block
    uintE* virtual_degree_ptr = (uintE*)edge_start;
    *virtual_degree_ptr = degree;  // update the virtual degree
    // block_offsets are unchanged
    nghs_start = edge_start + (new_blocks - 1) * sizeof(uintE) +
                 sizeof(uintE);  // update ngh_start
    parallel_for(0, new_blocks, 2, [&](size_t i) {
      size_t start = i * PARALLEL_DEGREE;
      size_t end = start + std::min<size_t>(PARALLEL_DEGREE, degree - start);
      uchar* finger = nghs_start + bytes_imap[i];
      // update block offsets with the distance from the start
      if (i > 0) {
        block_offsets[i - 1] = finger - edge_start;
      }

      // write the edge offset for this block
      uintE* block_offset = (uintE*)finger;
      *block_offset = start;
      size_t current_offset = sizeof(uintE);

      auto first_nw = U[start];
      current_offset = compressFirstEdge(finger, current_offset, source,
                                         std::get<0>(first_nw));
      current_offset =
          compressWeight<W>(finger, current_offset, std::get<1>(first_nw));
      uintE last_ngh = std::get<0>(first_nw);
      for (size_t j = start + 1; j < end; j++) {
        auto nw = U[j];
        uintE difference = std::get<0>(nw) - last_ngh;
        current_offset = compressEdge(finger, current_offset, difference);
        current_offset =
            compressWeight<W>(finger, current_offset, std::get<1>(nw));
        last_ngh = std::get<0>(nw);
      }
    });
  }
}

template <class W, class P>
inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                   const uintE& degree, std::tuple<uintE, W>* tmp_space,
                   bool par = true) {
  using uintEW = std::tuple<uintE, W>;
  uintE virtual_degree = *((uintE*)edge_start);
  size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;

  uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));
  uchar* nghs_start = edge_start + (num_blocks - 1) * sizeof(uintE) +
                      sizeof(uintE);  // block offs + virtual_degree

  size_t block_cts_stack[100];
  parlay::sequence<size_t> alloc;
  size_t* block_cts = block_cts_stack;
  if (num_blocks > 100) {
    alloc = parlay::sequence<size_t>::uninitialized(num_blocks + 1);
    block_cts = alloc.begin();
  }

  parallel_for(0, num_blocks, 2, [&](size_t i) {
    uchar* finger = (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
    uintE* block_deg_ptr = (uintE*)finger;
    uintE start_offset = *block_deg_ptr;
    uintE end_offset = (i == (num_blocks - 1))
                           ? degree
                           : (*((uintE*)(edge_start + block_offsets[i])));
    uintE block_deg = end_offset - start_offset;
    finger += sizeof(uintE);

    // TODO(laxmand): why did I write this to compute block_cts in one pass,
    // and then compact in the second pass?
    // A) Uncompress and filter edges into tmp
    uintEW tmp[PARALLEL_DEGREE];
    size_t ct = 0;
    debug(size_t final_off =
              finger -
              ((i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start););
    if (block_deg > 0) {
      uintE ngh = eatFirstEdge(finger, source);
      W wgh = eatWeight<W>(finger);
      if (pred(source, ngh, wgh)) {
        tmp[ct++] = std::make_tuple(ngh, wgh);
      }
      for (size_t j = 1; j < block_deg; j++) {
        ngh += eatEdge(finger);
        wgh = eatWeight<W>(finger);
        if (pred(source, ngh, wgh)) {
          tmp[ct++] = std::make_tuple(ngh, wgh);
        }
      }
      debug(final_off = finger - ((i > 0) ? (edge_start + block_offsets[i - 1])
                                          : nghs_start););
    }
    // B) write the number of live edges in this block to block_cts
    block_cts[i] = ct;

    // C) Recompress inside this block
    size_t offset = 0;
    if (ct > 0 && ct < block_deg) {
      uchar* finger2 =
          (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
      uchar* write_finger = finger2 + sizeof(uintE);
      offset =
          compressFirstEdge(write_finger, offset, source, std::get<0>(tmp[0]));
      offset = compressWeight<W>(write_finger, offset, std::get<1>(tmp[0]));
      uintE last_ngh = std::get<0>(tmp[0]);
      for (size_t j = 1; j < ct; j++) {
        const auto& e = tmp[j];
        uintE difference = std::get<0>(e) - last_ngh;
        offset = compressEdge(write_finger, offset, difference);
        offset = compressWeight<W>(write_finger, offset, std::get<1>(e));
        last_ngh = std::get<0>(e);
      }
    }

    assert(offset <= final_off);
  });

  // 2. Scan block_cts to get offsets within blocks
  block_cts[num_blocks] = 0;
  auto scan_cts = gbbs::make_slice(block_cts, num_blocks + 1);
  size_t deg_remaining = parlay::scan_inplace(scan_cts);

  parallel_for(0, num_blocks, 1000, [&](size_t i) {
    uchar* finger = (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
    uintE* block_deg_ptr = (uintE*)finger;
    *block_deg_ptr = scan_cts[i];
  });

  // Can comment out this call to avoid repacking; this can make algorithms,
  // e.g. set-cover no longer theoreticaly efficient
  if (deg_remaining < (virtual_degree / 10)) {
    repack<W>(source, deg_remaining, edge_start, tmp_space, par);
  }

  return deg_remaining;
}

template <class W>
inline void decode_block(uchar* finger, std::tuple<uintE, W>* out, size_t start,
                         size_t end, const uintE& source) {
  if (end - start > 0) {
    uintE ngh = eatFirstEdge(finger, source);
    {
      W wgh = eatWeight<W>(finger);
      out[start] = std::make_tuple(ngh, wgh);
    }
    for (size_t i = start + 1; i < end; i++) {
      // Eat the next 'edge', which is a difference, and reconstruct edge.
      ngh += eatEdge(finger);
      W wgh = eatWeight<W>(finger);
      out[i] = std::make_tuple(ngh, wgh);
    }
  }
}

template <class W, class P, class O>
inline void filter_sequential(P pred, uchar* edge_start, const uintE& source,
                              const uintE& degree, O& out) {
  uintE virtual_degree = *((uintE*)edge_start);
  size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
  uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));
  uchar* nghs_start = edge_start + (num_blocks - 1) * sizeof(uintE) +
                      sizeof(uintE);  // block offs + virtual_degree

  size_t k = 0;
  for (size_t i = 0; i < num_blocks; i++) {
    uchar* finger = (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
    uintE start_offset = *((uintE*)finger);
    uintE end_offset = (i == (num_blocks - 1))
                           ? degree
                           : (*((uintE*)(edge_start + block_offsets[i])));
    finger += sizeof(uintE);
    if (start_offset < end_offset) {
      uintE ngh = eatFirstEdge(finger, source);
      W wgh = eatWeight<W>(finger);
      if (pred(source, ngh, wgh)) {
        out(k++, std::make_tuple(ngh, wgh));
      }
      for (size_t edgeID = start_offset + 1; edgeID < end_offset; edgeID++) {
        ngh += eatEdge(finger);
        wgh = eatWeight<W>(finger);
        if (pred(source, ngh, wgh)) {
          out(k++, std::make_tuple(ngh, wgh));
        }
      }
    }
  }
}

template <class W, class P, class O>
inline void filter(P pred, uchar* edge_start, const uintE& source,
                   const uintE& degree, std::tuple<uintE, W>* tmp, O& out) {
  if (degree <= PD_PACK_THRESHOLD && degree > 0) {
    filter_sequential<W, P, O>(pred, edge_start, source, degree, out);
  } else if (degree > 0) {
    uintE virtual_degree = *((uintE*)edge_start);
    size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));
    uchar* nghs_start = edge_start + (num_blocks - 1) * sizeof(uintE) +
                        sizeof(uintE);  // block offs + virtual_degree

    size_t tmp_size = degree / kTemporarySpaceConstant;
    size_t blocks_per_iter = tmp_size / PARALLEL_DEGREE;
    size_t blocks_finished = 0, out_off = 0;

    while (blocks_finished < num_blocks) {
      size_t start_block = blocks_finished;
      size_t end_block = std::min(start_block + blocks_per_iter, num_blocks);
      size_t total_blocks = end_block - start_block;

      uchar* first_finger = (start_block > 0)
                                ? (edge_start + block_offsets[start_block - 1])
                                : nghs_start;
      uintE first_offset = *((uintE*)first_finger);
      size_t last_offset = 0;

      parallel_for(start_block, end_block, 1, [&](size_t i) {
        uchar* finger =
            (i > 0) ? (edge_start + block_offsets[i - 1]) : nghs_start;
        uintE start_offset = *((uintE*)finger) - first_offset;
        uintE end_offset =
            ((i == (num_blocks - 1))
                 ? degree
                 : (*((uintE*)(edge_start + block_offsets[i])))) -
            first_offset;
        if (i == (end_block - 1)) {
          last_offset = end_offset;
        }
        finger += sizeof(uintE);
        decode_block(finger, tmp, start_offset, end_offset, source);
      });

      // filter edges into tmp2
      auto pd = [&](const std::tuple<uintE, W>& nw) {
        return pred(source, std::get<0>(nw), std::get<1>(nw));
      };
      uintE k = parlay::filterf(tmp, last_offset, pd, out, out_off);
      out_off += k;

      blocks_finished += total_blocks;
    }
  }
}

template <class W, class I>
inline long sequentialCompressEdgeSet(uchar* edgeArray, size_t current_offset,
                                      uintT degree, uintE source, I& it,
                                      size_t encoded_degree = PARALLEL_DEGREE) {
  if (degree > 0) {
    size_t start_offset = current_offset;
    size_t num_blocks = 1 + (degree - 1) / encoded_degree;
    uintE* vertex_ctr = (uintE*)edgeArray;
    *vertex_ctr = degree;
    uintE* block_offsets = (uintE*)(edgeArray + sizeof(uintE));
    current_offset +=
        sizeof(uintE) +
        (num_blocks - 1) * sizeof(uintE);  // virtual deg + block_offs
    for (size_t i = 0; i < num_blocks; i++) {
      size_t o = i * encoded_degree;
      size_t end = std::min<size_t>(encoded_degree, degree - o);

      if (i > 0)
        block_offsets[i - 1] =
            current_offset -
            start_offset;  // store offset for all chunks but the first
      uintE* block_deg = (uintE*)(edgeArray + current_offset);
      *block_deg = o;
      current_offset += sizeof(uintE);

      std::tuple<uintE, W> lst = (i == 0) ? it.cur() : it.next();
      uintE last_ngh = std::get<0>(lst);

      uchar* test_fing = edgeArray + current_offset;
      current_offset =
          compressFirstEdge(edgeArray, current_offset, source, last_ngh);
      auto decf = eatFirstEdge(test_fing, source);
      if (decf != last_ngh) {
        std::cout << "# first enc: " << source << " and " << last_ngh
                  << " got back " << decf << "\n";
      }
      current_offset =
          compressWeight<W>(edgeArray, current_offset, std::get<1>(lst));
      for (size_t edgeI = 1; edgeI < end; edgeI++) {
        std::tuple<uintE, W> nxt = it.next();
        uintE difference = std::get<0>(nxt) - last_ngh;

        test_fing = edgeArray + current_offset;
        current_offset = compressEdge(edgeArray, current_offset, difference);
        auto dec = eatEdge(test_fing);
        if (dec != difference) {
          std::cout << "# src = " << source << " enc = " << std::get<0>(nxt)
                    << "\n";
        }
        current_offset =
            compressWeight<W>(edgeArray, current_offset, std::get<1>(nxt));
        last_ngh = std::get<0>(nxt);
      }
    }
  }
  return current_offset;
}
}  // namespace bytepd_amortized
}  // namespace gbbs
