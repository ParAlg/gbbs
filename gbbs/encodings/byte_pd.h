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

#include "gbbs/bridge.h"
#include "gbbs/macros.h"

namespace gbbs {
namespace bytepd {

inline size_t get_virtual_degree(uintE d, uchar* ngh_arr) { return d; }

template <class W>
inline std::tuple<uintE, W> get_ith_neighbor(uchar* edge_start, uintE source,
                                             uintE degree, size_t i) {}

uintE get_num_blocks(uchar* edge_start, uintE degree);

// Read default weight (expects gbbs::empty)
template <class W,
          typename std::enable_if<!std::is_same<W, intE>::value, int>::type = 0>
__attribute__((always_inline)) inline W eatWeight(uchar*& start) {
  return (W)gbbs::empty();
}

template <class W,
          typename std::enable_if<std::is_same<W, intE>::value, int>::type = 0>
__attribute__((always_inline)) inline void print_weight(W& wgh) {
  std::cout << wgh << "\n";
}

template <class W,
          typename std::enable_if<!std::is_same<W, intE>::value, int>::type = 0>
inline void print_weight(W& wgh) {}

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

inline intE eatFirstEdge(uchar*& start, const uintE source) {
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

/*
  Compresses the first edge, writing target-source and a sign bit.
*/
long compressFirstEdge(uchar* start, long offset, uintE source, uintE target);

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

long compressEdge(uchar* start, long curOffset, uintE e);

template <class W>
struct iter {
  uchar* base;
  uchar* finger;
  uintE src;
  uintT degree;

  uintE num_blocks;
  uintE cur_chunk;

  std::tuple<uintE, W> last_edge;
  uintE proc;

  iter(uchar* _base, uintT _degree, uintE _src)
      : base(_base), degree(_degree), src(_src) {
    num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    finger = base + (num_blocks - 1) * sizeof(uintE);
    cur_chunk = 0;
    if (degree == 0) {
      return;
    }

    std::get<0>(last_edge) = eatFirstEdge(finger, src);
    std::get<1>(last_edge) = eatWeight<W>(finger);
    proc = 1;
  }

  inline std::tuple<uintE, W> cur() { return last_edge; }

  inline std::tuple<uintE, W> next() {
    if (proc == PARALLEL_DEGREE) {
      cur_chunk++;
      std::get<0>(last_edge) = eatFirstEdge(finger, src);
      std::get<1>(last_edge) = eatWeight<W>(finger);
      proc = 1;
    } else {
      std::get<0>(last_edge) += eatEdge(finger);
      std::get<1>(last_edge) = eatWeight<W>(finger);
      proc++;
    }
    return last_edge;
  }

  inline bool has_next() {
    return (cur_chunk * PARALLEL_DEGREE + proc) < degree;
  }
};

template <class W>
struct simple_iter {
  uchar* base;
  uchar* finger;
  uintE src;
  uintT degree;

  uintE num_blocks;
  uintE cur_chunk;

  std::tuple<uintE, W> last_edge;
  uintE proc;

  simple_iter(uchar* _base, uintT _degree, uintE _src)
      : base(_base), degree(_degree), src(_src) {
    num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    finger = base + (num_blocks - 1) * sizeof(uintE);
    cur_chunk = 0;
    if (degree == 0) {
      return;
    }

    std::get<0>(last_edge) = eatFirstEdge(finger, src);
    std::get<1>(last_edge) = eatWeight<W>(finger);
    proc = 1;
  }

  inline std::tuple<uintE, W> cur() { return last_edge; }

  inline std::tuple<uintE, W> next() {
    if (proc == PARALLEL_DEGREE) {
      cur_chunk++;
      std::get<0>(last_edge) = eatFirstEdge(finger, src);
      std::get<1>(last_edge) = eatWeight<W>(finger);
      proc = 1;
    } else {
      std::get<0>(last_edge) += eatEdge(finger);
      std::get<1>(last_edge) = eatWeight<W>(finger);
      proc++;
    }
    return last_edge;
  }

  inline bool has_next() {
    return (cur_chunk * PARALLEL_DEGREE + proc) < degree;
  }
};

// Decode unweighted edges
template <
    class W, class T,
    typename std::enable_if<std::is_same<W, gbbs::empty>::value, int>::type = 0>
inline void decode(T t, uchar* edge_start, const uintE& source,
                   const uintT& degree, const bool par = true) {
  if (degree > 0) {
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)
        edge_start;  // use beginning of edgeArray for offsets into edge list
    uchar* first_finger = edge_start + (num_blocks - 1) * sizeof(uintE);
    // do first chunk
    size_t block_end = std::min<long>(PARALLEL_DEGREE, degree);

    // Eat first edge, which is compressed specially
    {
      uintE ngh = eatFirstEdge(first_finger, source);
      if (!t(source, ngh, gbbs::empty(), 0)) return;
      for (uintE edge_id = 1; edge_id < block_end; edge_id++) {
        // Eat the next 'edge', which is a difference, and reconstruct edge.
        ngh += eatEdge(first_finger);
        if (!t(source, ngh, gbbs::empty(), edge_id)) return;
      }
    }
    // do remaining chunks in parallel
    parallel_for(1, num_blocks, 1, [&](size_t i) {
      size_t o = i * PARALLEL_DEGREE;
      size_t end = std::min<long>(o + PARALLEL_DEGREE, degree);
      uchar* finger = edge_start + block_offsets[i - 1];
      // Eat first edge, which is compressed specially
      uintE ngh = eatFirstEdge(finger, source);
      if (!t(source, ngh, gbbs::empty(), o)) end = 0;
      for (size_t edge_id = o + 1; edge_id < end; edge_id++) {
        // Eat the next 'edge', which is a difference, and reconstruct edge.
        ngh += eatEdge(finger);
        if (!t(source, ngh, gbbs::empty(), edge_id)) break;
      }
    });
  }
}

// Decode weighted edges
template <class W, class T,
          typename std::enable_if<!std::is_same<W, gbbs::empty>::value,
                                  int>::type = 0>
inline void decode(T t, uchar* edge_start, const uintE& source,
                   const uintT& degree, const bool par = true) {
  if (degree > 0) {
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)
        edge_start;  // use beginning of edgeArray for offsets into edge list
    uchar* first_finger = edge_start + (num_blocks - 1) * sizeof(uintE);
    // do first chunk
    size_t block_end = std::min<long>(PARALLEL_DEGREE, degree);

    // Eat first edge, which is compressed specially
    {
      uintE ngh = eatFirstEdge(first_finger, source);
      W first_weight = eatWeight<W>(first_finger);
      if (!t(source, ngh, first_weight, 0)) return;
      for (size_t edge_id = 1; edge_id < block_end; edge_id++) {
        // Eat the next 'edge', which is a difference, and reconstruct edge.
        ngh += eatEdge(first_finger);
        W weight = eatWeight<W>(first_finger);
        if (!t(source, ngh, weight, edge_id)) return;
      }
    }
    // do remaining chunks in parallel
    parallel_for(1, num_blocks, 1, [&](size_t i) {
      size_t o = i * PARALLEL_DEGREE;
      size_t end = std::min<long>(o + PARALLEL_DEGREE, degree);
      uchar* finger = edge_start + block_offsets[i - 1];
      // Eat first edge, which is compressed specially
      uintE ngh = eatFirstEdge(finger, source);
      W weight = eatWeight<W>(finger);
      if (!t(source, ngh, weight, o)) end = 0;
      for (size_t edge_id = o + 1; edge_id < end; edge_id++) {
        // Eat the next 'edge', which is a difference, and reconstruct edge.
        ngh += eatEdge(finger);
        weight = eatWeight<W>(finger);
        if (!t(source, ngh, weight, edge_id)) break;
      }
    });
  }
}

template <class W, class T>
void decode_block(T t, uchar* edge_start, const uintE& source,
                  const uintT& degree, uintE block_num) {
  if (degree > 0) {
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)
        edge_start;  // use beginning of edgeArray for offsets into edge list

    size_t offset = block_num * PARALLEL_DEGREE;
    size_t end = std::min<long>(offset + PARALLEL_DEGREE, degree);
    uchar* finger = (block_num == 0)
                        ? (edge_start + (num_blocks - 1) * sizeof(uintE))
                        : (edge_start + block_offsets[block_num - 1]);
    // Eat first edge, which is compressed specially
    uintE ngh = eatFirstEdge(finger, source);
    W weight = eatWeight<W>(finger);
    t(ngh, weight, offset);
    for (size_t edge_id = offset + 1; edge_id < end; edge_id++) {
      // Eat the next 'edge', which is a difference, and reconstruct edge.
      ngh += eatEdge(finger);
      weight = eatWeight<W>(finger);
      t(ngh, weight, edge_id);
    }
  }
}

template <class W, class T>
void decode_block_cond(T t, uchar* edge_start, const uintE& source,
                       const uintT& degree, uintE block_num) {
  if (degree > 0) {
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)
        edge_start;  // use beginning of edgeArray for offsets into edge list

    size_t offset = block_num * PARALLEL_DEGREE;
    size_t end = std::min<long>(offset + PARALLEL_DEGREE, degree);
    uchar* finger = (block_num == 0)
                        ? (edge_start + (num_blocks - 1) * sizeof(uintE))
                        : (edge_start + block_offsets[block_num - 1]);
    // Eat first edge, which is compressed specially
    uintE ngh = eatFirstEdge(finger, source);
    W weight = eatWeight<W>(finger);
    if (!t(ngh, weight, offset)) return;
    for (size_t edge_id = offset + 1; edge_id < end; edge_id++) {
      // Eat the next 'edge', which is a difference, and reconstruct edge.
      ngh += eatEdge(finger);
      weight = eatWeight<W>(finger);
      if (!t(ngh, weight, edge_id)) break;
    }
  }
}

template <class W, class T>
inline void decode_block_seq(T t, uchar* edge_start, const uintE& source,
                             const uintT& degree, uintE block_size,
                             uintE block_num) {
  assert(false);
}

uintE get_block_degree(uchar* edge_start, uintE degree, uintE block_num);

// Map_reduce
template <class W, class M, class Monoid>
inline typename Monoid::T map_reduce(uchar* edge_start, const uintE& source,
                                     const uintT& degree, M& m, Monoid& reduce,
                                     const bool par = true) {
  using E = typename Monoid::T;
  if (degree > 0) {
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)edge_start;

    E stk[1000];
    E* block_outputs;
    parlay::sequence<E> alloc;
    if (num_blocks > 1000) {
      alloc = parlay::sequence<E>::uninitialized(num_blocks);
      block_outputs = alloc.begin();
    } else {
      block_outputs = (E*)stk;
    }

    parallel_for(0, num_blocks, 1, [&](size_t i) {
      size_t start = i * PARALLEL_DEGREE;
      size_t end = std::min<long>(start + PARALLEL_DEGREE, degree);
      uchar* finger = edge_start + ((i == 0) ? (num_blocks - 1) * sizeof(uintE)
                                             : block_offsets[i - 1]);

      // Eat first edge, which is compressed specially
      uintE ngh = eatFirstEdge(finger, source);
      W wgh = eatWeight<W>(finger);
      E cur = reduce.identity;
      cur = reduce.f(cur, m(source, ngh, wgh));
      for (size_t j = start + 1; j < end; j++) {
        ngh += eatEdge(finger);
        wgh = eatWeight<W>(finger);
        cur = reduce.f(cur, m(source, ngh, wgh));
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

template <class W>
inline size_t compute_block_size(std::tuple<uintE, W>* edges, size_t start,
                                 size_t end, const uintE& source) {
  uchar tmp[16];
  size_t total_bytes = 0;
  auto nw = edges[start];
  total_bytes += compressFirstEdge(tmp, 0, source, std::get<0>(nw));
  total_bytes += compressWeight<W>(tmp, 0, std::get<1>(nw));
  uintE last_ngh = std::get<0>(nw);
  for (size_t i = start + 1; i < end; i++) {
    nw = edges[i];
    uintE difference = std::get<0>(nw) - last_ngh;
    total_bytes += compressEdge(tmp, 0, difference);
    total_bytes += compressWeight<W>(tmp, 0, std::get<1>(nw));
    last_ngh = std::get<0>(nw);
  }
  return total_bytes;
}

template <class W>
inline size_t compute_size_in_bytes(std::tuple<uintE, W>* edges,
                                    const uintE& source, const uintE& d) {
  if (d > 0) {
    size_t num_blocks = 1 + (d - 1) / PARALLEL_DEGREE;
    uintE stk[100];
    uintE* block_bytes;
    parlay::sequence<uintE> alloc;
    if (num_blocks > 100) {
      alloc = parlay::sequence<uintE>::uninitialized(num_blocks);
      block_bytes = alloc.begin();
    } else {
      block_bytes = (uintE*)stk;
    }
    parallel_for(0, num_blocks, 1, [&](size_t i) {
      // calculate size in bytes of this block
      uintE start = i * PARALLEL_DEGREE;
      uintE end = start + std::min<uintE>(PARALLEL_DEGREE, d - start);
      block_bytes[i] = compute_block_size(edges, start, end, source);
    });
    auto bytes_imap = gbbs::make_slice(block_bytes, num_blocks);
    size_t total_space = parlay::scan_inplace(bytes_imap);

    // add in space for storing offsets to the start of each block
    total_space += sizeof(uintE) * (num_blocks - 1);
    return total_space;
  } else {
    return 0;
  }
}

template <class W, class I>
inline long sequentialCompressEdgeSet(uchar* edgeArray, size_t current_offset,
                                      uintT degree, uintE source, I& it,
                                      size_t encoded_degree = PARALLEL_DEGREE) {
  if (degree > 0) {
    size_t start_offset = current_offset;
    size_t num_blocks = 1 + (degree - 1) / encoded_degree;
    uintE* block_offsets = (uintE*)edgeArray;
    current_offset += (num_blocks - 1) * sizeof(uintE);  // block_offs
    for (size_t i = 0; i < num_blocks; i++) {
      size_t o = i * encoded_degree;
      size_t end = std::min<size_t>(encoded_degree, degree - o);

      if (i > 0)
        block_offsets[i - 1] =
            current_offset -
            start_offset;  // store offset for all chunks but the first

      std::tuple<uintE, W> lst = (i == 0) ? it.cur() : it.next();
      uintE last_ngh = std::get<0>(lst);
      current_offset =
          compressFirstEdge(edgeArray, current_offset, source, last_ngh);
      current_offset =
          compressWeight<W>(edgeArray, current_offset, std::get<1>(lst));
      for (size_t edgeI = 1; edgeI < end; edgeI++) {
        std::tuple<uintE, W> nxt = it.next();
        uintE difference = std::get<0>(nxt) - last_ngh;
        current_offset = compressEdge(edgeArray, current_offset, difference);
        current_offset =
            compressWeight<W>(edgeArray, current_offset, std::get<1>(nxt));
        last_ngh = std::get<0>(nxt);
      }
    }
  }
  return current_offset;
}

// Packs out the block sequentially.
template <class W, class P>
inline size_t packBlock(uchar* edge_start, const uintE& source,
                        const uintE& degree, P& pred) {
  uchar* finger = edge_start;
  size_t ct = 0;

  std::tuple<uintE, W> tmp[PARALLEL_DEGREE];
  uintE ngh = eatFirstEdge(finger, source);
  W wgh = eatWeight<W>(finger);
  if (pred(source, ngh, wgh)) {
    tmp[ct++] = std::make_tuple(ngh, wgh);
  }
  for (size_t i = 1; i < degree; i++) {
    ngh += eatEdge(finger);
    wgh = eatWeight<W>(finger);
    if (pred(source, ngh, wgh)) {
      tmp[ct++] = std::make_tuple(ngh, wgh);
    }
  }

  if (ct == 0) {
    return 0;
  }

  size_t offset = 0;  // the write offset
  offset = compressFirstEdge(edge_start, offset, source, std::get<0>(tmp[0]));
  offset = compressWeight<W>(edge_start, offset, std::get<1>(tmp[0]));
  uintE last_ngh = std::get<0>(tmp[0]);
  for (size_t i = 1; i < ct; i++) {
    const auto& e = tmp[i];
    uintE difference = std::get<0>(e) - last_ngh;
    offset = compressEdge(edge_start, offset, difference);
    offset = compressWeight<W>(edge_start, offset, std::get<1>(e));
    last_ngh = std::get<0>(e);
  }
  assert(offset <= (finger - edge_start));

  //  size_t ctt = map_reduce<W>(edge_start, source, ct, 0, map_f, red_f, true);

  return ct;
}

template <class W, class P, class O>
inline size_t filterBlock(P pred, uchar* edge_start, const uintE& source,
                          const uintE& degree, O& out) {
  size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
  size_t n_offsets = num_blocks - 1;
  uintE edges_offset = (n_offsets) * sizeof(uintE);
  uintE* block_offsets = (uintE*)edge_start;
  size_t k = 0;
  for (size_t i = 0; i < num_blocks; i++) {
    uchar* finger =
        edge_start + ((i == 0) ? edges_offset : block_offsets[i - 1]);
    uintE start = i * PARALLEL_DEGREE;
    uintE end = std::min<uintE>(start + PARALLEL_DEGREE, degree);

    uintE ngh = eatFirstEdge(finger, source);
    W wgh = eatWeight<W>(finger);
    if (pred(source, ngh, wgh)) {
      out(k++, std::make_tuple(ngh, wgh));
    }
    for (size_t j = start + 1; j < end; j++) {
      ngh += eatEdge(finger);
      wgh = eatWeight<W>(finger);
      if (pred(source, ngh, wgh)) {
        out(k++, std::make_tuple(ngh, wgh));
      }
    }
  }
  return k;
}

template <class W>
inline void compress_edges(uchar* edgeArray, const uintE& source,
                           const uintE& d, std::tuple<uintE, W>* edges,
                           uchar* last_finger, bool par = true) {
  size_t num_blocks = 1 + (d - 1) / PARALLEL_DEGREE;
  uintE stk[101];
  uintE* block_bytes;
  parlay::sequence<uintE> alloc;
  if (num_blocks > 100) {
    alloc = parlay::sequence<uintE>::uninitialized(num_blocks + 1);
    block_bytes = alloc.begin();
  } else {
    block_bytes = (uintE*)stk;
  }

  parallel_for(0, num_blocks, 1, [&](size_t i) {
    uintE start = i * PARALLEL_DEGREE;
    uintE end = std::min<uintE>(start + PARALLEL_DEGREE, d);
    block_bytes[i] = compute_block_size(edges, start, end, source);
  });
  block_bytes[num_blocks] = 0;

  uintE edges_offset = (num_blocks - 1) * sizeof(uintE);
  uintE* block_offsets = (uintE*)edgeArray;

  auto bytes_imap = gbbs::make_slice(block_bytes, num_blocks + 1);
  uintE total_space = parlay::scan_inplace(bytes_imap);

  if (total_space > (last_finger - edgeArray)) {
    std::cout << "# Space error!"
              << "\n";
    exit(0);
  }
  if (total_space == (last_finger - edgeArray)) {
    std::cout << "# d = " << d << " to exactly the same space"
              << "\n";
  }

  parallel_for(1, num_blocks, 1, [&](size_t i) {
    // store offset to start of this block from edgeArray
    block_offsets[i - 1] = edges_offset + bytes_imap[i];
  });

  parallel_for(0, num_blocks, 1, [&](size_t i) {
    uintE start = i * PARALLEL_DEGREE;
    uintE end = std::min<uintE>(start + PARALLEL_DEGREE, d);
    long write_offset = (i == 0) ? edges_offset : block_offsets[i - 1];
    auto nw = edges[start];
    write_offset =
        compressFirstEdge(edgeArray, write_offset, source, std::get<0>(nw));
    write_offset = compressWeight<W>(edgeArray, write_offset, std::get<1>(nw));
    uintE last_edge = std::get<0>(nw);
    for (size_t j = start + 1; j < end; j++) {
      nw = edges[j];
      write_offset =
          compressEdge(edgeArray, write_offset, std::get<0>(nw) - last_edge);
      write_offset =
          compressWeight<W>(edgeArray, write_offset, std::get<1>(nw));
      last_edge = std::get<0>(nw);
    }
  });
  if (num_blocks > 100) {
    gbbs::free_array(block_bytes, num_blocks + 1);
  }
}

template <class W, class P>
inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                   const uintE& degree, std::tuple<uintE, W>* tmp,
                   bool parallel = true) {
  if (degree <= PARALLEL_DEGREE) {
    // Pack small degree vertices sequentially
    // return packSequential<W, P>(edge_start, source, degree, pred);
    return packBlock<W, P>(edge_start, source, degree, pred);
  } else {
    return degree;

    // Pack large degree vertices in parallel
    // Decode nghs into temporary space, filter, then recompress.
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)edge_start;

    std::tuple<uintE, W> tmp_A[200000];
    std::tuple<uintE, W>* our_tmp = (std::tuple<uintE, W>*)tmp_A;
    if (degree > PD_PACK_THRESHOLD) {
      our_tmp = tmp;
    }

    uchar* mf = 0;
    parallel_for(0, num_blocks, 1, [&](size_t i) {
      size_t start = i * PARALLEL_DEGREE;
      size_t end = std::min<long>(start + PARALLEL_DEGREE, degree);
      uchar* finger = edge_start + ((i == 0) ? (num_blocks - 1) * sizeof(uintE)
                                             : block_offsets[i - 1]);
      uchar* last_finger = decode_block(finger, our_tmp, start, end, source);
      if (i == num_blocks - 1) {
        mf = last_finger;
      }
    });

    // filter edges into tmp2
    std::tuple<uintE, W>* tmp2 = our_tmp + degree;
    auto pd = [&](const std::tuple<uintE, W>& nw) {
      return pred(source, std::get<0>(nw), std::get<1>(nw));
    };
    uintE k = parlay::filterf(our_tmp, tmp2, degree, pd);
    if (k == degree || k == 0) {
      return k;
    }

    for (size_t i = 0; i < k; i++) {
      assert(std::get<0>(tmp2[i]) < 3072627);
    }

    // pack blocks out in parallel
    compress_edges<W>(edge_start, source, k, tmp2, mf);

    //    ct = map_reduce<W>(edge_start, source, k, 0, map_f, red_f, true);
    //    assert(ct == k);

    return k;
  }
}

template <class W, class P, class O>
inline size_t filter(P pred, uchar* edge_start, const uintE& source,
                     const uintE& degree, std::tuple<uintE, W>* tmp, O& out) {
  return filterBlock<W>(pred, edge_start, source, degree, out);
  if (degree <= PD_PACK_THRESHOLD && degree > 0) {
    // Filter small degree vertices sequentially
    uchar* finger = edge_start;
    uintE ngh = eatFirstEdge(finger, source);
    W wgh = eatWeight<W>(finger);
    uintE k = 0;
    if (pred(source, ngh, wgh)) {
      out(k++, std::make_tuple(ngh, wgh));
    }
    for (size_t i = 1; i < degree; i++) {
      ngh += eatEdge(finger);
      wgh = eatWeight<W>(finger);
      if (pred(source, ngh, wgh)) {
        out(k++, std::make_tuple(ngh, wgh));
      }
    }
    return k;
  } else {
    // Pack large degree vertices in parallel
    // Decode nghs into temporary space, then filter
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)edge_start;

    parallel_for(0, num_blocks, 1, [&](size_t i) {
      size_t start = i * PARALLEL_DEGREE;
      size_t end = std::min<long>(start + PARALLEL_DEGREE, degree);
      uchar* finger = edge_start + ((i == 0) ? (num_blocks - 1) * sizeof(uintE)
                                             : block_offsets[i - 1]);
      decode_block(finger, tmp, start, end, source);
    });

    // filter edges into tmp2
    std::tuple<uintE, W>* tmp2 = tmp + degree;
    auto pd = [&](const std::tuple<uintE, W>& nw) {
      return pred(source, std::get<0>(nw), std::get<1>(nw));
    };
    uintE k = parlay::filterf(tmp, tmp2, degree, pd);
    parallel_for(0, k, [&](size_t i) { out(i, tmp2[i]); });
    return k;
  }
  return 0;
}

}  // namespace bytepd
}  // namespace gbbs
