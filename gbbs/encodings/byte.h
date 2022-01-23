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
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <type_traits>

#include "gbbs/macros.h"

namespace gbbs {
namespace byte {

inline size_t get_virtual_degree(uintE d, uchar* nghArr) { return d; }

// Read an empty weight (noop)
template <class W, typename std::enable_if<std::is_same<W, gbbs::empty>::value,
                                           int>::type = 0>
inline gbbs::empty eatWeight(uchar*& start) {
  return gbbs::empty();
}

// Read an integer weight. Weights can be negative, so read using signed VarInt8
// logic.
template <class W,
          typename std::enable_if<std::is_same<W, intE>::value, int>::type = 0>
inline intE eatWeight(uchar*& start) {
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

// Read a float weight.
template <class W,
          typename std::enable_if<std::is_same<W, float>::value, int>::type = 0>
inline float eatWeight(uchar*& start) {
  float* float_start = (float*)start;
  start += sizeof(float);
  return *float_start;
}

// Eats the first (specially) encoded edge, which stores a sign bit. Should
// ideally be inlined, but let the compiler figure it out.
inline intE eatFirstEdge(uchar*& start, uintE source) {
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
  return (fb & 0x40) ? source - edgeRead : source + edgeRead;
}

// Reads positive edges after the first one.
inline uintE eatEdge(uchar*& start) {
  uintE edgeRead = 0;
  int shiftAmount = 0;

  while (1) {
    uchar b = *start++;
    edgeRead += ((b & 0x7f) << shiftAmount);
    if (LAST_BIT_SET(b))
      shiftAmount += EDGE_SIZE_PER_BYTE;
    else
      break;
  }
  return edgeRead;
}

template <class W, class T>
inline void decode(T t, uchar* edgeStart, const uintE& source,
                   const uintT& degree) {
  if (degree > 0) {
    uintE ngh = eatFirstEdge(edgeStart, source);
    W wgh = eatWeight<W>(edgeStart);
    if (!t(source, ngh, wgh, 0)) {
      return;
    }
    for (size_t i = 1; i < degree; i++) {
      ngh += eatEdge(edgeStart);
      wgh = eatWeight<W>(edgeStart);
      if (!t(source, ngh, wgh, i)) {
        break;
      }
    }
  }
}

template <class W>
inline std::tuple<uintE, W> get_ith_neighbor(uchar* edge_start, uintE source,
                                             uintE degree, size_t i) {
  if (degree == 0) {
    assert(false);  // Trying to access edge i = " << i << " from zero deg vtx
  }
  uintE ngh = eatFirstEdge(edge_start, source);
  W wgh = eatWeight<W>(edge_start);
  for (size_t k = 1; k < i; k++) {
    ngh += eatEdge(edge_start);
    wgh = eatWeight<W>(edge_start);
  }
  return std::make_tuple(ngh, wgh);
}

template <class W, class T>
inline void decode_block_seq(T t, uchar* edge_start, const uintE& source,
                             const uintT& degree, uintE block_size,
                             uintE block_num) {
  assert(false);  // Unimplemented
}

template <class W, class M, class Monoid>
inline typename Monoid::T map_reduce(uchar* edge_start, const uintE& source,
                                     const uintT& degree, M& m, Monoid& reduce,
                                     const bool par = true) {
  using E = typename Monoid::T;
  if (degree > 0) {
    uintE ngh = eatFirstEdge(edge_start, source);
    W wgh = eatWeight<W>(edge_start);
    E cur = m(source, ngh, wgh);
    for (size_t i = 1; i < degree; i++) {
      ngh += eatEdge(edge_start);
      wgh = eatWeight<W>(edge_start);
      cur = reduce.f(cur, m(source, ngh, wgh));
    }
    return cur;
  }
  return reduce.identity;
}

/*
  Compresses the first edge, writing target-source and a sign bit.
*/
long compressFirstEdge(uchar* start, long offset, long source, long target);

template <class W, typename std::enable_if<std::is_same<W, gbbs::empty>::value,
                                           int>::type = 0>
inline long compressWeight(uchar* start, long offset, W weight) {
  return offset;
}

template <class W,
          typename std::enable_if<std::is_same<W, intE>::value, int>::type = 0>
inline long compressWeight(uchar* start, long offset, W weight) {
  return compressFirstEdge(start, offset, 0, weight);
}

template <class W,
          typename std::enable_if<std::is_same<W, float>::value, int>::type = 0>
inline long compressWeight(uchar* start, long offset, W weight) {
  float* float_start = (float*)(start + offset);
  *float_start = weight;
  return offset + sizeof(float);
}

inline long compressEdge(uchar* start, long curOffset, uintE diff) {
  uchar curByte = diff & 0x7f;
  int bytesUsed = 0;
  while ((curByte > 0) || (diff > 0)) {
    bytesUsed++;
    uchar toWrite = curByte;
    diff = diff >> 7;
    // Check to see if there's any bits left to represent
    curByte = diff & 0x7f;
    if (diff > 0) {
      toWrite |= 0x80;
    }
    start[curOffset] = toWrite;
    curOffset++;
  }
  return curOffset;
}

template <class W, class P>
inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                   const uintE& degree, std::tuple<uintE, W>* tmp) {
  size_t new_deg = 0;
  if (degree > 0) {
    uchar* finger = edge_start;  // read-finger
    uintE last_ngh;              // the last neighbor written
    size_t offset = 0;           // write offset (from edge_start)
    // Consume the specially encoded first value
    uintE ngh = eatFirstEdge(finger, source);
    W wgh = eatWeight<W>(finger);
    if (pred(source, ngh, wgh)) {
      offset = compressFirstEdge(edge_start, offset, source, ngh);
      offset = compressWeight<W>(edge_start, offset, wgh);
      last_ngh = ngh;
      new_deg++;
    }
    for (size_t i = 1; i < degree; i++) {
      ngh += eatEdge(finger);
      wgh = eatWeight<W>(finger);
      if (pred(source, ngh, wgh)) {
        if (new_deg == 0) {
          offset = compressFirstEdge(edge_start, offset, source, ngh);
          offset = compressWeight<W>(edge_start, offset, wgh);
        } else {
          offset = compressEdge(edge_start, offset, ngh - last_ngh);
          offset = compressWeight<W>(edge_start, offset, wgh);
        }
        last_ngh = ngh;
        new_deg++;
      }
    }
  }
  return new_deg;
}

template <class W, class P, class O>
inline size_t filter(P pred, uchar* edge_start, const uintE& source,
                     const uintE& degree, std::tuple<uintE, W>* tmp, O& out) {
  if (degree > 0) {
    uchar* finger = edge_start;  // read-finger
    // Consume the specially encoded first value
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
  }
  return 0;
}

template <class W, class I>
inline long sequentialCompressEdgeSet(uchar* edgeArray, long currentOffset,
                                      uintE deg, uintE src, I& it) {
  if (deg > 0) {
    // Compress the first edge whole, which is signed difference coded
    uintE prev = std::get<0>(it.cur());
    W w = std::get<1>(it.cur());
    currentOffset = compressFirstEdge(edgeArray, currentOffset, src, prev);
    currentOffset = compressWeight<W>(edgeArray, currentOffset, w);
    for (size_t eid = 1; eid < deg; eid++) {
      // Store difference between cur and prev edge.
      uintE difference = std::get<0>(it.next()) - prev;
      currentOffset = compressEdge(edgeArray, currentOffset, difference);
      currentOffset =
          compressWeight<W>(edgeArray, currentOffset, std::get<1>(it.cur()));
      prev = std::get<0>(it.cur());
    }
  }
  return currentOffset;
}

template <class W>
struct iter {
  uchar* finger;
  uintE src;
  std::tuple<uintE, W> last_edge;
  uintE degree;
  uintE proc;

  iter(uchar* _finger, uintT _degree, uintE _src)
      : finger(_finger), src(_src), degree(_degree), proc(0) {
    if (degree > 0) {
      std::get<0>(last_edge) = eatFirstEdge(finger, src);
      std::get<1>(last_edge) = eatWeight<W>(finger);
      proc = 1;
    }
  }

  inline std::tuple<uintE, W> cur() { return last_edge; }

  inline std::tuple<uintE, W> next() {
    std::get<0>(last_edge) += eatEdge(finger);
    std::get<1>(last_edge) = eatWeight<W>(finger);
    proc++;
    return last_edge;
  }

  inline bool has_next() { return proc < degree; }
};

template <class W>
inline size_t intersect(uchar* l1, uchar* l2, uintE l1_size, uintE l2_size,
                        uintE l1_src, uintE l2_src) {
  if (l1_size == 0 || l2_size == 0) return 0;
  auto it_1 = iter<W>(l1, l1_size, l1_src);
  auto it_2 = iter<W>(l2, l2_size, l2_src);
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
  auto it_1 = iter<W>(l1, l1_size, l1_src);
  auto it_2 = iter<W>(l2, l2_size, l2_src);
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

}  // namespace byte
}  // namespace gbbs
