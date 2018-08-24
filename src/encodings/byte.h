// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

#pragma once

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <type_traits>

#include "../macros.h"
#include "../oldlib/utils.h"

namespace byte {

  inline size_t get_virtual_degree(uintE d, uchar* nghArr) {
    return d;
  }

  template <class W, typename std::enable_if<!std::is_same<W, intE>::value, int>::type=0>
  inline W eatWeight(uchar* &start) {
    return (W)pbbs::empty();
  }

  template <class W, typename std::enable_if<std::is_same<W, intE>::value, int>::type=0>
  inline W eatWeight(uchar* &start) {
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

  inline intE eatFirstEdge(uchar* &start, uintE source) {
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

  /*
    Reads any edge of an out-edge list after the first edge.
  */
  inline uintE eatEdge(uchar* &start) {
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
    inline void decode(T t, uchar* edgeStart, const uintE &source, const uintT &degree, const bool par=true) {
    if (degree > 0) {
      uintE ngh = eatFirstEdge(edgeStart,source);
      W wgh = eatWeight<W>(edgeStart);
      if (!t(source, ngh, wgh, 0)) {
        return;
      }
      for (size_t i=1; i < degree; i++) {
        ngh += eatEdge(edgeStart);
        wgh = eatWeight<W>(edgeStart);
        if (!t(source, ngh, wgh, i)) {
          break;
        }
      }
    }
  }

  template <class W, class T>
  inline void decode_block_seq(T t, uchar* edge_start, const uintE &source,
                               const uintT &degree, uintE block_size,
                               uintE block_num) {
    assert(false);
  }

  template <class W, class E, class M, class R>
  inline auto map_reduce(uchar* edge_start, const uintE &source,
      const uintT &degree, E id, M& m, R& r, const bool par=true) {
    if (degree > 0) {
      uintE ngh = eatFirstEdge(edge_start, source);
      W wgh = eatWeight<W>(edge_start);
      E cur = m(source, ngh, wgh);
      for (size_t i=1; i < degree; i++) {
        ngh += eatEdge(edge_start);
        wgh = eatWeight<W>(edge_start);
        cur = r(cur, m(source, ngh, wgh));
      }
      return cur;
    }
    return id;
  }

  /*
    Compresses the first edge, writing target-source and a sign bit.
  */
  long compressFirstEdge(uchar *start, long offset, long source, long target) {
    uchar* saveStart = start;
    long saveOffset = offset;

    long diff = target - source;
    long preCompress = diff;
    int bytesUsed = 0;
    uchar firstByte = 0;
    uintE toCompress = abs(preCompress);
    firstByte = toCompress & 0x3f; // 0011|1111
    if (preCompress < 0) {
      firstByte |= 0x40;
    }
    toCompress = toCompress >> 6;
    if (toCompress > 0) {
      firstByte |= 0x80;
    }
    start[offset] = firstByte;
    offset++;

    uchar curByte = toCompress & 0x7f;
    while ((curByte > 0) || (toCompress > 0)) {
      bytesUsed++;
      uchar toWrite = curByte;
      toCompress = toCompress >> 7;
      // Check to see if there's any bits left to represent
      curByte = toCompress & 0x7f;
      if (toCompress > 0) {
        toWrite |= 0x80;
      }
      start[offset] = toWrite;
      offset++;
    }
    return offset;
  }


  template <class W, typename std::enable_if<!std::is_same<W, intE>::value, int>::type=0>
  inline long compressWeight(uchar* start, long offset, W weight) {
    return offset;
  }

  template <class W, typename std::enable_if<std::is_same<W, intE>::value, int>::type=0>
  inline long compressWeight(uchar* start, long offset, W weight) {
    return compressFirstEdge(start, offset, 0, weight);
  }

  long compressEdge(uchar *start, long curOffset, uintE diff) {
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
  inline size_t pack(P pred, uchar* edge_start, const uintE &source, const uintE &degree, tuple<uintE, W>* tmp) {
    size_t new_deg = 0;
    if (degree > 0) {
      uchar* finger = edge_start; // read-finger
      uintE last_ngh; // the last neighbor written
      size_t offset = 0; // write offset (from edge_start)
      // Consume the specially encoded first value
      uintE ngh = eatFirstEdge(finger, source);
      W wgh = eatWeight<W>(finger);
      if (pred(source, ngh, wgh)) {
        offset = compressFirstEdge(edge_start, offset, source, ngh);
        offset = compressWeight<W>(edge_start, offset, wgh);
        last_ngh = ngh;
        new_deg++;
      }
      for (size_t i=1; i<degree; i++) {
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
  inline size_t filter(P pred, uchar* edge_start, const uintE &source, const uintE &degree, tuple<uintE, W>* tmp, O& out) {
    if (degree > 0) {
      uchar* finger = edge_start; // read-finger
      // Consume the specially encoded first value
      uintE ngh = eatFirstEdge(finger, source);
      W wgh = eatWeight<W>(finger);
      uintE k = 0;
      if (pred(source, ngh, wgh)) {
        out(k++, make_tuple(ngh, wgh));
      }
      for (size_t i=1; i<degree; i++) {
        ngh += eatEdge(finger);
        wgh = eatWeight<W>(finger);
        if (pred(source, ngh, wgh)) {
          out(k++, make_tuple(ngh, wgh));
        }
      }
      return k;
    }
    return 0;
  }

  template <class W, class I>
  long sequentialCompressEdgeSet(uchar *edgeArray, long currentOffset, uintE deg, uintE src, I& it) {
    if (deg > 0) {
      // Compress the first edge whole, which is signed difference coded
      uintE prev = get<0>(it.cur()); W w = get<1>(it.cur());
      long old_off = currentOffset;
      currentOffset = compressFirstEdge(edgeArray, currentOffset, src, prev);
      currentOffset = compressWeight<W>(edgeArray, currentOffset, w);
      for (size_t eid = 1; eid < deg; eid++) {
        // Store difference between cur and prev edge.
        uintE difference = get<0>(it.next()) - prev;
        currentOffset = compressEdge(edgeArray, currentOffset, difference);
        currentOffset = compressWeight<W>(edgeArray, currentOffset, get<1>(it.cur()));
        prev = get<0>(it.cur());
      }
    }
    return currentOffset;
  }

  template <class W>
  struct iter {
    uchar* finger;
    uintE src;
    tuple<uintE, W> last_edge;
    uintE degree;
    uintE proc;

    iter(uchar* _finger, uintT _degree, uintE _src) : finger(_finger), degree(_degree), src(_src), proc(0) {
      if (degree > 0) {
        get<0>(last_edge) = eatFirstEdge(finger, src);
        get<1>(last_edge) = eatWeight<W>(finger);
        proc = 1;
      }
    }

    inline tuple<uintE, W> cur() {
      return last_edge;
    }

    inline tuple<uintE, W> next() {
      get<0>(last_edge) += eatEdge(finger);
      get<1>(last_edge) = eatWeight<W>(finger);
      proc++;
      return last_edge;
    }

    inline bool has_next() {
      return proc < degree;
    }
  };

  template <class W>
  size_t intersect(uchar* l1, uchar* l2, uintE l1_size, uintE l2_size, uintE l1_src, uintE l2_src) {
    if (l1_size == 0 || l2_size == 0) return 0;
    auto it_1 = iter<W>(l1, l1_size, l1_src);
    auto it_2 = iter<W>(l2, l2_size, l2_src);
    size_t i=0, j = 0, ct = 0;
    while (i < l1_size && j < l2_size) {
      uintE e1 = get<0>(it_1.cur());
      uintE e2 = get<0>(it_2.cur());
      if (e1 == e2) { i++, j++, it_1.next(), it_2.next(), ct++; }
      else if (e1 < e2) { i++, it_1.next(); }
      else { j++, it_2.next(); }
    }
    return ct;
  }


}; // namespace byte
