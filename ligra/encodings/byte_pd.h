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

#include "pbbslib/binary_search.h"
#include "pbbslib/seq.h"
#include "pbbslib/sequence_ops.h"
#include "pbbslib/utilities.h"

#include "ligra/macros.h"

namespace bytepd {

template <class W>
inline std::tuple<uintE, W> get_ith_neighbor(uchar* edge_start, uintE source,
                                             uintE degree, size_t i) {}

// Read default weight (expects pbbslib::empty)
template <class W,
          typename std::enable_if<!std::is_same<W, intE>::value, int>::type = 0>
inline W eatWeight(uchar*& start) {
  return (W)pbbslib::empty();
}

// Read integer weight
template <class W,
          typename std::enable_if<std::is_same<W, intE>::value, int>::type = 0>
inline W eatWeight(uchar*& start) {
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

/*
  Reads any edge of an out-edge list after the first edge.
*/
inline uintE eatEdge(uchar*& start) {
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
    typename std::enable_if<std::is_same<W, pbbslib::empty>::value, int>::type = 0>
inline void decode(T t, uchar* edge_start, const uintE& source,
                   const uintT& degree, const bool par = true) {
  if (degree > 0) {
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)
        edge_start;  // use beginning of edgeArray for offsets into edge list
    uchar* finger = edge_start + (num_blocks - 1) * sizeof(uintE);
    // do first chunk
    size_t end = std::min<long>(PARALLEL_DEGREE, degree);

    // Eat first edge, which is compressed specially
    uintE ngh = eatFirstEdge(finger, source);
    if (!t(source, ngh, pbbslib::empty(), 0)) return;
    for (uintE edgeID = 1; edgeID < end; edgeID++) {
      // Eat the next 'edge', which is a difference, and reconstruct edge.
      ngh += eatEdge(finger);
      if (!t(source, ngh, pbbslib::empty(), edgeID)) return;
    }
    // do remaining chunks in parallel
    par_for(1, num_blocks, 1, [&] (size_t i) {
      size_t o = i * PARALLEL_DEGREE;
      size_t end = std::min<long>(o + PARALLEL_DEGREE, degree);
      uchar* finger = edge_start + block_offsets[i - 1];
      // Eat first edge, which is compressed specially
      uintE ngh = eatFirstEdge(finger, source);
      if (!t(source, ngh, pbbslib::empty(), o)) end = 0;
      for (size_t edgeID = o + 1; edgeID < end; edgeID++) {
        // Eat the next 'edge', which is a difference, and reconstruct edge.
        ngh += eatEdge(finger);
        if (!t(source, ngh, pbbslib::empty(), edgeID)) break;
      }
    }, par);
  }
}

// Decode weighted edges
template <class W, class T,
          typename std::enable_if<!std::is_same<W, pbbslib::empty>::value,
                                  int>::type = 0>
inline void decode(T t, uchar* edge_start, const uintE& source,
                   const uintT& degree, const bool par = true) {
  if (degree > 0) {
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)
        edge_start;  // use beginning of edgeArray for offsets into edge list
    uchar* finger = edge_start + (num_blocks - 1) * sizeof(uintE);
    // do first chunk
    size_t end = std::min<long>(PARALLEL_DEGREE, degree);

    // Eat first edge, which is compressed specially
    uintE ngh = eatFirstEdge(finger, source);
    W weight = eatWeight<W>(finger);
    if (!t(source, ngh, weight, 0)) return;
    for (size_t edgeID = 1; edgeID < end; edgeID++) {
      // Eat the next 'edge', which is a difference, and reconstruct edge.
      ngh += eatEdge(finger);
      W weight = eatWeight<W>(finger);
      if (!t(source, ngh, weight, edgeID)) return;
    }
    // do remaining chunks in parallel
    par_for(1, num_blocks, 1, [&] (size_t i) {
      size_t o = i * PARALLEL_DEGREE;
      size_t end = std::min<long>(o + PARALLEL_DEGREE, degree);
      uchar* finger = edge_start + block_offsets[i - 1];
      // Eat first edge, which is compressed specially
      uintE ngh = eatFirstEdge(finger, source);
      W weight = eatWeight<W>(finger);
      if (!t(source, ngh, weight, o)) end = 0;
      for (size_t edgeID = o + 1; edgeID < end; edgeID++) {
        // Eat the next 'edge', which is a difference, and reconstruct edge.
        ngh += eatEdge(finger);
        W weight = eatWeight<W>(finger);
        if (!t(source, ngh, weight, edgeID)) break;
      }
    }, par);
  }
}

// Map_reduce
template <class W, class E, class M, class R>
inline E map_reduce(uchar* edge_start, const uintE& source, const uintT& degree,
                    E id, M& m, R& r, const bool par = true) {
  if (degree > 0) {
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)edge_start;

    E stk[1000];
    E* block_outputs;
    if (num_blocks > 1000) {
      block_outputs = pbbslib::new_array_no_init<E>(num_blocks);
    } else {
      block_outputs = (E*)stk;
    }

    par_for(0, num_blocks, 2, [&] (size_t i) {
      size_t start = i * PARALLEL_DEGREE;
      size_t end = std::min<long>(start + PARALLEL_DEGREE, degree);
      uchar* finger = edge_start + ((i == 0) ? (num_blocks - 1) * sizeof(uintE)
                                             : block_offsets[i - 1]);

      // Eat first edge, which is compressed specially
      uintE ngh = eatFirstEdge(finger, source);
      W wgh = eatWeight<W>(finger);
      E cur = m(source, ngh, wgh);
      for (size_t j = start + 1; j < end; j++) {
        ngh += eatEdge(finger);
        W wgh = eatWeight<W>(finger);
        cur = r(cur, m(source, ngh, wgh));
      }
      block_outputs[i] = cur;
    }, par);

    auto im = pbbslib::make_sequence(block_outputs, num_blocks);
    E res = pbbslib::reduce(im, r);
    if (num_blocks > 1000) {
      pbbslib::free_array(block_outputs);
    }
    return res;
  } else {
    return id;
  }
}

template <class W>
uchar* decode_block(uchar* finger, std::tuple<uintE, W>* out, size_t start,
                    size_t end, const uintE& source) {
  uintE ngh = eatFirstEdge(finger, source);
  W wgh = eatWeight<W>(finger);
  out[start] = std::make_tuple(ngh, wgh);
  for (size_t i = start + 1; i < end; i++) {
    // Eat the next 'edge', which is a difference, and reconstruct edge.
    ngh += eatEdge(finger);
    W wgh = eatWeight<W>(finger);
    out[i] = std::make_tuple(ngh, wgh);
  }
  return finger;
}

template <class W, class T>
inline void decode_block_seq(T t, uchar* edge_start, const uintE& source,
                             const uintT& degree, uintE block_size,
                             uintE block_num) {
  assert(false);
}

inline size_t get_virtual_degree(uintE d, uchar* ngh_arr) { return d; }

#define SEQ_THRESH 10

// Represents the sequence from [...[|cur_block|...|cur_block+num_blocks|]...]
struct seq_info {
  uchar* edge_start;
  uintE degree;
  uintE source_id;
  uintE total_blocks;

  uintE start;
  uintE end;

  seq_info(uchar* es, uintE d, uintE sid, uintE tb, uintE s, uintE e)
      : edge_start(es),
        degree(d),
        source_id(sid),
        total_blocks(tb),
        start(s),
        end(e) {}

  uchar* get_start_of_block(const uintE& block_id) {
    if (total_blocks == 1) {
      return edge_start;
    }
    uintE* offs = (uintE*)edge_start;
    return edge_start + offs[block_id - 1];
  }

  uintE get_pivot() {
    uintE* offs = (uintE*)edge_start;
    uintE mid_block = (end + start) / 2;
    uchar* finger = edge_start + offs[mid_block];
    return eatFirstEdge(finger, source_id);
  }

  uintE pivot_block() { return (end + start) / 2; }

  uintE binary_search(uintE pivot) {
    uintE* offs = (uintE*)edge_start;
    auto start_f = [&](size_t i) {
      uchar* finger = edge_start + offs[start + i];
      return eatFirstEdge(finger, source_id);
    };
    auto start_im = pbbslib::make_sequence<uintE>(size(), start_f);
    uintE ind =
        pbbslib::binary_search(start_im, pivot, std::greater<uintE>());  // check
    // ind is the first block index (from start) <= our pivot.
    //    decode_block<

    return ind;
  }

  seq_info cut(uintE l, uintE r) {
    return seq_info(edge_start, degree, source_id, total_blocks, l, r);
  }

  uintE size() { return end - start; }
};

inline uintE seq_intersect_full(seq_info u, seq_info v) {
  assert(false);  // not implemented
  return 0;
}

inline uintE seq_intersect(seq_info u, seq_info v) {
  assert(false);  // not implemented
  return 0;
  //  decode_block<pbbslib::empty>(finger, (std::tuple<uintE, pbbslib::empty>*)ngh_u,
  //  0,
}

inline uintE intersect(seq_info u, seq_info v) {
  // Might need to swap here
  uintE nA = u.size();
  uintE nB = v.size();
  if (nA + nB < SEQ_THRESH) {
    return seq_intersect_full(u, v);
  } else if (nA == 1) {  // merge base
    return seq_intersect(u, v);
  } else if (nB == 0) {
    return 0;
  } else {  // (large, large)
    uintE pivot = u.get_pivot();
    uintE mU = u.pivot_block();
    uintE mV = v.binary_search(pivot);
    uintE lA = 0, rA = 0;
    par_do(true, [&]() { lA = intersect(u.cut(0, mU), v.cut(0, mV)); },
           [&]() { rA = intersect(u.cut(mU, mU), v.cut(0, mV)); });
    return lA + rA;
  }
}

/*
  Compresses the first edge, writing target-source and a sign bit.
*/
long compressFirstEdge(uchar* start, long offset, uintE source, uintE target) {
  intE preCompress = (intE)target - source;
  int bytesUsed = 0;
  uchar firstByte = 0;
  intE toCompress = std::abs(preCompress);
  firstByte = toCompress & 0x3f;  // 0011|1111
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

long compressEdge(uchar* start, long curOffset, uintE e) {
  uchar curByte = e & 0x7f;
  int bytesUsed = 0;
  while ((curByte > 0) || (e > 0)) {
    bytesUsed++;
    uchar toWrite = curByte;
    e = e >> 7;
    // Check to see if there's any bits left to represent
    curByte = e & 0x7f;
    if (e > 0) {
      toWrite |= 0x80;
    }
    start[curOffset] = toWrite;
    curOffset++;
  }
  return curOffset;
}

template <class W>
size_t compute_block_size(std::tuple<uintE, W>* edges, size_t start, size_t end,
                          const uintE& source) {
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
size_t compute_size_in_bytes(std::tuple<uintE, W>* edges, const uintE& source,
                             const uintE& d) {
  if (d > 0) {
    size_t num_blocks = 1 + (d - 1) / PARALLEL_DEGREE;
    uintE stk[100];
    uintE* block_bytes;
    if (num_blocks > 100) {
      block_bytes = pbbslib::new_array_no_init<uintE>(num_blocks);
    } else {
      block_bytes = (uintE*)stk;
    }
    par_for(0, num_blocks, 10, [&] (size_t i) {
      // calculate size in bytes of this block
      uintE start = i * PARALLEL_DEGREE;
      uintE end = start + std::min<uintE>(PARALLEL_DEGREE, d - start);
      block_bytes[i] = compute_block_size(edges, start, end, source);
    });
    auto bytes_imap = pbbslib::make_sequence(block_bytes, num_blocks);
    size_t total_space = pbbslib::scan_add_inplace(bytes_imap);

    // add in space for storing offsets to the start of each block
    total_space += sizeof(uintE) * (num_blocks - 1);
    if (num_blocks > 100) {
      pbbslib::free_array(block_bytes);
    }
    return total_space;
  } else {
    return 0;
  }
}

template <class W, class I>
long sequentialCompressEdgeSet(uchar* edgeArray, size_t current_offset,
                               uintT degree, uintE source, I& it) {
  if (degree > 0) {
    size_t start_offset = current_offset;
    size_t num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)edgeArray;
    current_offset += (num_blocks - 1) * sizeof(uintE);  // block_offs
    for (long i = 0; i < num_blocks; i++) {
      size_t o = i * PARALLEL_DEGREE;
      size_t end = std::min<size_t>(PARALLEL_DEGREE, degree - o);

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

/*
  Takes:
    1. The edge array of chars to write into
    2. The current offset into this array
    3. The vertices degree
    4. The vertices vertex number
    5. The array of saved out-edges we're compressing
  Returns:
    The new offset into the edge array
*/
long sequentialCompressEdgeSet(uchar* edgeArray, size_t current_offset,
                               uintT degree, uintE source, uintE* edges) {
  if (degree > 0) {
    long start_offset = current_offset;
    long num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)
        edgeArray;  // use beginning of edgeArray for offsets into edge list
    current_offset += (num_blocks - 1) * sizeof(uintE);
    for (long i = 0; i < num_blocks; i++) {
      long o = i * PARALLEL_DEGREE;
      long end = std::min<long>(PARALLEL_DEGREE, degree - o);
      uintE* edges_i = edges + o;
      if (i > 0)
        block_offsets[i - 1] =
            current_offset -
            start_offset;  // store offset for all chunks but the first
      // Compress the first edge whole, which is signed difference coded
      current_offset =
          compressFirstEdge(edgeArray, current_offset, source, edges_i[0]);
      for (uintT edgeI = 1; edgeI < end; edgeI++) {
        // Store difference between cur and prev edge.
        uintE difference = edges_i[edgeI] - edges_i[edgeI - 1];
        current_offset = compressEdge(edgeArray, current_offset, difference);
      }
    }
  }
  return current_offset;
}

// Packs out the block sequentially.
template <class W, class P>
size_t packBlock(uchar* edge_start, const uintE& source, const uintE& degree,
                 P& pred) {
  uchar* finger = edge_start;

  //  auto map_f = [&] (const uintE& src, const uintE& ngh, const W& w) {
  //    assert(src < 3072627);
  //    assert(ngh < 3072627);
  //    return 1;
  //  };
  //  auto red_f = [&] (size_t l, size_t r) { return l + r; };
  //  map_reduce<W>(edge_start, source, degree, 0, map_f, red_f, true);
  //
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
      W wgh = eatWeight<W>(finger);
      if (pred(source, ngh, wgh)) {
        out(k++, std::make_tuple(ngh, wgh));
      }
    }
  }
  return k;
}

template <class W>
void compress_edges(uchar* edgeArray, const uintE& source, const uintE& d,
                    std::tuple<uintE, W>* edges, uchar* last_finger,
                    bool par = true) {
  size_t num_blocks = 1 + (d - 1) / PARALLEL_DEGREE;
  uintE stk[101];
  uintE* block_bytes;
  if (num_blocks > 100) {
    block_bytes = pbbslib::new_array_no_init<uintE>(num_blocks + 1);
  } else {
    block_bytes = (uintE*)stk;
  }

  par_for(0, num_blocks, 1, [&] (size_t i) {
    uintE start = i * PARALLEL_DEGREE;
    uintE end = std::min<uintE>(start + PARALLEL_DEGREE, d);
    block_bytes[i] = compute_block_size(edges, start, end, source);
  });
  block_bytes[num_blocks] = 0;

  uintE edges_offset = (num_blocks - 1) * sizeof(uintE);
  uintE* block_offsets = (uintE*)edgeArray;

  auto bytes_imap = pbbslib::make_sequence(block_bytes, num_blocks + 1);
  uintE total_space = pbbslib::scan_add_inplace(bytes_imap);

  if (total_space > (last_finger - edgeArray)) {
    std::cout << "# Space error!"
              << "\n";
    exit(0);
  }
  if (total_space == (last_finger - edgeArray)) {
    std::cout << "# d = " << d << " to exactly the same space"
              << "\n";
  }

  par_for(1, num_blocks, 1, [&] (size_t i) {
    // store offset to start of this block from edgeArray
    block_offsets[i - 1] = edges_offset + bytes_imap[i];
  });

  par_for(0, num_blocks, 1, [&] (size_t i) {
    uintE start = i * PARALLEL_DEGREE;
    uintE end = std::min<uintE>(start + PARALLEL_DEGREE, d);
    long write_offset = (i == 0) ? edges_offset : block_offsets[i - 1];
    auto nw = edges[start];
    write_offset =
        compressFirstEdge(edgeArray, write_offset, source, std::get<0>(nw));
    write_offset = compressWeight<W>(edgeArray, write_offset, std::get<1>(nw));
    uintE last_edge = std::get<0>(nw);
    for (size_t i = start + 1; i < end; i++) {
      nw = edges[i];
      write_offset =
          compressEdge(edgeArray, write_offset, std::get<0>(nw) - last_edge);
      write_offset =
          compressWeight<W>(edgeArray, write_offset, std::get<1>(nw));
      last_edge = std::get<0>(nw);
    }
  });
  if (num_blocks > 100) {
    pbbslib::free_array(block_bytes);
  }
}

template <class W, class P>
inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                   const uintE& degree, std::tuple<uintE, W>* tmp) {
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

    //    auto map_f = [&] (const uintE& src, const uintE& ngh, const W& w) {
    //      assert(src < 3072627);
    //      assert(ngh < 3072627);
    //      return 1;
    //    };
    //    auto red_f = [&] (size_t l, size_t r) { return l + r; };
    //    size_t ct = map_reduce<W>(edge_start, source, degree, 0, map_f, red_f,
    //    true);

    uchar* mf = 0;
    par_for(0, num_blocks, 1, [&] (size_t i) {
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
    uintE k = pbbslib::filterf(our_tmp, tmp2, degree, pd);
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
      W wgh = eatWeight<W>(finger);
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

    par_for(0, num_blocks, 1, [&] (size_t i) {
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
    uintE k = pbbslib::filterf(tmp, tmp2, degree, pd);
    par_for(0, k, 2000, [&] (size_t i) { out(i, tmp2[i]); });
    return k;
  }
  return 0;
}

/*
  Compresses the edge set in parallel.
*/
uintE* parallelCompressEdges(uintE* edges, uintT* offsets, long n, long m,
                             uintE* Degrees) {
  std::cout << "parallel compressing, (n,m) = (" << n << "," << m << ")"
            << "\n";
  uintE** edgePts = pbbslib::new_array_no_init<uintE*>(n);
  long* charsUsedArr = pbbslib::new_array_no_init<long>(n+1);
  auto charsUsed = pbbslib::make_sequence(charsUsedArr, n+1);
//  long* compressionStarts = pbbslib::new_array_no_init<long>(n + 1);
  {
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { charsUsed[i] = ceil((Degrees[i] * 9) / 8) + 4; });
  }
  charsUsed[n] = 0;
  long toAlloc = pbbslib::scan_inplace(charsUsed, pbbslib::addm<long>());
  uintE* iEdges = pbbslib::new_array_no_init<uintE>(toAlloc);

  {
    par_for(0, n, [&] (size_t i) {
      edgePts[i] = iEdges + charsUsed[i];
      long charsUsed =
          sequentialCompressEdgeSet((uchar*)(iEdges + charsUsed[i]), 0,
                                    Degrees[i], i, edges + offsets[i]);
      charsUsed[i] = charsUsed;
    });
  }

  // produce the total space needed for all compressed lists in chars.
  auto compressionStarts = pbbslib::scan(charsUsedArr, pbbslib::addm<long>());
  size_t totalSpace = compressionStarts[n];
  pbbslib::free_array(charsUsedArr);

  uchar* finalArr = pbbslib::new_array_no_init<uchar>(totalSpace);
  std::cout << "total space requested is : " << totalSpace << "\n";
  float avgBitsPerEdge = (float)totalSpace * 8 / (float)m;
  std::cout << "Average bits per edge: " << avgBitsPerEdge << "\n";

  {
    par_for(0, n, [&] (size_t i) {
      long o = compressionStarts[i];
      memcpy(finalArr + o, (uchar*)(edgePts[i]), compressionStarts[i + 1] - o);
      offsets[i] = o;
    });
  }
  offsets[n] = totalSpace;
  pbbslib::free_array(iEdges);
  pbbslib::free_array(edgePts);
  std::cout << "finished compressing, bytes used = " << totalSpace << "\n";
  std::cout << "would have been, " << (m * 4) << "\n";
  return ((uintE*)finalArr);
}

long sequentialCompressWeightedEdgeSet(uchar* edgeArray, long current_offset,
                                       const uintE& degree, const uintE& source,
                                       std::tuple<uintE, intE>* edges) {
  if (degree > 0) {
    uchar* our_start = edgeArray + current_offset;
    long start_offset = current_offset;
    long num_blocks = 1 + (degree - 1) / PARALLEL_DEGREE;
    uintE* block_offsets = (uintE*)our_start;
    current_offset += (num_blocks - 1) * sizeof(uintE);

    for (long i = 0; i < num_blocks; i++) {
      long start = i * PARALLEL_DEGREE;
      long end = start + std::min<long>(PARALLEL_DEGREE, degree - start);
      if (i > 0) block_offsets[i - 1] = current_offset - start_offset;

      auto nw = edges[start];
      current_offset =
          compressFirstEdge(edgeArray, current_offset, source, std::get<0>(nw));
      current_offset =
          compressWeight<intE>(edgeArray, current_offset, std::get<1>(nw));

      uintE last_ngh = std::get<0>(nw);
      for (size_t i = start + 1; i < end; i++) {
        auto nw = edges[i];
        uintE difference = std::get<0>(nw) - last_ngh;
        current_offset = compressEdge(edgeArray, current_offset, difference);
        current_offset =
            compressWeight<intE>(edgeArray, current_offset, std::get<1>(nw));
        last_ngh = std::get<0>(nw);
      }
    }
  }
  return current_offset;
}

uchar* parallelCompressWeightedEdges(std::tuple<uintE, intE>* edges,
                                     uintT* offsets, long n, long m,
                                     uintE* Degrees) {
  std::cout << "parallel compressing, (n,m) = (" << n << "," << m << ")"
            << "\n";
  auto bytes_used = pbbslib::new_array_no_init<size_t>(n + 1);
  auto bytes_seq = pbbslib::make_sequence<size_t>(bytes_used, n+1);

  par_for(0, n, [&] (size_t i) {
    bytes_used[i] = compute_size_in_bytes(edges + offsets[i], i, Degrees[i]);
  });
  bytes_used[n] = 0;
  size_t total_bytes = pbbslib::scan_inplace(bytes_used, pbbslib::addm<size_t>());

  uchar* edges_c = pbbslib::new_array_no_init<uchar>(total_bytes);

  par_for(0, n, [&] (size_t i) {
    sequentialCompressWeightedEdgeSet(edges_c, bytes_used[i], Degrees[i], i,
                                      edges + offsets[i]);
  });

  float avgBitsPerEdge = (float)total_bytes * 8 / (float)m;
  std::cout << "Average bits per edge: " << avgBitsPerEdge << "\n";

  size_t end = n+1;
  par_for(0, end, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { offsets[i] = bytes_used[i]; });
  std::cout << "finished compressing, bytes used = " << total_bytes << "\n";
  std::cout << "would have been, " << (m * 8) << "\n";
  return edges_c;
}

};  // namespace bytepd
