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

#include "byte.h"
#include "byte_pd.h"
#include "byte_pd_amortized.h"

namespace gbbs {

struct byte_decode {
  template <class W, class T>
  static inline void decode_block(T t, uchar* edge_start, const uintE& source,
                                  const uintT& degree, uintE block_num) {
    assert(false);  // Unimplemented for byte decoder
    exit(0);
    bytepd_amortized::decode_block<W, T>(t, edge_start, source, degree,
                                         block_num);
  }

  static inline size_t get_virtual_degree(uintE d, uchar* nghArr) { return d; }

  template <class W>
  static inline auto iter(uchar* edge_start, uintE degree, uintE id)
      -> byte::iter<W> {
    return byte::iter<W>(edge_start, degree, id);
  }

  template <class W>
  static inline size_t intersect(uchar* l1, uchar* l2, uintE l1_size,
                                 uintE l2_size, uintE l1_src, uintE l2_src) {
    return byte::intersect<W>(l1, l2, l1_size, l2_size, l1_src, l2_src);
  }

  template <class W, class F>
  static inline size_t intersect_f(uchar* l1, uchar* l2, uintE l1_size,
                                   uintE l2_size, uintE l1_src, uintE l2_src,
                                   const F& f) {
    return byte::intersect_f<W>(l1, l2, l1_size, l2_size, l1_src, l2_src, f);
  }

  template <class W, class I>
  static inline long sequentialCompressEdgeSet(uchar* edgeArray,
                                               size_t current_offset,
                                               uintT degree, uintE source,
                                               I& it) {
    return byte::sequentialCompressEdgeSet(edgeArray, current_offset, degree,
                                           source, it);
  }

  template <class W, class P, class O>
  static inline void filter(P pred, uchar* edge_start, const uintE& source,
                            const uintE& degree, std::tuple<uintE, W>* tmp,
                            O& out) {
    byte::filter(pred, edge_start, source, degree, tmp, out);
    return;
  }

  template <class W, class P>
  static inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                            const uintE& degree,
                            std::tuple<uintE, W>* tmp_space, bool par = true) {
    return byte::pack(pred, edge_start, source, degree, tmp_space);
  }

  template <class W, class M, class Monoid>
  static inline typename Monoid::T map_reduce(uchar* edge_start,
                                              const uintE& source,
                                              const uintT& degree, M& m,
                                              Monoid& reduce,
                                              const bool par = true) {
    return byte::map_reduce<W, M, Monoid>(edge_start, source, degree, m, reduce,
                                          par);
  }

  template <class W, class T>
  __attribute__((always_inline)) static inline void decode(
      T t, uchar* edge_start, const uintE& source, const uintT& degree,
      const bool& parallel) {
    return byte::decode<W, T>(t, edge_start, source, degree);
  }

  template <class W>
  static inline std::tuple<uintE, W> get_ith_neighbor(uchar* edge_start,
                                                      uintE source,
                                                      uintE degree, size_t i) {
    return byte::get_ith_neighbor<W>(edge_start, source, degree, i);
  }

  static inline uintE get_num_blocks(uchar* edge_start, uintE degree) {
    return 1;  // single block in byte-compressed format
  }
};

struct bytepd_decode {
  template <class W>
  static inline size_t intersect(uchar* l1, uchar* l2, uintE l1_size,
                                 uintE l2_size, uintE l1_src, uintE l2_src) {
    exit(-1);  // TODO
    return 1;
  }

  template <class W, class F>
  static inline size_t intersect_f(uchar* l1, uchar* l2, uintE l1_size,
                                   uintE l2_size, uintE l1_src, uintE l2_src,
                                   const F& f) {
    exit(-1);  // TODO
    return 1;
  }

  template <class W>
  static inline auto iter(uchar* edge_start, uintE degree, uintE id)
      -> bytepd::iter<W> {
    return bytepd::iter<W>(edge_start, degree, id);
  }

  template <class W, class I>
  static inline long sequentialCompressEdgeSet(uchar* edgeArray,
                                               size_t current_offset,
                                               uintT degree, uintE source,
                                               I& it) {
    return bytepd::sequentialCompressEdgeSet(edgeArray, current_offset, degree,
                                             source, it);
  }

  template <class W, class P, class O>
  static inline void filter(P pred, uchar* edge_start, const uintE& source,
                            const uintE& degree, std::tuple<uintE, W>* tmp,
                            O& out) {
    return bytepd::filter(pred, edge_start, source, degree, tmp, out);
  }

  template <class W, class P>
  static inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                            const uintE& degree,
                            std::tuple<uintE, W>* tmp_space, bool par = true) {
    return bytepd::pack(pred, edge_start, source, degree, tmp_space, par);
  }

  template <class W, class M, class Monoid>
  static inline typename Monoid::T map_reduce(uchar* edge_start,
                                              const uintE& source,
                                              const uintT& degree, M& m,
                                              Monoid& reduce,
                                              const bool par = true) {
    return bytepd::map_reduce<W>(edge_start, source, degree, m, reduce, par);
  }

  template <class W, class T>
  __attribute__((always_inline)) static inline void decode(
      T& t, uchar* edge_start, const uintE& source, const uintT& degree,
      const bool parallel = true) {
    return bytepd::decode<W, T>(t, edge_start, source, degree, parallel);
  }

  static inline size_t get_virtual_degree(uintE d, uchar* nghArr) {
    return bytepd::get_virtual_degree(d, nghArr);
  }

  template <class W, class T>
  static inline void decode_block(T t, uchar* edge_start, const uintE& source,
                                  const uintT& degree, uintE block_num) {
    bytepd::decode_block<W, T>(t, edge_start, source, degree, block_num);
  }

  template <class W>
  static inline std::tuple<uintE, W> get_ith_neighbor(uchar* edge_start,
                                                      uintE source,
                                                      uintE degree, size_t i) {
    return bytepd::get_ith_neighbor<W>(edge_start, source, degree, i);
  }

  static inline uintE get_num_blocks(uchar* edge_start, uintE degree) {
    return bytepd::get_num_blocks(edge_start, degree);
  }

  static inline uintE get_block_degree(uchar* edge_start, uintE degree,
                                       uintE block_num) {
    return bytepd::get_block_degree(edge_start, degree, block_num);
  }
};

struct bytepd_amortized_decode {
  template <class W>
  static inline size_t intersect(uchar* l1, uchar* l2, uintE l1_size,
                                 uintE l2_size, uintE l1_src, uintE l2_src) {
    return bytepd_amortized::intersect<W>(l1, l2, l1_size, l2_size, l1_src,
                                          l2_src);
  }

  template <class W, class F>
  static inline size_t intersect_f(uchar* l1, uchar* l2, uintE l1_size,
                                   uintE l2_size, uintE l1_src, uintE l2_src,
                                   const F& f) {
    return bytepd_amortized::intersect_f<W>(l1, l2, l1_size, l2_size, l1_src,
                                            l2_src, f);
  }

  template <class W>
  static inline auto iter(uchar* edge_start, uintE degree, uintE id)
      -> bytepd_amortized::iter<W> {
    return bytepd_amortized::iter<W>(edge_start, degree, id);
  }

  template <class W, class I>
  static inline long sequentialCompressEdgeSet(uchar* edgeArray,
                                               size_t current_offset,
                                               uintT degree, uintE source,
                                               I& it) {
    return bytepd_amortized::sequentialCompressEdgeSet(
        edgeArray, current_offset, degree, source, it);
  }

  template <class W, class P, class O>
  static inline void filter(P pred, uchar* edge_start, const uintE& source,
                            const uintE& degree, std::tuple<uintE, W>* tmp,
                            O& out) {
    return bytepd_amortized::filter(pred, edge_start, source, degree, tmp, out);
  }

  template <class W, class P>
  static inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                            const uintE& degree,
                            std::tuple<uintE, W>* tmp_space, bool par = true) {
    return bytepd_amortized::pack(pred, edge_start, source, degree, tmp_space,
                                  par);
  }

  template <class W, class M, class Monoid>
  static inline typename Monoid::T map_reduce(uchar* edge_start,
                                              const uintE& source,
                                              const uintT& degree, M& m,
                                              Monoid& reduce,
                                              const bool par = true) {
    return bytepd_amortized::map_reduce<W>(edge_start, source, degree, m,
                                           reduce, par);
  }

  template <class W, class T>
  __attribute__((always_inline)) static inline void decode(
      T& t, uchar* edge_start, const uintE& source, const uintT& degree,
      const bool parallel = true) {
    return bytepd_amortized::decode<W, T>(t, edge_start, source, degree,
                                          parallel);
  }

  static inline size_t get_virtual_degree(uintE d, uchar* nghArr) {
    return bytepd_amortized::get_virtual_degree(d, nghArr);
  }

  template <class W, class T>
  static inline void decode_block(T t, uchar* edge_start, const uintE& source,
                                  const uintT& degree, uintE block_num) {
    bytepd_amortized::decode_block<W, T>(t, edge_start, source, degree,
                                         block_num);
  }

  template <class W>
  static inline std::tuple<uintE, W> get_ith_neighbor(uchar* edge_start,
                                                      uintE source,
                                                      uintE degree, size_t i) {
    return bytepd_amortized::get_ith_neighbor<W>(edge_start, source, degree, i);
  }

  static inline uintE get_num_blocks(uchar* edge_start, uintE degree) {
    return bytepd_amortized::get_num_blocks(edge_start, degree);
  }

  static inline uintE get_block_degree(uchar* edge_start, uintE degree,
                                       uintE block_num) {
    return bytepd_amortized::get_block_degree(edge_start, degree, block_num);
  }
};

}  // namespace gbbs
