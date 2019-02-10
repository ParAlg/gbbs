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
//#include "byte_pd.h"
#include "byte_pd_amortized.h"

template <class W>
struct byte_decode {
  using iter_type = encodings::byte::iter<W>;

  template <class T>
  static inline void decode_block_seq(T t, uchar* edge_start,
                                      const uintE& source, const uintT& degree,
                                      uintE block_size, uintE block_num) {
    assert(false);  // Unimplemented for byte decoder
    return encodings::bytepd_amortized::decode_block_seq<W, T>(
        t, edge_start, source, degree, block_size, block_num);
  }

  static inline size_t get_virtual_degree(uintE d, uchar* nghArr) { return d; }

  static inline auto iter(uchar* edge_start, uintE degree, uintE id)
      -> encodings::byte::iter<W> {
    return encodings::byte::iter<W>(edge_start, degree, id);
  }

  static inline size_t intersect(uchar* l1, uchar* l2, uintE l1_size,
                                 uintE l2_size, uintE l1_src, uintE l2_src) {
    return encodings::byte::intersect<W>(l1, l2, l1_size, l2_size, l1_src,
                                         l2_src);
  }

  template <class F>
  static inline size_t intersect_f(uchar* l1, uchar* l2, uintE l1_size,
                                   uintE l2_size, uintE l1_src, uintE l2_src,
                                   const F& f) {
    return encodings::byte::intersect_f<W>(l1, l2, l1_size, l2_size, l1_src,
                                           l2_src, f);
  }

  template <class I>
  static inline long sequentialCompressEdgeSet(uchar* edgeArray,
                                               size_t current_offset,
                                               uintT degree, uintE source,
                                               I& it) {
    return encodings::byte::sequentialCompressEdgeSet(edgeArray, current_offset,
                                                      degree, source, it);
  }

  template <class P, class O>
  static inline void filter(P pred, uchar* edge_start, const uintE& source,
                            const uintE& degree, std::tuple<uintE, W>* tmp,
                            O& out) {
    return encodings::byte::filter(pred, edge_start, source, degree, tmp, out);
  }

  template <class P>
  static inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                            const uintE& degree,
                            std::tuple<uintE, W>* tmp_space, bool par = true) {
    return encodings::byte::pack(pred, edge_start, source, degree, tmp_space);
  }

  template <class E, class M, class Monoid>
  static inline E map_reduce(uchar* edge_start, const uintE& source,
                             const uintT& degree, M& m, Monoid& reduce,
                             const bool par = true) {
    return encodings::byte::map_reduce<W, E, M, Monoid>(edge_start, source, degree, m, reduce,
                                          par);
  }

  template <class T>
  static inline void decode(T t, uchar* edge_start, const uintE& source,
                            const uintT& degree, const bool par = true) {
    return encodings::byte::decode<W, T>(t, edge_start, source, degree, par);
  }

  static inline std::tuple<uintE, W> get_ith_neighbor(uchar* edge_start,
                                                      uintE source,
                                                      uintE degree, size_t i) {
    return encodings::byte::get_ith_neighbor<W>(edge_start, source, degree, i);
  }
};

// TODO(laxmand): fully deprecate byte-pd from the repo. We should not support
// the byte-pd format as there is no way to implement Pack(..) properly (this is
// why we wrote the bytepd_amortized format).
//
// template <class W>
// struct bytepd_decode {
//  template <class I>
//  inline long sequentialCompressEdgeSet(uchar* edgeArray, size_t
//  current_offset,
//                                        uintT degree, uintE source, I& it) {
//    return encodings::bytepd::sequentialCompressEdgeSet(
//        edgeArray, current_offset, degree, source, it);
//  }
//
//  template <class P, class O>
//  inline void filter(P pred, uchar* edge_start, const uintE& source,
//                     const uintE& degree, std::tuple<uintE, W>* tmp, O& out) {
//    return encodings::bytepd::filter(pred, edge_start, source, degree, tmp,
//                                     out);
//  }
//
//  template <class P>
//  inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
//                     const uintE& degree, std::tuple<uintE, W>* tmp_space,
//                     bool par = true) {
//    return encodings::bytepd::pack(pred, edge_start, source, degree,
//    tmp_space,
//                                   par);
//  }
//
//  template <class E, class M, class R>
//  inline E map_reduce(uchar* edge_start, const uintE& source,
//                         const uintT& degree, E id, M& m, R& r,
//                         const bool par = true) {
//    return encodings::bytepd::map_reduce(edge_start, source, degree, id, m, r,
//                                         par);
//  }
//
//  template <class T>
//  inline void decode(T t, uchar* edge_start, const uintE& source,
//                     const uintT& degree, const bool par = true) {
//    return encodings::bytepd::decode(t, edge_start, source, degree, par);
//  }
//
//  static inline std::tuple<uintE, W> get_ith_neighbor(uchar* edge_start,
//                                                 uintE source, uintE degree,
//                                                 size_t i) {
//    assert(
//        false);  // currently not supported in this format; use
//        bytepd_amortized
//    return i;
//  }
//};

template <class W>
struct bytepd_amortized_decode {
  using iter_type = encodings::bytepd_amortized::iter<W>;

  static inline size_t intersect(uchar* l1, uchar* l2, uintE l1_size,
                                 uintE l2_size, uintE l1_src, uintE l2_src) {
    return encodings::bytepd_amortized::intersect<W>(l1, l2, l1_size, l2_size,
                                                     l1_src, l2_src);
  }

  template <class F>
  static inline size_t intersect_f(uchar* l1, uchar* l2, uintE l1_size,
                                   uintE l2_size, uintE l1_src, uintE l2_src,
                                   const F& f) {
    return encodings::bytepd_amortized::intersect_f<W>(l1, l2, l1_size, l2_size,
                                                       l1_src, l2_src, f);
  }

  static inline auto iter(uchar* edge_start, uintE degree, uintE id)
      -> encodings::bytepd_amortized::iter<W> {
    return encodings::bytepd_amortized::iter<W>(edge_start, degree, id);
  }

  template <class I>
  static inline long sequentialCompressEdgeSet(uchar* edgeArray,
                                               size_t current_offset,
                                               uintT degree, uintE source,
                                               I& it) {
    return encodings::bytepd_amortized::sequentialCompressEdgeSet(
        edgeArray, current_offset, degree, source, it);
  }

  template <class P, class O>
  static inline void filter(P pred, uchar* edge_start, const uintE& source,
                            const uintE& degree, std::tuple<uintE, W>* tmp,
                            O& out) {
    return encodings::bytepd_amortized::filter(pred, edge_start, source, degree,
                                               tmp, out);
  }

  template <class P>
  static inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                            const uintE& degree,
                            std::tuple<uintE, W>* tmp_space, bool par = true) {
    return encodings::bytepd_amortized::pack(pred, edge_start, source, degree,
                                             tmp_space, par);
  }

  template <class E, class M, class Monoid>
  static inline E map_reduce(uchar* edge_start, const uintE& source,
                             const uintT& degree, M& m, Monoid& reduce,
                             const bool par = true) {
    return encodings::bytepd_amortized::map_reduce<W, E>(edge_start, source,
                                                         degree, m, reduce, par);
  }

  template <class T>
  static inline void decode(T t, uchar* edge_start, const uintE& source,
                            const uintT& degree, const bool par = true) {
    return encodings::bytepd_amortized::decode<W, T>(t, edge_start, source,
                                                     degree, par);
  }

  static inline size_t get_virtual_degree(uintE d, uchar* nghArr) {
    return encodings::bytepd_amortized::get_virtual_degree(d, nghArr);
  }

  template <class T>
  static inline void decode_block_seq(T t, uchar* edge_start,
                                      const uintE& source, const uintT& degree,
                                      uintE block_size, uintE block_num) {
    return encodings::bytepd_amortized::decode_block_seq<W, T>(
        t, edge_start, source, degree, block_size, block_num);
  }

  static inline std::tuple<uintE, W> get_ith_neighbor(uchar* edge_start,
                                                      uintE source,
                                                      uintE degree, size_t i) {
    return encodings::bytepd_amortized::get_ith_neighbor<W>(edge_start, source,
                                                            degree, i);
  }
};
