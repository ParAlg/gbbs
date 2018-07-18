#pragma once

#include "byte.h"
#include "byte_pd.h"
#include "byte_pd_amortized.h"

struct byte_decode {

  template <class W>
  static inline size_t intersect(uchar* l1, uchar* l2, uintE l1_size, uintE l2_size, uintE l1_src, uintE l2_src) {
    return byte::intersect<W>(l1, l2, l1_size, l2_size, l1_src, l2_src);
  }

  template <class W, class I>
  static inline long sequentialCompressEdgeSet(uchar* edgeArray, size_t current_offset,
                                        uintT degree, uintE source, I& it) {
    return byte::sequentialCompressEdgeSet(edgeArray, current_offset, degree,
                                           source, it);
  }

  template <class W, class P, class O>
  static inline void filter(P pred, uchar* edge_start, const uintE& source,
                     const uintE& degree, tuple<uintE, W>* tmp, O& out) {
    return byte::filter(pred, edge_start, source, degree, tmp, out);
  }

  template <class W, class P>
  static inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                     const uintE& degree, tuple<uintE, W>* tmp_space,
                     bool par = true) {
    return byte::pack(pred, edge_start, source, degree, tmp_space, par);
  }

  template <class W, class E, class M, class R>
  static inline auto map_reduce(uchar* edge_start, const uintE& source,
                         const uintT& degree, E id, M& m, R& r,
                         const bool par = true) {
    return byte::map_reduce(edge_start, source, degree, id, m, r, par);
  }

  template <class W, class T>
  static inline void decode(T t, uchar* edge_start, const uintE& source,
                     const uintT& degree, const bool par = true) {
    return byte::decode<W, T>(t, edge_start, source, degree, par);
  }

  template <class W>
  static inline tuple<uintE, W> get_ith_neighbor(uchar* edge_start, uintE source,
      uintE degree, size_t i) {
    assert(false); // currently not supported in this format; use bytepd_amortized
    return i;
  }
};

struct bytepd_decode {
  template <class W, class I>
  inline long sequentialCompressEdgeSet(uchar* edgeArray, size_t current_offset,
                                        uintT degree, uintE source, I& it) {
    return bytepd::sequentialCompressEdgeSet(edgeArray, current_offset, degree,
                                             source, it);
  }

  template <class W, class P, class O>
  inline void filter(P pred, uchar* edge_start, const uintE& source,
                     const uintE& degree, tuple<uintE, W>* tmp, O& out) {
    return bytepd::filter(pred, edge_start, source, degree, tmp, out);
  }

  template <class W, class P>
  inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                     const uintE& degree, tuple<uintE, W>* tmp_space,
                     bool par = true) {
    return bytepd::pack(pred, edge_start, source, degree, tmp_space, par);
  }

  template <class W, class E, class M, class R>
  inline auto map_reduce(uchar* edge_start, const uintE& source,
                         const uintT& degree, E id, M& m, R& r,
                         const bool par = true) {
    return bytepd::map_reduce(edge_start, source, degree, id, m, r, par);
  }

  template <class T>
  inline void decode(T t, uchar* edge_start, const uintE& source,
                     const uintT& degree, const bool par = true) {
    return bytepd::decode(t, edge_start, source, degree, par);
  }

  template <class W>
  static inline tuple<uintE, W> get_ith_neighbor(uchar* edge_start, uintE source,
      uintE degree, size_t i) {
    assert(false); // currently not supported in this format; use bytepd_amortized
    return i;
  }
};

struct bytepd_amortized_decode {

  template <class W>
  static inline size_t intersect(uchar* l1, uchar* l2, uintE l1_size, uintE l2_size, uintE l1_src, uintE l2_src) {
    return bytepd_amortized::intersect<W>(l1, l2, l1_size, l2_size, l1_src, l2_src);
  }

  template <class W>
  static inline auto iter(uchar* edge_start, uintE degree, uintE id) {
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
                            const uintE& degree, tuple<uintE, W>* tmp, O& out) {
    return bytepd_amortized::filter(pred, edge_start, source, degree, tmp, out);
  }

  template <class W, class P>
  static inline size_t pack(P& pred, uchar* edge_start, const uintE& source,
                            const uintE& degree, tuple<uintE, W>* tmp_space,
                            bool par = true) {
    return bytepd_amortized::pack(pred, edge_start, source, degree, tmp_space,
                                  par);
  }

  template <class W, class E, class M, class R>
  static inline auto map_reduce(uchar* edge_start, const uintE& source,
                                const uintT& degree, E id, M& m, R& r,
                                const bool par = true) {
    return bytepd_amortized::map_reduce<W, E>(edge_start, source, degree, id, m, r,
                                        par);
  }

  template <class W, class T>
  static inline void decode(T t, uchar* edge_start, const uintE& source,
                            const uintT& degree, const bool par = true) {
    return bytepd_amortized::decode<W, T>(t, edge_start, source, degree, par);
  }

  static inline size_t get_virtual_degree(uintE d, uchar* nghArr) {
    return bytepd_amortized::get_virtual_degree(d, nghArr);
  }

  template <class W, class T>
  static inline void decode_block_seq(T t, uchar* edge_start,
                                      const uintE& source, const uintT& degree,
                                      uintE block_size, uintE block_num) {
    return bytepd_amortized::decode_block_seq<W, T>(
        t, edge_start, source, degree, block_size, block_num);
  }

  template <class W>
  static inline tuple<uintE, W> get_ith_neighbor(uchar* edge_start, uintE source,
      uintE degree, size_t i) {
    return bytepd_amortized::get_ith_neighbor<W>(edge_start, source, degree, i);
  }

};
