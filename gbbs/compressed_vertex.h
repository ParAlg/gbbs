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

// Functors for representing a compressed vertex. The underlying storage
// is a pointer to a uchar* and a degree (two pointers and two degrees in the
// case of directed/asymmetric graphs).
//
// The classes have two templates, W and C:
// W : the weight type of the underlying graph (pbbslib::empty if unweighted)
// C : the compression format used. See encodings/decoders.h.
//
// To avoid duplication, a lot of the implementation is factored out into
// cvertex; the functions in compressed_symmetric/compressed_asymmetric just
// call out to the appropriate cvertex or compression-struct method.
#pragma once

#include <tuple>
#include <utility>

#include "encodings/decoders.h"
#include "macros.h"
#include "pbbslib/monoid.h"

namespace gbbs {

namespace cvertex {

template <class W, class C, class F, class G, class VS>
__attribute__((always_inline)) inline void decodeNghsBreakEarly(const uintE vtx_id, const uintE d, uchar* nghArr, const VS& vs,
                                 F& f, const G& g, bool parallel=false) {
  auto T = [&](const uintE& src, const uintE& target, const W& weight,
               const uintT& edgeNumber) __attribute__((always_inline)) {
    if (vs.isIn(target)) {
      auto m = f.update(target, src, weight);
      g(src, m);
    }
    return f.cond(src);
  };
  C::template decode<W>(T, nghArr, vtx_id, d, parallel);
}

template <class W, class C, class F, class G>
inline void decodeNghs(uintE vtx_id, uintE d, uchar* nghArr, F& f, G& g,
                       bool parallel = true) {
  auto T = [&](const uintE& src, const uintE& target, const W& weight,
               const uintT& edgeNumber) {
    if (f.cond(target)) {
      auto m = f.updateAtomic(src, target, weight);
      g(target, m);
    }
    return true;
  };
  C::template decode<W>(T, nghArr, vtx_id, d, parallel);
}

template <class W, class C, class F, class G, class H>
inline void decodeNghsSparse(uintE vtx_id, uintE d, uchar* nghArr, uintT o,
                             F& f, G& g, H& h, bool parallel = true) {
  auto T = [&](const uintE& src, const uintE& target, const W& weight,
               const uintT& edgeNumber) {
    if (f.cond(target)) {
      auto m = f.updateAtomic(src, target, weight);
      g(target, o + edgeNumber, m);
    } else {
      h(target, o + edgeNumber);
    }
    return true;
  };
  C::template decode<W>(T, nghArr, vtx_id, d, parallel);
}

template <class W, class C, class F, class G>
inline size_t decodeNghsSparseSeq(uintE vtx_id, uintE d, uchar* nghArr, uintT o,
                                  F& f, G& g) {
  size_t k = 0;
  auto T = [&](const uintE& src, const uintE& target, const W& weight,
               const uintT& edgeNumber) {
    if (f.cond(target)) {
      auto m = f.updateAtomic(src, target, weight);
      if (g(target, o + k, m)) {
        k++;
      }
    }
    return true;
  };
  C::template decode<W>(T, nghArr, vtx_id, d, /* parallel = */ false);
  return k;
}

template <class W, class C, class F, class G>
inline size_t decodeNghsSparseBlock(uintE vtx_id, uintE d, uchar* nghArr,
                                    uintT o, uintE block_size, uintE block_num,
                                    F& f, G& g) {
  size_t k = 0;
  auto T = [&](const uintE& src, const uintE& target, const W& weight) {
    if (f.cond(target)) {
      auto m = f.updateAtomic(src, target, weight);
      if (g(target, o + k, m)) {
        k++;
      }
    }
  };
  C::template decode_block_seq<W>(T, nghArr, vtx_id, d, block_size, block_num);
  return k;
}

template <class W, class C, class F, class G>
inline size_t decode_block(uintE vtx_id, uintE d, uchar* nghArr, uintT o,
                           uintE block_num, F& f, G& g) {
  size_t k = 0;
  auto T = [&](const uintE& target, const W& weight, const size_t& edge_id) {
    if (f.cond(target)) {
      auto m = f.updateAtomic(vtx_id, target, weight);
      if (g(target, o + k, m)) {
        k++;
      }
    }
  };
  C::template decode_block<W>(T, nghArr, vtx_id, d, block_num);
  return k;
}

template <class W, class C, class F>
inline void mapNghs(uintE vtx_id, uintE d, uchar* nghArr, F& f,
                    bool parallel = true) {
  auto T = [&](const uintE& src, const uintE& target, const W& weight,
               const uintT& edgeNumber) {
    f(src, target, weight);
    return true;
  };
  C::template decode<W>(T, nghArr, vtx_id, d, parallel);
}

template <class W, class C, class F>
inline void mapNghsWithIndex(uintE vtx_id, uintE d, uchar* nghArr, F& f,
                             bool parallel = true) {
  auto T = [&](const uintE& src, const uintE& target, const W& weight,
               const uintT& edgeNumber) {
    f(src, target, weight, edgeNumber);
    return true;
  };
  C::template decode<W>(T, nghArr, vtx_id, d, parallel);
}

template <class W, class C, class F, class G>
inline void copyNghs(uintE vtx_id, uintE d, uchar* nghArr, uintT o, F& f,
                     G& g, bool parallel=true) {
  auto T = [&](const uintE& src, const uintE& target, const W& weight,
               const uintT& edgeNumber) {
    auto val = f(src, target, weight);
    g(target, o + edgeNumber, val);
    return true;
  };
  C::template decode<W>(T, nghArr, vtx_id, d, parallel);
}

template <class W, class C, class F>
inline size_t countNghs(uintE vtx_id, uintE d, uchar* nghArr, F& f,
                        bool parallel = true) {
  auto monoid = pbbs::addm<size_t>();
  return C::template map_reduce<W, size_t>(nghArr, vtx_id, d, f, monoid, parallel);
}

template <class W, class C, class E, class M, class Monoid>
inline E reduceNghs(uintE vtx_id, uintE d, uchar* nghArr, M& m, Monoid& r) {
  return C::template map_reduce<W, E>(nghArr, vtx_id, d, m, r);
}

template <class W, class C, class P, class O>
inline void filterNghs(uintE vtx_id, uintE d, uchar* nghArr, P& pred,
                       std::tuple<uintE, W>* tmp, O& out) {
  C::template filter<W, P, O>(pred, nghArr, vtx_id, d, tmp, out);
}

template <class W, class C, class P>
inline size_t packNghs(uintE vtx_id, uintE d, uchar* nghArr, P& pred,
                       std::tuple<uintE, W>* tmp) {
  return C::template pack<W>(pred, nghArr, vtx_id, d, tmp);
}

inline size_t calculateTemporarySpace(uintE deg) {
#if defined(PD) || defined(AMORTIZEDPD)
  return (deg <= PD_PACK_THRESHOLD) ? 0 : (deg / kTemporarySpaceConstant);
#else
  return 0;
#endif
}

template <class C>
inline size_t getVirtualDegree(uintE d, uchar* nghArr) {
  return C::get_virtual_degree(d, nghArr);
}
}  // namespace cvertex

template <class W, class C>
struct compressed_symmetric_vertex {
  using vertex = compressed_symmetric_vertex<W, C>;
  using edge_type = uchar;

  edge_type* neighbors;
  uintE degree;

  compressed_symmetric_vertex(edge_type* n, vertex_data& vdata) {
    neighbors = n + vdata.offset;
    degree = vdata.degree;
  }

  edge_type* getInNeighbors() { return neighbors; }
  edge_type* getOutNeighbors() { return neighbors; }
  uintE getInNeighbor(intT j) {
    assert(false);
    return -1;
  }  // should not be called
  uintE getOutNeighbor(intT j) {
    assert(false);
    return -1;
  }  // should not be called
  uintE getInDegree() { return degree; }
  uintE getOutDegree() { return degree; }
  uintE getInVirtualDegree() {
    return cvertex::getVirtualDegree<C>(degree, getInNeighbors());
  }
  uintE getOutVirtualDegree() {
    return cvertex::getVirtualDegree<C>(degree, getOutNeighbors());
  }

  constexpr static uintE getInternalBlockSize() {
    return PARALLEL_DEGREE;
  }
  uintE getNumInBlocks() { return C::get_num_blocks(neighbors, degree); }
  uintE getNumOutBlocks() { return getNumInBlocks(); }
  inline uintE in_block_degree(uintE block_num) {
    return C::get_block_degree(neighbors, degree, block_num);
  }
  inline uintE out_block_degree(uintE block_num) {
    return in_block_degree(block_num);
  }

  void flipEdges() {}
  void clear() {}

  auto getOutIter(uintE id)  {
    return C::template iter<W>(getOutNeighbors(), getOutDegree(), id);
  }
  auto getInIter(uintE id) {
    return C::template iter<W>(getInNeighbors(), getInDegree(), id);
  }

  template <class VS, class F, class G>
 __attribute__((always_inline)) inline void decodeInNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g,
                                    bool parallel = 0) {
    cvertex::decodeNghsBreakEarly<W, C, F, G, VS>(
        vtx_id, getInDegree(), getInNeighbors(), vertexSubset, f, g, parallel);
  }

  template <class VS, class F, class G>
 __attribute__((always_inline)) inline void decodeOutNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g, bool parallel = 0) {
    cvertex::decodeNghsBreakEarly<W, C, F, G, VS>(
        vtx_id, getOutDegree(), getOutNeighbors(), vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeInNgh(uintE vtx_id, F& f, G& g) {
    cvertex::decodeNghs<W, C, F, G>(vtx_id, getInDegree(), getInNeighbors(), f,
                                    g);
  }

  template <class F, class G>
  inline void decodeOutNgh(uintE vtx_id, F& f, G& g) {
    cvertex::decodeNghs<W, C, F, G>(vtx_id, getOutDegree(), getOutNeighbors(),
                                    f, g);
  }

  template <class F, class G, class H>
  inline void decodeInNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    cvertex::decodeNghsSparse<W, C, F, G, H>(vtx_id, getInDegree(),
                                             getInNeighbors(), o, f, g, h, parallel);
  }

  template <class F, class G, class H>
  inline void decodeOutNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    cvertex::decodeNghsSparse<W, C, F, G, H>(vtx_id, getOutDegree(),
                                             getOutNeighbors(), o, f, g, h, parallel);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseSeq(uintE vtx_id, uintT o, F& f, G& g) {
    return cvertex::decodeNghsSparseSeq<W, C, F, G>(vtx_id, getInDegree(),
                                                    getInNeighbors(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(uintE vtx_id, uintT o, F& f, G& g) {
    return cvertex::decodeNghsSparseSeq<W, C, F, G>(vtx_id, getOutDegree(),
                                                    getOutNeighbors(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                               G& g) {
    return cvertex::decode_block<W, C, F, G>(
        vtx_id, getOutDegree(), getOutNeighbors(), o, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                              G& g) {
    return cvertex::decode_block<W, C, F, G>(
        vtx_id, getInDegree(), getInNeighbors(), o, block_num, f, g);
  }


  template <class F, class G>
  inline size_t decodeOutNghSparseBlock(uintE vtx_id, uintT o, uintE block_size,
                                        uintE block_num, F& f, G& g) {
    return cvertex::decodeNghsSparseBlock<W, C, F, G>(
        vtx_id, getOutDegree(), getOutNeighbors(), o, block_size, block_num, f,
        g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseBlock(uintE vtx_id, uintT o, uintE block_size,
                                       uintE block_num, F& f, G& g) {
    return cvertex::decodeNghsSparseBlock<W, C, F, G>(
        vtx_id, getInDegree(), getInNeighbors(), o, block_size, block_num, f,
        g);
  }

  template <class F, class G>
  inline void copyInNgh(uintE vtx_id, uintT o, F& f, G& g) {
    cvertex::copyNghs<W, C, F, G>(vtx_id, getInDegree(), getInNeighbors(), o, f,
                                  g);
  }

  template <class F, class G>
  inline void copyOutNgh(uintE vtx_id, uintT o, F& f, G& g) {
    cvertex::copyNghs<W, C, F, G>(vtx_id, getOutDegree(), getOutNeighbors(), o,
                                  f, g);
  }

  template <class F>
  inline void mapInNgh(uintE vtx_id, F& f, bool parallel = 1) {  // TODO
    cvertex::mapNghs<W, C, F>(vtx_id, getInDegree(), getInNeighbors(), f,
                              parallel);
  }

  template <class F>
  inline void mapOutNgh(uintE vtx_id, F& f, bool parallel = 1) {  // TODO
    cvertex::mapNghs<W, C, F>(vtx_id, getOutDegree(), getOutNeighbors(), f,
                              parallel);
  }

  template <class F>
  inline void mapOutNghWithIndex(uintE vtx_id, F& f, bool parallel = 1) {  // TODO
    cvertex::mapNghsWithIndex<W, C, F>(
        vtx_id, getOutDegree(), getOutNeighbors(), f, parallel);
  }

  inline std::tuple<uintE, W> get_ith_out_neighbor(uintE vtx_id, size_t i) {
    return C::template get_ith_neighbor<W>(getOutNeighbors(), vtx_id, getOutDegree(), i);
  }

  inline std::tuple<uintE, W> get_ith_in_neighbor(uintE vtx_id, size_t i) {
    return C::template get_ith_neighbor<W>(getInNeighbors(), vtx_id, getInDegree(), i);
  }

  template <class F>
  inline size_t countInNgh(uintE vtx_id, F& f, bool parallel = 1) {
    return cvertex::countNghs<W, C, F>(vtx_id, getInDegree(), getInNeighbors(),
                                       f, parallel);
  }

  template <class F>
  inline size_t countOutNgh(uintE vtx_id, F& f, bool parallel = 1) {
    return cvertex::countNghs<W, C, F>(vtx_id, getOutDegree(),
                                       getOutNeighbors(), f, parallel);
  }

  template <class E, class M, class Monoid>
  inline E reduceInNgh(uintE vtx_id, M& m, Monoid& r) {
    return cvertex::reduceNghs<W, C, E, M, Monoid>(vtx_id, getInDegree(),
                                              getInNeighbors(), m, r);
  }

  template <class E, class M, class Monoid>
  inline E reduceOutNgh(uintE vtx_id, M& m, Monoid& r) {
    return cvertex::reduceNghs<W, C, E, M, Monoid>(vtx_id, getOutDegree(),
                                              getOutNeighbors(), m, r);
  }

  template <class P>
  inline size_t packInNgh(uintE vtx_id, P& pred, std::tuple<uintE, W>* tmp) {
    size_t deg = cvertex::packNghs<W, C, P>(vtx_id, getInDegree(),
                                            getInNeighbors(), pred, tmp);
    return deg;
  }

  template <class P, class O>
  inline void filterOutNgh(uintE vtx_id, P& p, O& out,
                           std::tuple<uintE, W>* tmp) {
    cvertex::filterNghs<W, C, P, O>(vtx_id, getOutDegree(), getOutNeighbors(),
                                    p, tmp, out);
  }

  template <class P, class O>
  inline void filterInNgh(uintE vtx_id, P& p, O& out,
                          std::tuple<uintE, W>* tmp) {
    cvertex::filterNghs<W, C, P, O>(vtx_id, getInDegree(), getInNeighbors(), p,
                                    tmp, out);
  }

  template <class P>
  inline size_t packOutNgh(uintE vtx_id, P& pred, std::tuple<uintE, W>* tmp) {
    uintE orig_degree = getOutDegree();
    if (orig_degree > 0) {
      size_t deg = cvertex::packNghs<W, C, P>(vtx_id, orig_degree,
                                              getOutNeighbors(), pred, tmp);
      return deg;
    }
    return orig_degree;
  }

  inline size_t intersect(compressed_symmetric_vertex<W, C>* other, long our_id,
                          long other_id) {
    return C::template intersect<W>(getOutNeighbors(), other->getOutNeighbors(),
        getOutDegree(), other->getOutDegree(), our_id, other_id);
  }

  template <class F>
  inline size_t intersect_f(compressed_symmetric_vertex<W, C>* other, long our_id,
                            long other_id, const F& f) {
    return C::template intersect_f<W>(getOutNeighbors(),
        other->getOutNeighbors(), getOutDegree(), other->getOutDegree(), our_id,
        other_id, f);
  }

  template <class F>
  inline size_t intersect_f_par(compressed_symmetric_vertex<W, C>* other, long our_id,
                            long other_id, const F& f) {
    return C::template intersect_f<W>(getOutNeighbors(),
        other->getOutNeighbors(), getOutDegree(), other->getOutDegree(), our_id,
        other_id, f);
  }


  inline size_t calculateOutTemporarySpace() {
    return cvertex::calculateTemporarySpace(getOutDegree());
  }

  inline size_t calculateOutTemporarySpaceBytes() {
    return cvertex::calculateTemporarySpace(getOutDegree()) * sizeof(std::tuple<uintE, W>);
  }

  inline size_t calculateInTemporarySpace() {
    return cvertex::calculateTemporarySpace(getInDegree());
  }
};

template <class W, class C>
struct compressed_asymmetric_vertex {
  using vertex = compressed_symmetric_vertex<W, C>;
  using edge_type = uchar;

  edge_type* inNeighbors;
  edge_type* outNeighbors;
  uintE outDegree;
  uintE inDegree;

  compressed_asymmetric_vertex(edge_type* out_neighbors, vertex_data& out_data, edge_type* in_neighbors, vertex_data& in_data) {
    inNeighbors = in_neighbors + in_data.offset;
    outNeighbors = out_neighbors + out_data.offset;

    inDegree = in_data.degree;
    outDegree = out_data.degree;
  }

  edge_type* getInNeighbors() { return inNeighbors; }
  edge_type* getOutNeighbors() { return outNeighbors; }
  uintE getInNeighbor(uintE j) {
    assert(false);
    return -1;
  }  // should not be called
  uintE getOutNeighbor(uintE j) {
    assert(false);
    return -1;
  }  // should not be called
  uintE getInDegree() { return inDegree; }
  uintE getOutDegree() { return outDegree; }
  uintE getInVirtualDegree() {
    return cvertex::getVirtualDegree<C>(inDegree, getInNeighbors());
  }
  uintE getOutVirtualDegree() {
    return cvertex::getVirtualDegree<C>(outDegree, getOutNeighbors());
  }

  constexpr static uintE getInternalBlockSize() {
    return PARALLEL_DEGREE;
  }
  uintE getNumInBlocks() { return C::get_num_blocks(inNeighbors, inDegree); }
  uintE getNumOutBlocks() { return C::get_num_blocks(outNeighbors, outDegree); }
  inline uintE in_block_degree(uintE block_num) {
    return C::get_block_degree(inNeighbors, inDegree, block_num);
  }
  inline uintE out_block_degree(uintE block_num) {
    return C::get_block_degree(outNeighbors, outDegree, block_num);
  }

  void flipEdges() {
    std::swap(inNeighbors, outNeighbors);
    std::swap(inDegree, outDegree);
  }
  void clear() {}

  auto getOutIter(uintE vtx_id)  {
    return C::template iter<W>(getOutNeighbors(), getOutDegree(), vtx_id);
  }
  auto getInIter(uintE vtx_id)  {
    return C::template iter<W>(getInNeighbors(), getInDegree(), vtx_id);
  }

  template <class VS, class F, class G>
 __attribute__((always_inline)) inline void decodeInNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g,
                                    bool parallel = 0) {
    cvertex::decodeNghsBreakEarly<W, C, F, G, VS>(
        vtx_id, getInDegree(), getInNeighbors(), vertexSubset, f, g, parallel);
  }

  template <class VS, class F, class G>
 __attribute__((always_inline)) inline void decodeOutNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g, bool parallel = 0) {
    cvertex::decodeNghsBreakEarly<W, C, F, G, VS>(
        vtx_id, getOutDegree(), getOutNeighbors(), vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeInNgh(uintE vtx_id, F& f, G& g) {
    cvertex::decodeNghs<W, C, F, G>(vtx_id, getInDegree(), getInNeighbors(), f,
                                    g);
  }

  template <class F, class G>
  inline void decodeOutNgh(uintE vtx_id, F& f, G& g) {
    cvertex::decodeNghs<W, C, F, G>(vtx_id, getOutDegree(), getOutNeighbors(),
                                    f, g);
  }

  template <class F, class G, class H>
  inline void decodeInNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    cvertex::decodeNghsSparse<W, C, F, G, H>(vtx_id, getInDegree(),
                                             getInNeighbors(), o, f, g, h, parallel);
  }

  template <class F, class G, class H>
  inline void decodeOutNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    cvertex::decodeNghsSparse<W, C, F, G, H>(vtx_id, getOutDegree(),
                                             getOutNeighbors(), o, f, g, h, parallel);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseSeq(uintE vtx_id, uintT o, F& f, G& g) {
    return cvertex::decodeNghsSparseSeq<W, C, F, G>(vtx_id, getInDegree(),
                                                    getInNeighbors(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(uintE vtx_id, uintT o, F& f, G& g) {
    return cvertex::decodeNghsSparseSeq<W, C, F, G>(vtx_id, getOutDegree(),
                                                    getOutNeighbors(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                               G& g) {
    return cvertex::decode_block<W, C, F, G>(
        vtx_id, getOutDegree(), getOutNeighbors(), o, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                              G& g) {
    return cvertex::decode_block<W, C, F, G>(
        vtx_id, getInDegree(), getInNeighbors(), o, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseBlock(uintE vtx_id, uintT o, uintE block_size,
                                        uintE block_num, F& f, G& g) {
    return cvertex::decodeNghsSparseBlock<W, C, F, G>(
        vtx_id, getOutDegree(), getOutNeighbors(), o, block_size, block_num, f,
        g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseBlock(uintE vtx_id, uintT o, uintE block_size,
                                       uintE block_num, F& f, G& g) {
    return cvertex::decodeNghsSparseBlock<W, C, F, G>(
        vtx_id, getInDegree(), getInNeighbors(), o, block_size, block_num, f,
        g);
  }

  inline std::tuple<uintE, W> get_ith_out_neighbor(uintE vtx_id, size_t i) {
    return C::template get_ith_neighbor<W>(getOutNeighbors(), vtx_id, getOutDegree(), i);
  }

  inline std::tuple<uintE, W> get_ith_in_neighbor(uintE vtx_id, size_t i) {
    return C::template get_ith_neighbor<W>(getInNeighbors(), vtx_id, getInDegree(), i);
  }

  template <class F, class G>
  inline void copyInNgh(uintE vtx_id, uintT o, F& f, G& g) {
    cvertex::copyNghs<W, C, F, G>(vtx_id, getInDegree(), getInNeighbors(), o, f,
                                  g);
  }

  template <class F, class G>
  inline void copyOutNgh(uintE vtx_id, uintT o, F& f, G& g) {
    cvertex::copyNghs<W, C, F, G>(vtx_id, getOutDegree(), getOutNeighbors(), o,
                                  f, g);
  }

  template <class F>
  inline void mapInNgh(uintE vtx_id, F& f, bool parallel = 1) {
    cvertex::mapNghs<W, C, F>(vtx_id, getInDegree(), getInNeighbors(), f,
                              parallel);
  }

  template <class F>
  // TODO: parallel false?
  inline void mapOutNgh(uintE vtx_id, F& f, bool parallel = 1) {
    cvertex::mapNghs<W, C, F>(vtx_id, getOutDegree(), getOutNeighbors(), f,
                              parallel);
  }

  template <class F>
  inline size_t countInNgh(uintE vtx_id, F& f, bool parallel = 1) {
    return cvertex::countNghs<W, C, F>(vtx_id, getInDegree(), getInNeighbors(),
                                       f, parallel);
  }

  template <class F>
  inline size_t countOutNgh(uintE vtx_id, F& f, bool parallel = 1) {
    return cvertex::countNghs<W, C, F>(vtx_id, getOutDegree(),
                                       getOutNeighbors(), f, parallel);
  }

  template <class E, class M, class Monoid>
  inline E reduceInNgh(uintE vtx_id, M& m, Monoid& r) {
    return cvertex::reduceNghs<W, C, E, M, Monoid>(vtx_id, getInDegree(),
                                              getInNeighbors(), m, r);
  }

  template <class E, class M, class Monoid>
  inline E reduceOutNgh(uintE vtx_id, M& m, Monoid& r) {
    return cvertex::reduceNghs<W, C, E, M, Monoid>(vtx_id, getOutDegree(),
                                              getOutNeighbors(), m, r);
  }

  template <class P, class O>
  inline void filterOutNgh(uintE vtx_id, P& p, O& out,
                           std::tuple<uintE, W>* tmp) {
    cvertex::filterNghs<W, C, P, O>(vtx_id, getOutDegree(), getOutNeighbors(),
                                    p, tmp, out);
  }

  template <class P, class O>
  inline void filterInNgh(uintE vtx_id, P& p, O& out,
                          std::tuple<uintE, W>* tmp) {
    cvertex::filterNghs<W, C, P, O>(vtx_id, getInDegree(), getInNeighbors(), p,
                                    tmp, out);
  }

  template <class P>
  inline size_t packInNgh(uintE vtx_id, P& pred, std::tuple<uintE, W>* tmp) {
    size_t deg = cvertex::packNghs<W, C, P>(vtx_id, getInDegree(),
                                            getInNeighbors(), pred, tmp);
    return deg;
  }

  template <class P>
  inline size_t packOutNgh(uintE vtx_id, P& pred, std::tuple<uintE, W>* tmp) {
    size_t deg = cvertex::packNghs<W, C, P>(vtx_id, getOutDegree(),
                                            getOutNeighbors(), pred, tmp);
    return deg;
  }

  inline size_t intersect(compressed_asymmetric_vertex<W, C>* other, long our_id,
                          long other_id) {
    return C::template intersect<W>(getOutNeighbors(), other->getOutNeighbors(),
        getOutDegree(), other->getOutDegree(), our_id, other_id);
  }

  template <class F>
  inline size_t intersect_f(compressed_asymmetric_vertex<W, C>* other,
                            long our_id, long other_id, const F& f) {
    return C::template intersect_f<W>(getOutNeighbors(),
        other->getOutNeighbors(), getOutDegree(), other->getOutDegree(), our_id,
        other_id, f);
  }

  template <class F>
  inline size_t intersect_f_par(compressed_asymmetric_vertex<W, C>* other,
                            long our_id, long other_id, const F& f) {
    return C::template intersect_f<W>(getOutNeighbors(),
        other->getOutNeighbors(), getOutDegree(), other->getOutDegree(), our_id,
        other_id, f);
  }


  inline size_t calculateOutTemporarySpace() {
    return cvertex::calculateTemporarySpace(getOutDegree());
  }

  inline size_t calculateOutTemporarySpaceBytes() {
    return cvertex::calculateTemporarySpace(getOutDegree()) * sizeof(std::tuple<uintE, W>);
  }

  inline size_t calculateInTemporarySpace() {
    return cvertex::calculateTemporarySpace(getInDegree());
  }
};

// This is us manually partially applying the functors. Generates two compressed
// vertex classes per encoding. Classes are prefixed with:
// csv for "compressed_symmetric_vertex"
// cav for "compressed_asymmetric_vertex"
template <class W>
struct csv_bytepd_amortized
    : compressed_symmetric_vertex<W, bytepd_amortized_decode> {
  using inner = compressed_symmetric_vertex<W, bytepd_amortized_decode>;
  using inner::inner;
};

template <class W>
struct cav_bytepd_amortized
    : compressed_asymmetric_vertex<W, bytepd_amortized_decode> {
  using inner = compressed_asymmetric_vertex<W, bytepd_amortized_decode>;
  using inner::inner;
};

template <class W>
struct cav_byte
    : compressed_asymmetric_vertex<W, byte_decode> {
  using inner = compressed_asymmetric_vertex<W, byte_decode>;
  using inner::inner;
};

template <class W>
struct csv_byte
    : compressed_symmetric_vertex<W, byte_decode> {
  using inner = compressed_symmetric_vertex<W, byte_decode>;
  using inner::inner;
};

}  // namespace gbbs
