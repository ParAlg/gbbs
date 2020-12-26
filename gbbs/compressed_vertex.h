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

template <class W, class C>
struct compressed_neighbors {

  uchar* neighbors;
  uintE degree;
  uintE id;

  template <class F>
  inline void mapNghs(F& f, bool parallel = true) {
    auto T = [&](const uintE& src, const uintE& target, const W& weight,
                 const uintT& edgeNumber) {
      f(src, target, weight);
      return true;
    };
    C::template decode<W>(T, neighbors, id, degree, parallel);
  }

  template <class F>
  inline void mapNghsWithIndex(F& f, bool parallel = true) {
    auto T = [&](const uintE& src, const uintE& target, const W& weight,
                 const uintT& edgeNumber) {
      f(src, target, weight, edgeNumber);
      return true;
    };
    C::template decode<W>(T, neighbors, id, degree, parallel);
  }

  template <class F, class G>
  inline void copyNghs(uintT offset, F& f, G& g, bool parallel=true) {
    auto T = [&](const uintE& src, const uintE& target, const W& weight,
                 const uintT& edgeNumber) {
      auto val = f(src, target, weight);
      g(target, offset + edgeNumber, val);
      return true;
    };
    C::template decode<W>(T, neighbors, id, degree, parallel);
  }

  template <class F>
  inline size_t countNghs(F& f, bool parallel = true) {
    auto monoid = pbbs::addm<size_t>();
    return C::template map_reduce<W>(neighbors, id, degree, f, monoid, parallel);
  }

  template <class M, class Monoid>
  inline typename Monoid::T reduceNghs(M& m, Monoid& r) {
    return C::template map_reduce<W>(neighbors, id, degree, m, r);
  }

  template <class P, class O>
  inline void filterNghs(P& pred, std::tuple<uintE, W>* tmp, O& out) {
    C::template filter<W, P, O>(pred, neighbors, id, degree, tmp, out);
  }

  template <class P>
  inline size_t packNghs(P& pred, std::tuple<uintE, W>* tmp) {
    return C::template pack<W>(pred, neighbors, id, degree, tmp);
  }

  inline size_t intersect(compressed_neighbors<W, C>* other) {
    return C::template intersect<W>(neighbors, other->neighbors,
        degree, other->degree, id, other->id);
  }

  template <class F>
  inline size_t intersect_f(compressed_neighbors<W, C>* other, const F& f) {
    return C::template intersect<W>(neighbors, other->neighbors,
        degree, other->degree, id, other->id, f);
  }

  template <class F>
  inline size_t intersect_f_par(compressed_neighbors<W, C>* other, const F& f) {
    return intersect_f(other, f);
  }

  inline size_t calculateTemporarySpace() {
  #if defined(PD) || defined(AMORTIZEDPD)
    return (degree <= PD_PACK_THRESHOLD) ? 0 : (degree / kTemporarySpaceConstant);
  #else
    return 0;
  #endif
  }

  inline size_t calculateTemporarySpaceBytes() {
    return calculateTemporarySpace() * sizeof(std::tuple<uintE, W>);
  }

  inline std::tuple<uintE, W> get_ith_neighbor(size_t i) {
    return C::template get_ith_neighbor<W>(neighbors, id, degree, i);
  }

  template <class C>
  inline size_t getVirtualDegree() {
    return C::get_virtual_degree(degree, neighbors);
  }

  // ======== Internal primitives used by EdgeMap implementations =======

  template <class F, class G, class VS>
  __attribute__((always_inline)) inline void decodeBreakEarly(const VS& vs,
                                   F& f, const G& g, bool parallel=false) {
    auto T = [&](const uintE& src, const uintE& target, const W& weight,
                 const uintT& edgeNumber) __attribute__((always_inline)) {
      if (vs.isIn(target)) {
        auto m = f.update(target, src, weight);
        g(src, m);
      }
      return f.cond(src);
    };
    C::template decode<W>(T, neighbors, id, degree, parallel);
  }


  template <class F, class G>
  inline void decodeNghs(F& f, G& g, bool parallel = true) {
    auto T = [&](const uintE& src, const uintE& target, const W& weight,
                 const uintT& edgeNumber) {
      if (f.cond(target)) {
        auto m = f.updateAtomic(src, target, weight);
        g(target, m);
      }
      return true;
    };
    C::template decode<W>(T, neighbors, id, degree, parallel);
  }

  template <class F, class G>
  inline size_t decodeNghsSparseSeq(uintT offset, F& f, G& g) {
    size_t k = 0;
    auto T = [&](const uintE& src, const uintE& target, const W& weight,
                 const uintT& edgeNumber) {
      if (f.cond(target)) {
        auto m = f.updateAtomic(src, target, weight);
        if (g(target, offset + k, m)) {
          k++;
        }
      }
      return true;
    };
    C::template decode<W>(T, neighbors, id, degree, /* parallel = */ false);
    return k;
  }

  template <class F, class G>
  inline size_t decodeNghsSparseBlock(uintT offset, uintE block_size, uintE block_num,
                                      F& f, G& g) {
    size_t k = 0;
    auto T = [&](const uintE& src, const uintE& target, const W& weight) {
      if (f.cond(target)) {
        auto m = f.updateAtomic(src, target, weight);
        if (g(target, offset + k, m)) {
          k++;
        }
      }
    };
    C::template decode_block_seq<W>(T, neighbors, id, degree, block_size, block_num);
    return k;
  }

  template <class F, class G>
  inline size_t decode_block(uintT offset, uintE block_num, F& f, G& g) {
    size_t k = 0;
    auto T = [&](const uintE& target, const W& weight, const size_t& edge_id) {
      if (f.cond(target)) {
        auto m = f.updateAtomic(vtx_id, target, weight);
        if (g(target, offset + k, m)) {
          k++;
        }
      }
    };
    C::template decode_block<W>(T, neighbors, id, degree, block_num);
    return k;
  }

};  // struct compressed_neighbors


template <class W, class C>
struct compressed_symmetric_vertex {
  using vertex = compressed_symmetric_vertex<W, C>;
  using edge_type = uchar;

  edge_type* neighbors;
  uintE degree;
  uintE id;

  compressed_symmetric_vertex(edge_type* n, vertex_data& vdata) {
    neighbors = n + vdata.offset;
    degree = vdata.degree;
  }

  compressed_neighbors<W, C> in_neighbors() {
    return compressed_neighbors<W, C>(neighbors, degree, id); }
  compressed_neighbors<W, C> out_neighbors() {
    return in_neighbors(); }

  uintE in_degree() { return degree; }
  uintE out_degree() { return degree; }

  constexpr static uintE getInternalBlockSize() {
    return PARALLEL_DEGREE;
  }

  // TODO: used?
  void clear() {}
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

  compressed_neighbors<W, C> in_neighbors() {
    return compressed_neighbors<W, C>(inNeighbors, inDegree, id); }
  compressed_neighbors<W, C> out_neighbors() {
    return compressed_neighbors<W, C>(outNeighbors, outDegree, id); }

  uintE in_degree() { return inDegree; }
  uintE out_degree() { return outDegree; }

  constexpr static uintE getInternalBlockSize() {
    return PARALLEL_DEGREE;
  }

  // TODO: used?
  void clear() {}
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
