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

#include "pbbslib/sequence_ops.h"
#include "macros.h"

namespace intersection {

template <template <typename W> class vertex, class W>
inline size_t intersect(vertex<W>* A, vertex<W>* B, uintE a, uintE b) {
  uintT i = 0, j = 0, nA = A->getOutDegree(), nB = B->getOutDegree();
  auto nghA = A->getOutNeighbors();
  auto nghB = B->getOutNeighbors();
  size_t ans = 0;
  while (i < nA && j < nB) {
    if (std::get<0>(nghA[i]) == std::get<0>(nghB[j]))
      i++, j++, ans++;
    else if (std::get<0>(nghA[i]) < std::get<0>(nghB[j]))
      i++;
    else
      j++;
  }
  return ans;
}

template <template <typename W> class vertex, class W, class F>
inline size_t intersect_f(vertex<W>* A, vertex<W>* B, uintE a, uintE b,
                          const F& f) {
  uintT i = 0, j = 0, nA = A->getOutDegree(), nB = B->getOutDegree();
  auto nghA = A->getOutNeighbors();
  auto nghB = B->getOutNeighbors();
  size_t ans = 0;
  while (i < nA && j < nB) {
    if (std::get<0>(nghA[i]) == std::get<0>(nghB[j])) {
      f(a, b, std::get<0>(nghA[i]));
      i++, j++, ans++;
    } else if (std::get<0>(nghA[i]) < std::get<0>(nghB[j])) {
      i++;
    } else {
      j++;
    }
  }
  return ans;
}

constexpr const size_t _bs_merge_base = 32;
constexpr const size_t _seq_merge_thresh = 2048;

template <class SeqA, class SeqB, class F>
size_t seq_merge_full(SeqA& A, SeqB& B, F& f) {
  using T = typename SeqA::value_type;
  size_t nA = A.size(), nB = B.size();
  size_t i = 0, j = 0;
  size_t ct = 0;
  while (i < nA && j < nB) {
    T& a = A[i];
    T& b = B[j];
    if (a == b) {
      f(a);
      i++;
      j++;
      ct++;
    } else if (a < b) {
      i++;
    } else {
      j++;
    }
  }
  return ct;
}

template <class SeqA, class SeqB, class F>
size_t seq_merge(const SeqA& A, const SeqB& B, const F& f) {
  using T = typename SeqA::value_type;
  size_t nA = A.size();
  size_t ct = 0;
  for (size_t i=0; i < nA; i++) {
    const T& a = A[i];
    size_t mB = pbbslib::binary_search(B, a, std::less<T>());
    const T& b = B[mB];
    if (a == b) {
      f(a);
      ct++;
    }
  }
  return ct;
}

template <class SeqA, class SeqB, class F>
size_t merge(const SeqA& A, const SeqB& B, const F& f) {
  using T = typename SeqA::value_type;
  size_t nA = A.size();
  size_t nB = B.size();
  size_t nR = nA + nB;
  if (nR < _seq_merge_thresh) { // handles (small, small) using linear-merge
    return intersection::seq_merge_full(A, B, f);
  } else if (nB < nA) {
    return intersection::merge(B, A, f);
  } else if (nA < _bs_merge_base) {
    return intersection::seq_merge(A, B, f);
  } else {
    size_t mA = nA/2;
    size_t mB = pbbslib::binary_search(B, A[mA], std::less<T>());
    size_t m_left = 0;
    size_t m_right = 0;
    par_do([&] () { m_left = intersection::merge(A.slice(0, mA), B.slice(0, mB), f);},
     [&] () { m_right = intersection::merge(A.slice(mA, nA), B.slice(mB, nB), f);});
    return m_left + m_right;
  }
}

template <template <typename W> class vertex, class W, class F>
inline size_t intersect_f_par(vertex<W>* A, vertex<W>* B, uintE a, uintE b,
                          const F& f) {
  uintT nA = A->getOutDegree(), nB = B->getOutDegree();
  uintE* nghA = (uintE*)(A->getOutNeighbors());
  uintE* nghB = (uintE*)(B->getOutNeighbors());

  // Will not work if W is not pbbslib::empty, should assert.
  auto seqA = pbbslib::make_sequence<uintE>(nghA, nA);
  auto seqB = pbbslib::make_sequence<uintE>(nghB, nB);

  auto merge_f = [&] (uintE ngh) {
    f(a, b, ngh);
  };
  return intersection::merge(seqA, seqB, merge_f);
}

};  // namespace intersection

namespace vertex_ops {

// allocate temporary space for vertices with degree > kAllocThreshold
const constexpr size_t kAllocThreshold = 10000;
static constexpr uintE kBlockSize = 1024;

template <template <typename W> class vertex, class W, class F, class G,
          class VS>
inline void decodeNghsBreakEarly(vertex<W>* v, uintE vtx_id,
                                 std::tuple<uintE, W>* nghs, uintE d,
                                 VS& vertexSubset, F& f, G& g,
                                 bool parallel = 0) {
  if (!parallel || d < 1000) {
    for (size_t j = 0; j < d; j++) {
      auto nw = nghs[j];
      uintE ngh = std::get<0>(nw);
      if (vertexSubset.isIn(ngh)) {
        auto m = f.update(ngh, vtx_id, std::get<1>(nw));
        g(vtx_id, m);
      }
      if (!f.cond(vtx_id)) break;
    }
  } else {
    size_t b_size = 2048;
    size_t n_blocks = d/b_size + 1;
    par_for(0, n_blocks, 1, [&] (size_t b) {
      if (f.cond(vtx_id)) {
       size_t start = b*b_size;
       size_t end = std::min((b+1)*b_size, static_cast<size_t>(d));
       for (size_t j=start; j<end; j++) {
         if (!f.cond(vtx_id)) break;
         auto nw = nghs[j];
         uintE ngh = std::get<0>(nw);
         if (vertexSubset.isIn(ngh)) {
           auto m = f.updateAtomic(ngh, vtx_id, std::get<1>(nw));
           g(vtx_id, m);
         }
       }
      }
    });
  }
}

// Used by edgeMapDenseForward. For each out-neighbor satisfying cond, call
// updateAtomic.
template <template <typename W> class vertex, class W, class F, class G>
inline void decodeNghs(vertex<W>* v, uintE vtx_id, std::tuple<uintE, W>* nghs,
                       uintE d, F& f, G& g) {
  par_for(0, d, pbbslib::kSequentialForThreshold, [&] (size_t j) {
    auto nw = nghs[j];
    uintE ngh = std::get<0>(nw);
    if (f.cond(ngh)) {
      auto m = f.updateAtomic(vtx_id, ngh, std::get<1>(nw));
      g(ngh, m);
    }
  });
}

// Used by edgeMapSparse. For each out-neighbor satisfying cond, call
// updateAtomic.
template <template <typename W> class vertex, class W, class F, class G,
          class H>
inline void decodeNghsSparse(vertex<W>* v, uintE vtx_id,
                             std::tuple<uintE, W>* nghs, uintE d, uintT o, F& f,
                             G& g, H& h, bool parallel=true) {
  par_for(0, d, pbbslib::kSequentialForThreshold, [&] (size_t j) {
    auto nw = nghs[j];
    uintE ngh = std::get<0>(nw);
    if (f.cond(ngh)) {
      auto m = f.updateAtomic(vtx_id, ngh, std::get<1>(nw));
      g(ngh, o + j, m);
    } else {
      h(ngh, o + j);
    }
  }, parallel);
}

// Used by edgeMapSparse_no_filter. Sequentially decode the out-neighbors,
// and compactly write all neighbors satisfying g().
template <template <typename W> class vertex, class W, class F, class G>
inline size_t decodeNghsSparseSeq(vertex<W>* v, uintE vtx_id,
                                  std::tuple<uintE, W>* nghs, uintE d, uintT o,
                                  F& f, G& g) {
  size_t k = 0;
  for (size_t j = 0; j < d; j++) {
    auto nw = nghs[j];
    uintE ngh = std::get<0>(nw);
    if (f.cond(ngh)) {
      auto m = f.updateAtomic(vtx_id, ngh, std::get<1>(nw));
      bool wrote = g(ngh, o + k, m);
      if (wrote) {
        k++;
      }
    }
  }
  return k;
}

// Used by edgeMapBlocked. Sequentially decode nghs between
// [block_num*KBlockSize, block_num*kBlockSize + block_size)
// and compactly write all neighbors satisfying g().
template <template <typename W> class vertex, class W, class F, class G>
inline size_t decodeNghsSparseBlock(vertex<W>* v, uintE vtx_id,
                                    std::tuple<uintE, W>* nghs, uintE d,
                                    uintT o, uintE block_size, uintE block_num,
                                    F& f, G& g) {
  size_t k = 0;
  size_t start = kEMBlockSize * block_num;
  size_t end = start + block_size;
  for (size_t j = start; j < end; j++) {
    auto nw = nghs[j];
    uintE ngh = std::get<0>(nw);
    if (f.cond(ngh)) {
      auto m = f.updateAtomic(vtx_id, ngh, std::get<1>(nw));
      bool wrote = g(ngh, o + k, m);
      if (wrote) {
        k++;
      }
    }
  }
  return k;
}

template <class W, class F, class G>
inline size_t decode_block(uintE vtx_id, std::tuple<uintE, W>* nghs, uintE d,
                           uintT o, uintE block_num, F& f, G& g) {
  size_t k = 0;
  uintE start = vertex_ops::kBlockSize * block_num;
  uintE end = std::min(start + vertex_ops::kBlockSize, d);
  for (uintE j = start; j < end; j++) {
    auto nw = nghs[j];
    uintE ngh = std::get<0>(nw);
    if (f.cond(ngh)) {
      auto m = f.updateAtomic(vtx_id, ngh, std::get<1>(nw));
      bool wrote = g(ngh, o + k, m);
      if (wrote) {
        k++;
      }
    }
  }
  return k;
}

template <template <typename W> class vertex, class W>
inline std::tuple<uintE, W> get_ith_neighbor(std::tuple<uintE, W>* nghs,
                                             size_t i) {
  return nghs[i];
}

template <template <typename W> class vertex, class W, class F>
inline size_t countNghs(vertex<W>* v, long vtx_id, std::tuple<uintE, W>* nghs,
                        uintE d, F& f, bool parallel = true) {
  if (d == 0) return 0;
  auto im_f = [&](size_t i) -> size_t {
    auto nw = nghs[i];
    return f(vtx_id, std::get<0>(nw), std::get<1>(nw));
  };
  auto im = pbbslib::make_sequence<size_t>(d,im_f);
  return pbbslib::reduce_add(im);
}

template <template <typename W> class vertex, class W, class E, class M,
          class Monoid>
inline E reduceNghs(vertex<W>* v, uintE vtx_id, std::tuple<uintE, W>* nghs, uintE d, M& m, Monoid& reduce) {
  if (d == 0) return reduce.identity;
  auto im_f = [&](size_t i) {
    auto nw = nghs[i];
    return m(vtx_id, std::get<0>(nw), std::get<1>(nw));
  };
  auto im = pbbslib::make_sequence<E>(d, im_f);
  return pbbslib::reduce(im, reduce);
}

template <template <typename W> class vertex, class W, class F>
inline void mapNghs(vertex<W>* v, uintE vtx_id, std::tuple<uintE, W>* nghs,
                    uintE d, F& f, bool parallel) {
  par_for(0, d, pbbslib::kSequentialForThreshold, [&] (size_t j) {
    uintE ngh = v->getOutNeighbor(j);
    f(vtx_id, ngh, v->getOutWeight(j));
  }, parallel);
}

// Expects that out has enough space to hold the output of the filter
template <template <typename W> class vertex, class W, class P, class O>
inline void filterNghs(vertex<W>* v, uintE vtx_id, std::tuple<uintE, W>* nghs,
                       uintE d, P& p, O& out, std::tuple<uintE, W>* tmp) {
  if (d > 0) {
    if (d < vertex_ops::kAllocThreshold) {
      size_t k = 0;
      for (size_t i = 0; i < d; i++) {
        auto nw = nghs[i];
        if (p(vtx_id, std::get<0>(nw), std::get<1>(nw))) {
          out(k++, nw);
        }
      }
    } else {
      auto pc = [&](const std::tuple<uintE, W>& nw) {
        return p(vtx_id, std::get<0>(nw), std::get<1>(nw));
      };
      auto in_im = pbbslib::make_sequence(nghs, d);
      size_t k = pbbslib::filter_out(in_im, pbbslib::make_sequence(tmp, d), pc);
      par_for(0, k, pbbslib::kSequentialForThreshold, [&] (size_t i)
                      { out(i, tmp[i]); });
    }
  }
}

// Caller is responsible for setting the degree on v.
template <template <typename W> class vertex, class W, class Pred>
inline size_t packNghs(vertex<W>* v, uintE vtx_id, Pred& p,
                       std::tuple<uintE, W>* nghs, uintE d,
                       std::tuple<uintE, W>* tmp) {
  if (d < vertex_ops::kAllocThreshold) {
    uintE k = 0;
    for (size_t i = 0; i < d; i++) {
      auto nw = nghs[i];
      uintE ngh = std::get<0>(nw);
      W wgh = std::get<1>(nw);
      if (p(vtx_id, ngh, wgh)) {
        nghs[k++] = std::make_tuple(ngh, wgh);
      }
    }
    return k;
  } else {
    // copy to tmp
    par_for(0, d, pbbslib::kSequentialForThreshold, [&] (size_t i) { tmp[i] = nghs[i]; });
    auto pc = [&](const std::tuple<uintE, W>& nw) {
      return p(vtx_id, std::get<0>(nw), std::get<1>(nw));
    };
//    size_t k = pbbslib::filter_out(pbbslib::make_sequence(tmp, d), pbbslib::make_sequence(nghs, d), pc);
    size_t k = pbbslib::filterf(tmp, nghs, d, pc);
    return k;
  }
}

template <template <typename W> class vertex, class W, class F, class G>
inline void copyNghs(vertex<W>* v, uintE vtx_id, std::tuple<uintE, W>* nghs,
                     uintE d, uintT o, F& f, G& g) {
  par_for(0, d, pbbslib::kSequentialForThreshold, [&] (size_t j) {
    auto nw = nghs[j];
    uintE ngh = std::get<0>(nw);
    auto val = f(vtx_id, ngh, std::get<1>(nw));
    g(ngh, o + j, val);
  });
}

inline size_t calculateTemporarySpace(uintE deg) {
  return (deg < vertex_ops::kAllocThreshold) ? 0 : deg;
}

template <class W>
struct iter {
  std::tuple<uintE, W>* edges;
  uintE degree;
  uintE proc;
  std::tuple<uintE, W> last_edge;
  iter(std::tuple<uintE, W>* _e, uintE _d) : edges(_e), degree(_d), proc(0) {
    if (degree > 0) {
      last_edge = *edges;
      edges++;
      proc++;
    }
  }

  inline std::tuple<uintE, W> cur() { return last_edge; }

  inline std::tuple<uintE, W> next() {
    auto edge = *edges;
    edges++;
    proc++;
    last_edge = edge;
    return edge;
  }

  inline bool has_next() { return proc < degree; }
};

template <class W>
inline iter<W> get_iter(std::tuple<uintE, W>* edges, uintE degree) {
  return iter<W>(edges, degree);
}
}  // namespace vertex_ops

template <class W>
struct symmetric_vertex {
  using vertex = symmetric_vertex<W>;
  using edge_type = std::tuple<uintE, W>;

  edge_type* neighbors;
  uintE degree;

  symmetric_vertex(edge_type* n, vertex_data& vdata) {
    neighbors = (n + vdata.offset);
    degree = vdata.degree;
  }

  /* Unlikely to be necessary since neighbors are likeyl flat
   * allocated in a shared array */
  void clear() { pbbslib::free_array(neighbors); }

  edge_type* getInNeighbors() { return neighbors; }
  edge_type* getOutNeighbors() { return neighbors; }
  uintE getInNeighbor(uintE j) { return std::get<0>(neighbors[j]); }
  uintE getOutNeighbor(uintE j) { return std::get<0>(neighbors[j]); }
  W getInWeight(uintE j) { return std::get<1>(neighbors[j]); }
  W getOutWeight(uintE j) { return std::get<1>(neighbors[j]); }

  uintE getInDegree() { return degree; }
  uintE getOutDegree() { return degree; }
  uintE getInVirtualDegree() { return degree; }
  uintE getOutVirtualDegree() { return degree; }

  constexpr static uintE getInternalBlockSize() {
    return vertex_ops::kBlockSize;
  }
  uintE getNumInBlocks() {
    return pbbs::num_blocks(degree, vertex_ops::kBlockSize);
  }
  uintE getNumOutBlocks() { return getNumInBlocks(); }
  inline uintE in_block_degree(uintE block_num) {
    uintE block_start = block_num * vertex_ops::kBlockSize;
    uintE block_end = std::min(block_start + vertex_ops::kBlockSize, degree);
    return block_end - block_start;  // TODO: check
  }
  inline uintE out_block_degree(uintE block_num) {
    return in_block_degree(block_num);
  }

  void flipEdges() {}

  auto getInIter(uintE id) -> vertex_ops::iter<W> {
    return vertex_ops::get_iter(getInNeighbors(), getInDegree());
  }
  auto getOutIter(uintE id) -> vertex_ops::iter<W> {
    return vertex_ops::get_iter(getOutNeighbors(), getOutDegree());
  }

  inline size_t intersect(symmetric_vertex<W>* other, long our_id,
                          long other_id) {
    return intersection::intersect(this, other, our_id, other_id);
  }

  template <class F>
  inline size_t intersect_f(symmetric_vertex<W>* other, long our_id,
                            long other_id, const F& f) {
    return intersection::intersect_f(this, other, our_id, other_id, f);
  }

  template <class F>
  inline size_t intersect_f_par(symmetric_vertex<W>* other, long our_id,
                            long other_id, const F& f) {
    return intersection::intersect_f_par(this, other, our_id, other_id, f);
  }


  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g,
                                    bool parallel = 0) {
    vertex_ops::decodeNghsBreakEarly<symmetric_vertex, W, F, G, VS>(
        this, vtx_id, getInNeighbors(), getInDegree(), vertexSubset, f, g,
        parallel);
  }

  template <class VS, class F, class G>
  inline void decodeOutNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g,
                                     bool parallel = 0) {
    vertex_ops::decodeNghsBreakEarly<symmetric_vertex, W, F, G, VS>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), vertexSubset, f, g,
        parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(uintE vtx_id, F& f, G& g) {
    vertex_ops::decodeNghs<symmetric_vertex, W, F, G>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), f, g);
  }

  template <class F, class G>
  inline void decodeInNgh(uintE vtx_id, F& f, G& g) {
    vertex_ops::decodeNghs<symmetric_vertex, W, F, G>(
        this, vtx_id, getInNeighbors(), getInDegree(), f, g);
  }

  template <class F, class G, class H>
  inline void decodeOutNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    vertex_ops::decodeNghsSparse<symmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), o, f, g, h);
  }

  template <class F, class G, class H>
  inline void decodeInNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    vertex_ops::decodeNghsSparse<symmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), o, f, g, h, parallel);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(uintE vtx_id, uintT o, F& f, G& g) {
    return vertex_ops::decodeNghsSparseSeq<symmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseSeq(uintE vtx_id, uintT o, F& f, G& g) {
    return vertex_ops::decodeNghsSparseSeq<symmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                               G& g) {
    return vertex_ops::decode_block<W, F>(vtx_id, getOutNeighbors(),
                                          getOutDegree(), o, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                              G& g) {
    return decodeOutBlock(vtx_id, o, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseBlock(uintE vtx_id, uintT o, uintE block_size,
                                        uintE block_num, F& f, G& g) {
    return vertex_ops::decodeNghsSparseBlock<symmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), o, block_size,
        block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseBlock(uintE vtx_id, uintT o, uintE block_size,
                                       uintE block_num, F& f, G& g) {
    return vertex_ops::decodeNghsSparseBlock<symmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), o, block_size, block_num,
        f, g);
  }

  template <class F, class G>
  inline void copyOutNgh(uintE vtx_id, uintT o, F& f, G& g) {
    vertex_ops::copyNghs<symmetric_vertex, W>(this, vtx_id, getOutNeighbors(),
                                             getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline void copyInNgh(uintE vtx_id, uintT o, F& f, G& g) {
    vertex_ops::copyNghs<symmetric_vertex, W>(this, vtx_id, getInNeighbors(),
                                             getInDegree(), o, f, g);
  }

  inline edge_type get_ith_out_neighbor(uintE vtx_id, size_t i) {
    return vertex_ops::get_ith_neighbor<symmetric_vertex, W>(getOutNeighbors(),
                                                            i);
  }

  inline edge_type get_ith_in_neighbor(uintE vtx_id, size_t i) {
    return vertex_ops::get_ith_neighbor<symmetric_vertex, W>(getInNeighbors(),
                                                            i);
  }

  template <class F>
  inline size_t countOutNgh(uintE vtx_id, F& f, bool parallel = true) {
    return vertex_ops::countNghs<symmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), f, parallel);
  }

  template <class F>
  inline size_t countInNgh(uintE vtx_id, F& f, bool parallel = true) {
    return vertex_ops::countNghs<symmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), f, parallel);
  }

  template <class E, class M, class Monoid>
  inline E reduceOutNgh(uintE vtx_id, M& m, Monoid& reduce) {
    return vertex_ops::reduceNghs<symmetric_vertex, W, E, M, Monoid>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), m, reduce);
  }

  template <class E, class M, class Monoid>
  inline E reduceInNgh(uintE vtx_id, M& m, Monoid& reduce) {
    return vertex_ops::reduceNghs<symmetric_vertex, W, E, M, Monoid>(
        this, vtx_id, getInNeighbors(), getInDegree(), m, reduce);
  }

  template <class F>
  inline void mapOutNgh(uintE vtx_id, F& f, bool parallel = true) {
    vertex_ops::mapNghs<symmetric_vertex, W, F>(this, vtx_id, getOutNeighbors(),
                                               getOutDegree(), f, parallel);
  }

  template <class F>
  inline void mapInNgh(uintE vtx_id, F& f, bool parallel = true) {
    vertex_ops::mapNghs<symmetric_vertex, W, F>(this, vtx_id, getInNeighbors(),
                                               getInDegree(), f, parallel);
  }

  template <class P, class O>
  inline void filterOutNgh(uintE vtx_id, P& p, O& out,
                           edge_type* tmp) {
    vertex_ops::filterNghs<symmetric_vertex, W, P, O>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), p, out, tmp);
  }

  template <class P, class O>
  inline void filterInNgh(uintE vtx_id, P& p, O& out,
                          edge_type* tmp) {
    vertex_ops::filterNghs<symmetric_vertex, W, P, O>(
        this, vtx_id, getInNeighbors(), getInDegree(), p, out, tmp);
  }

  template <class P>
  inline size_t packOutNgh(uintE vtx_id, P& p, edge_type* tmp) {
    uintE k = vertex_ops::packNghs<symmetric_vertex, W, P>(
        this, vtx_id, p, getOutNeighbors(), getOutDegree(), tmp);
    return k;
  }

  template <class P>
  inline size_t packInNgh(uintE vtx_id, P& p, edge_type* tmp) {
    uintE k = vertex_ops::packNghs<symmetric_vertex, W, P>(
        this, vtx_id, p, getInNeighbors(), getInDegree(), tmp);
    setInDegree(k);
    return k;
  }

  inline size_t calculateOutTemporarySpace() {
    return vertex_ops::calculateTemporarySpace(getOutDegree());
  }

  inline size_t calculateInTemporarySpace() {
    return vertex_ops::calculateTemporarySpace(getInDegree());
  }
};

template <class W>
struct asymmetric_vertex {
  using vertex = asymmetric_vertex<W>;
  using edge_type = std::tuple<uintE, W>;

  edge_type* inNeighbors;
  edge_type* outNeighbors;

  uintE inDegree;
  uintE outDegree;

  asymmetric_vertex(edge_type* out_neighbors, vertex_data& out_data,
                    edge_type* in_neighbors, vertex_data& in_data) {
    inNeighbors = in_neighbors + in_data.offset;
    outNeighbors = out_neighbors + out_data.offset;

    inDegree = in_data.degree;
    outDegree = out_data.degree;
  }

  /* Unlikely to be necessary since neighbors are likeyl flat
   * allocated in a shared array */
  void clear() {
    pbbslib::free_array(inNeighbors);
    pbbslib::free_array(outNeighbors);
  }

  edge_type* getInNeighbors() { return inNeighbors; }
  edge_type* getOutNeighbors() { return outNeighbors; }
  uintE getInNeighbor(uintE j) { return std::get<0>(inNeighbors[j]); }
  uintE getOutNeighbor(uintE j) { return std::get<0>(outNeighbors[j]); }
  W getInWeight(uintE j) { return std::get<1>(inNeighbors[j]); }
  W getOutWeight(uintE j) { return std::get<1>(outNeighbors[j]); }

  uintE getInDegree() { return inDegree; }
  uintE getOutDegree() { return outDegree; }
  uintE getInVirtualDegree() { return inDegree; }
  uintE getOutVirtualDegree() { return outDegree; }

  constexpr static uintE getInternalBlockSize() {
    return vertex_ops::kBlockSize;
  }
  uintE getNumInBlocks() {
    return pbbs::num_blocks(inDegree, vertex_ops::kBlockSize);
  }
  uintE getNumOutBlocks() {
    return pbbs::num_blocks(outDegree, vertex_ops::kBlockSize);
  }
  inline uintE in_block_degree(uintE block_num) {
    uintE block_start = block_num * vertex_ops::kBlockSize;
    uintE block_end = std::min(block_start + vertex_ops::kBlockSize, inDegree);
    return block_end - block_start;
  }
  inline uintE out_block_degree(uintE block_num) {
    uintE block_start = block_num * vertex_ops::kBlockSize;
    uintE block_end = std::min(block_start + vertex_ops::kBlockSize, outDegree);
    return block_end - block_start;
  }

  void flipEdges() {
    std::swap(inNeighbors, outNeighbors);
    std::swap(inDegree, outDegree);
  }

  auto getInIter(uintE id) -> vertex_ops::iter<W> {
    return vertex_ops::get_iter(getInNeighbors(), getInDegree());
  }
  auto getOutIter(uintE id) -> vertex_ops::iter<W> {
    return vertex_ops::get_iter(getOutNeighbors(), getOutDegree());
  }

  inline size_t intersect(asymmetric_vertex<W>* other, long our_id,
                          long other_id) {
    return intersection::intersect(this, other, our_id, other_id);
  }

  template <class F>
  inline size_t intersect_f(asymmetric_vertex<W>* other, long our_id,
                            long other_id, const F& f) {
    return intersection::intersect_f(this, other, our_id, other_id, f);
  }

  template <class F>
  inline size_t intersect_f_par(asymmetric_vertex<W>* other, long our_id,
                            long other_id, const F& f) {
    return intersection::intersect_f_par(this, other, our_id, other_id, f);
  }


  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g,
                                    bool parallel = 0) {
    vertex_ops::decodeNghsBreakEarly<asymmetric_vertex, W, F, G, VS>(
        this, vtx_id, getInNeighbors(), getInDegree(), vertexSubset, f, g,
        parallel);
  }

  template <class VS, class F, class G>
  inline void decodeOutNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g,
                                     bool parallel = 0) {
    vertex_ops::decodeNghsBreakEarly<asymmetric_vertex, W, F, G, VS>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), vertexSubset, f, g,
        parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(uintE vtx_id, F& f, G& g) {
    vertex_ops::decodeNghs<asymmetric_vertex, W, F, G>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), f, g);
  }

  template <class F, class G>
  inline void decodeInNgh(uintE vtx_id, F& f, G& g) {
    vertex_ops::decodeNghs<asymmetric_vertex, W, F, G>(
        this, vtx_id, getInNeighbors(), getInDegree(), f, g);
  }

  template <class F, class G, class H>
  inline void decodeOutNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    vertex_ops::decodeNghsSparse<asymmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), o, f, g, h, parallel);
  }

  template <class F, class G, class H>
  inline void decodeInNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    vertex_ops::decodeNghsSparse<asymmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), o, f, g, h, parallel);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(uintE vtx_id, uintT o, F& f, G& g) {
    return vertex_ops::decodeNghsSparseSeq<asymmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseSeq(uintE vtx_id, uintT o, F& f, G& g) {
    return vertex_ops::decodeNghsSparseSeq<asymmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                               G& g) {
    return vertex_ops::decode_block<W, F>(vtx_id, getOutNeighbors(),
                                          getOutDegree(), o, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                              G& g) {
    return vertex_ops::decode_block<W, F>(vtx_id, getInNeighbors(),
                                          getInDegree(), o, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseBlock(uintE vtx_id, uintT o, uintE block_size,
                                        uintE block_num, F& f, G& g) {
    return vertex_ops::decodeNghsSparseBlock<asymmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), o, block_size,
        block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseBlock(uintE vtx_id, uintT o, uintE block_size,
                                       uintE block_num, F& f, G& g) {
    return vertex_ops::decodeNghsSparseBlock<asymmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), o, block_size, block_num,
        f, g);
  }

  template <class F, class G>
  inline void copyOutNgh(uintE vtx_id, uintT o, F& f, G& g) {
    vertex_ops::copyNghs<asymmetric_vertex, W>(this, vtx_id, getOutNeighbors(),
                                              getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline void copyInNgh(uintE vtx_id, uintT o, F& f, G& g) {
    vertex_ops::copyNghs<asymmetric_vertex, W>(this, vtx_id, getInNeighbors(),
                                              getInDegree(), o, f, g);
  }

  inline edge_type get_ith_out_neighbor(uintE vtx_id, size_t i) {
    return vertex_ops::get_ith_neighbor<asymmetric_vertex, W>(getOutNeighbors(),
                                                             i);
  }

  inline edge_type get_ith_in_neighbor(uintE vtx_id, size_t i) {
    return vertex_ops::get_ith_neighbor<asymmetric_vertex, W>(getInNeighbors(),
                                                             i);
  }

  template <class F>
  inline size_t countOutNgh(uintE vtx_id, F& f, bool parallel = true) {
    return vertex_ops::countNghs<asymmetric_vertex, W, F>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), f, parallel);
  }

  template <class F>
  inline size_t countInNgh(uintE vtx_id, F& f, bool parallel = true) {
    return vertex_ops::countNghs<asymmetric_vertex, W, F>(
        this, vtx_id, getInNeighbors(), getInDegree(), f, parallel);
  }

  template <class E, class M, class Monoid>
  inline E reduceOutNgh(uintE vtx_id, M& m, Monoid& reduce) {
    return vertex_ops::reduceNghs<asymmetric_vertex, W, E, M, Monoid>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), m, reduce);
  }

  template <class E, class M, class Monoid>
  inline E reduceInNgh(uintE vtx_id, M& m, Monoid& reduce) {
    return vertex_ops::reduceNghs<asymmetric_vertex, W, E, M, Monoid>(
        this, vtx_id, getInNeighbors(), getInDegree(), m, reduce);
  }

  template <class F>
  inline void mapOutNgh(uintE vtx_id, F& f, bool parallel = true) {
    vertex_ops::mapNghs<asymmetric_vertex, W, F>(this, vtx_id, getOutNeighbors(),
                                                getOutDegree(), f, parallel);
  }

  template <class F>
  inline void mapInNgh(uintE vtx_id, F& f, bool parallel = true) {
    vertex_ops::mapNghs<asymmetric_vertex, W, F>(this, vtx_id, getInNeighbors(),
                                                getInDegree(), f, parallel);
  }

  template <class P, class O>
  inline void filterOutNgh(uintE vtx_id, P& p, O& out,
                           edge_type* tmp) {
    vertex_ops::filterNghs<asymmetric_vertex, W, P, O>(
        this, vtx_id, getOutNeighbors(), getOutDegree(), p, out, tmp);
  }

  template <class P, class O>
  inline void filterInNgh(uintE vtx_id, P& p, O& out,
                          edge_type* tmp) {
    vertex_ops::filterNghs<asymmetric_vertex, W, P, O>(
        this, vtx_id, getInNeighbors(), getInDegree(), p, out, tmp);
  }

  /* Should not use if running an immutable graph algorithm */
  template <class P>
  inline size_t packOutNgh(uintE vtx_id, P& p, edge_type* tmp) {
    uintE k = vertex_ops::packNghs<asymmetric_vertex, W, P>(
        this, vtx_id, p, getOutNeighbors(), getOutDegree(), tmp);
    return k;
  }

  /* Should not use if running an immutable graph algorithm */
  template <class P>
  inline size_t packInNgh(uintE vtx_id, P& p, edge_type* tmp) {
    uintE k = vertex_ops::packNghs<asymmetric_vertex, W, P>(
        this, vtx_id, p, getInNeighbors(), getInDegree(), tmp);
    return k;
  }

  inline size_t calculateOutTemporarySpace() {
    return vertex_ops::calculateTemporarySpace(getOutDegree());
  }

  inline size_t calculateInTemporarySpace() {
    return vertex_ops::calculateTemporarySpace(getInDegree());
  }
};
