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

#include "../lib/sequence_ops.h"
#include "../lib/index_map.h"
#include "macros.h"

using namespace std;


namespace intersection {

  // TODO: parallelize
  template <template <typename W> class vertex, class W>
  inline long intersect(vertex<W>* A, vertex<W>* B, uintE a, uintE b) {
    uintT i=0,j=0,nA = A->getOutDegree(), nB = B->getOutDegree();
    auto nghA = A->getOutNeighbors(); auto nghB = B->getOutNeighbors();
    long ans=0;
    while (i < nA && j < nB) {
      if (get<0>(nghA[i]) == get<0>(nghB[j])) i++, j++, ans++;
      else if (get<0>(nghA[i]) < get<0>(nghB[j])) i++;
      else j++;
    }
    return ans;
  }

};

namespace vertex {

  // allocate temporary space for vertices with degree > alloc_threshold
  const constexpr size_t alloc_threshold = 10000;

  template <template <typename W> class vertex, class W, class F,
           class G, class VS>
  inline void decodeNghsBreakEarly(vertex<W>* v, uintE vtx_id,
      tuple<uintE, W>* nghs, uintE d, VS& vertexSubset, F &f, G &g,
      bool parallel=0) {

    if (!parallel || d < 1000) {
      for (size_t j=0; j<d; j++) {
        auto nw = nghs[j];
        uintE ngh = get<0>(nw);
        if (vertexSubset.isIn(ngh)) {
          auto m = f.update(ngh, vtx_id, get<1>(nw));
          g(vtx_id, m);
        }
        if(!f.cond(vtx_id)) break;
      }
    } else {
      parallel_for(size_t j=0; j<d; j++) {
        auto nw = nghs[j];
        uintE ngh = get<0>(nw);
        if (vertexSubset.isIn(ngh)) {
          auto m = f.updateAtomic(ngh, vtx_id, get<1>(nw));
          g(vtx_id, m);
        }
      }
    }
  }

  // Used by edgeMapDenseForward. For each out-neighbor satisfying cond, call
  // updateAtomic.
  template <template <typename W> class vertex, class W, class F, class G>
  inline void decodeNghs(vertex<W>* v, uintE vtx_id, tuple<uintE, W>* nghs,
      uintE d, F &f, G &g) {
    granular_for(j, 0, d, (d > 1000), {
      auto nw = nghs[j];
      uintE ngh = get<0>(nw);
      if (f.cond(ngh)) {
        auto m = f.updateAtomic(vtx_id,ngh,get<1>(nw));
        g(ngh, m);
      }
    });
  }

  // Used by edgeMapSparse. For each out-neighbor satisfying cond, call
  // updateAtomic.
  template <template <typename W> class vertex, class W, class F, class G>
  inline void decodeNghsSparse(vertex<W>* v, uintE vtx_id,
      tuple<uintE, W>* nghs, uintE d, uintT o, F &f, G &g) {
    granular_for(j, 0, d, (d > 1000), {
      auto nw = nghs[j];
      uintE ngh = get<0>(nw);
      if (f.cond(ngh)) {
        auto m = f.updateAtomic(vtx_id, ngh, get<1>(nw));
        g(ngh, o+j, m);
      } else {
        g(ngh, o+j);
      }
    });
  }

  // Used by edgeMapSparse_no_filter. Sequentially decode the out-neighbors,
  // and compactly write all neighbors satisfying g().
  template <template <typename W> class vertex, class W, class F, class G>
  inline size_t decodeNghsSparseSeq(vertex<W>* v, uintE vtx_id,
      tuple<uintE, W>* nghs, uintE d, uintT o, F &f, G &g) {
    size_t k = 0;
    for (size_t j=0; j<d; j++) {
      auto nw = nghs[j];
      uintE ngh = get<0>(nw);
      if (f.cond(ngh)) {
        auto m = f.updateAtomic(vtx_id, ngh, get<1>(nw));
        bool wrote = g(ngh, o+k, m);
        if (wrote) { k++; }
      }
    }
    return k;
  }

  // Used by edgeMapBlocked. Sequentially decode nghs between
  // [block_num*KBlockSize, block_num*kBlockSize + block_size)
  // and compactly write all neighbors satisfying g().
  template <template <typename W> class vertex, class W, class F, class G>
  inline size_t decodeNghsSparseBlock(vertex<W>* v, uintE vtx_id,
      tuple<uintE, W>* nghs, uintE d, uintT o, uintE block_size,
      uintE block_num, F &f, G &g) {
    size_t k = 0;
    size_t start = kEMBlockSize*block_num; size_t end = start + block_size;
    for (size_t j=start; j<end; j++) {
      auto nw = nghs[j];
      uintE ngh = get<0>(nw);
      if (f.cond(ngh)) {
        auto m = f.updateAtomic(vtx_id, ngh, get<1>(nw));
        bool wrote = g(ngh, o+k, m);
        if (wrote) { k++; }
      }
    }
    return k;
  }

  template <template <typename W> class vertex, class W>
  inline tuple<uintE, W> get_ith_neighbor(tuple<uintE, W>* nghs, size_t i) {
    return nghs[i];
  }

  template <template <typename W> class vertex, class W, class F>
  inline size_t countNghs(vertex<W>* v, long vtx_id, tuple<uintE, W>* nghs,
      uintE d, F& f, bool parallel=true) {
    if (d == 0) return 0;
    auto im = make_in_imap<size_t>(d, [&] (size_t i) {
      auto nw = nghs[i];
      return f(vtx_id, get<0>(nw), get<1>(nw));
    });
    return pbbs::reduce_add(im);
  }

  template <template <typename W> class vertex, class W, class E, class M, class R>
  inline E reduceNghs(vertex<W>* v, uintE vtx_id, tuple<uintE, W>* nghs,
      uintE d, E id, M& m, R& r) {
    if (d == 0) return id;
    auto im = make_in_imap<E>(d, [&] (size_t i) {
      auto nw = nghs[i];
      return m(vtx_id, get<0>(nw), get<1>(nw)); });
    return pbbs::reduce(im, r);
  }


  template <template <typename W> class vertex, class W, class F>
  inline void mapNghs(vertex<W>* v, uintE vtx_id, tuple<uintE, W>* nghs,
       uintE d, F& f, bool parallel) {
    granular_for(j, 0, d, (d > 1000 && parallel), {
      uintE ngh = v->getOutNeighbor(j);
      f(vtx_id, ngh, v->getOutWeight(j));
    });
  }

  // Expects that out has enough space to hold the output of the filter
  template <template <typename W> class vertex, class W, class P, class O>
  inline void filterNghs(vertex<W>* v, uintE vtx_id, tuple<uintE, W>* nghs,
       uintE d, P& p, O& out, tuple<uintE, W>* tmp) {
    if (d > 0) {
      if (d < vertex::alloc_threshold) {
        size_t k = 0;
        for (size_t i=0; i<d; i++) {
          auto nw = nghs[i];
          if (p(vtx_id, get<0>(nw), get<1>(nw))) {
            out(k++, nw);
          }
        }
      } else {
        auto pc = [&] (const tuple<uintE, W>& nw) {
          return p(vtx_id, get<0>(nw), get<1>(nw));
        };
        auto in_im = make_array_imap(nghs, d);
        auto s = pbbs::filter(in_im, pc, pbbs::no_flag, tmp);
        size_t k = s.size();
        granular_for(i, 0, k, (k > 2000), {
          out(i, tmp[i]);
        });
      }
    }
  }

  // Caller is responsible for setting the degree on v.
  template <template <typename W> class vertex, class W, class Pred>
  inline size_t packNghs(vertex<W>* v, uintE vtx_id, Pred& p,
                         tuple<uintE, W>* nghs, uintE d, tuple<uintE, W>* tmp) {
    if (d < vertex::alloc_threshold) {
      uintE k = 0;
      for (size_t i=0; i<d; i++) {
        auto nw = nghs[i];
        uintE ngh = get<0>(nw);
        W wgh = get<1>(nw);
        if (p(vtx_id, ngh, wgh)) {
          nghs[k++] = make_tuple(ngh, wgh);
        }
      }
      return k;
    } else {
      // copy to tmp
      parallel_for(size_t i=0; i<d; i++) {
        tmp[i] = nghs[i];
      }
      auto pc = [&] (const tuple<uintE, W>& nw) {
        return p(vtx_id, get<0>(nw), get<1>(nw));
      };
      size_t k = pbbs::filterf(tmp, nghs, d, pc);
      return k;
    }
  }

  template <template <typename W> class vertex, class W, class F, class G>
  inline void copyNghs(vertex<W>* v, uintE vtx_id, tuple<uintE, W>* nghs,
      uintE d, uintT o, F& f, G& g) {
    granular_for(j, 0, d, (d > 1000), {
      auto nw = nghs[j];
      uintE ngh = get<0>(nw);
      auto val = f(vtx_id, ngh, get<1>(nw));
      g(ngh, o+j, val);
    });
  }

  inline size_t calculateTemporarySpace(uintE deg) {
    return (deg < vertex::alloc_threshold) ? 0 : deg;
  }



  template <class W>
  struct iter {
    tuple<uintE, W>* edges;
    uintE degree;
    uintE proc;
    tuple<uintE, W> last_edge;
    iter(tuple<uintE, W>* _e, uintE _d) : edges(_e), degree(_d), proc(0) {
      if (degree > 0) {
        last_edge = *edges;
        edges++;
        proc++;
      }
    }

    inline tuple<uintE, W> cur() {
      return last_edge;
    }

    inline tuple<uintE, W> next() {
      auto edge = *edges;
      edges++;
      proc++;
      last_edge = edge;
      return edge;
    }

    inline bool has_next() {
      return proc < degree;
    }
  };

  template <class W>
  auto get_iter(tuple<uintE, W>* edges, uintE degree) {
    return iter<W>(edges, degree);
  }
}

template <class W>
struct symmetricVertex {
  using wghVtx = symmetricVertex<W>;
  tuple<uintE, W>* neighbors;
  uintE degree;
  symmetricVertex(tuple<uintE, W>* n, uintE d) : neighbors(n), degree(d) { }
  void del() { free(neighbors); }

  tuple<uintE, W>* getInNeighbors () { return neighbors; }
  tuple<uintE, W>* getOutNeighbors () { return neighbors; }
  uintE getInNeighbor(uintE j) { return get<0>(neighbors[j]); }
  uintE getOutNeighbor(uintE j) { return get<0>(neighbors[j]); }
  W getInWeight(uintE j) { return get<1>(neighbors[j]); }
  W getOutWeight(uintE j) { return get<1>(neighbors[j]); }
  void setInNeighbor(uintE j, uintE ngh) { get<0>(neighbors[j]) = ngh; }
  void setOutNeighbor(uintE j, uintE ngh) { get<0>(neighbors[j]) = ngh; }
  void setInWeight(uintE j, W wgh) { get<1>(neighbors[j]) = wgh; }
  void setOutWeight(uintE j, W wgh) { get<1>(neighbors[j]) = wgh; }
  void setInNeighbors(tuple<uintE, W>* _i) { neighbors = _i; }
  void setOutNeighbors(tuple<uintE, W>* _i) { neighbors = _i; }

  uintE getInDegree() { return degree; }
  uintE getOutDegree() { return degree; }
  uintE getInVirtualDegree() { return degree; }
  uintE getOutVirtualDegree() { return degree; }
  void setInDegree(uintE _d) { degree = _d; }
  void setOutDegree(uintE _d) { degree = _d; }
  void flipEdges() { }

  auto getInIter(uintE id) { return vertex::get_iter(getInNeighbors(), getInDegree()); }
  auto getOutIter(uintE id) { return vertex::get_iter(getOutNeighbors(), getOutDegree()); }

  inline size_t intersect(symmetricVertex<W>* other, long our_id, long other_id) {
    return intersection::intersect(this, other, our_id, other_id);
  }

  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(uintE vtx_id, VS& vertexSubset, F &f, G &g,
      bool parallel = 0) {
    vertex::decodeNghsBreakEarly<symmetricVertex, W, F, G, VS>(this, vtx_id,
        getInNeighbors(), getInDegree(), vertexSubset, f, g, parallel);
  }

  template <class VS, class F, class G>
  inline void decodeOutNghBreakEarly(uintE vtx_id, VS& vertexSubset, F &f, G &g,
      bool parallel = 0) {
    vertex::decodeNghsBreakEarly<symmetricVertex, W, F, G, VS>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(uintE vtx_id, F &f, G& g) {
     vertex::decodeNghs<symmetricVertex, W, F, G>(this, vtx_id,
         getOutNeighbors(), getOutDegree(), f, g);
  }

  template <class F, class G>
  inline void decodeInNgh(uintE vtx_id, F &f, G& g) {
     vertex::decodeNghs<symmetricVertex, W, F, G>(this, vtx_id,
         getInNeighbors(), getInDegree(), f, g);
  }

  template <class F, class G>
  inline void decodeOutNghSparse(uintE vtx_id, uintT o, F &f, G &g) {
    vertex::decodeNghsSparse<symmetricVertex, W, F>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline void decodeInNghSparse(uintE vtx_id, uintT o, F &f, G &g) {
    vertex::decodeNghsSparse<symmetricVertex, W, F>(this, vtx_id,
        getInNeighbors(), getInDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(uintE vtx_id, uintT o, F &f, G &g) {
    return vertex::decodeNghsSparseSeq<symmetricVertex, W, F>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseSeq(uintE vtx_id, uintT o, F &f, G &g) {
    return vertex::decodeNghsSparseSeq<symmetricVertex, W, F>(this, vtx_id,
        getInNeighbors(), getInDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseBlock(uintE vtx_id, uintT o, uintE block_size, uintE block_num, F &f, G &g) {
    return vertex::decodeNghsSparseBlock<symmetricVertex, W, F>(
       this, vtx_id, getOutNeighbors(), getOutDegree(), o, block_size, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseBlock(uintE vtx_id, uintT o, uintE block_size, uintE block_num, F &f, G &g) {
    return vertex::decodeNghsSparseBlock<symmetricVertex, W, F>(
       this, vtx_id, getInNeighbors(), getInDegree(), o, block_size, block_num, f, g);
  }

  template <class F, class G>
  inline void copyOutNgh(uintE vtx_id, uintT o, F& f, G& g) {
    vertex::copyNghs<symmetricVertex, W>(this, vtx_id, getOutNeighbors(),
        getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline void copyInNgh(uintE vtx_id, uintT o, F& f, G& g) {
    vertex::copyNghs<symmetricVertex, W>(this, vtx_id, getInNeighbors(),
        getInDegree(), o, f, g);
  }

  inline tuple<uintE, W> get_ith_out_neighbor(uintE vtx_id, size_t i) {
    return vertex::get_ith_neighbor<symmetricVertex, W>(getOutNeighbors(), i);
  }

  inline tuple<uintE, W> get_ith_in_neighbor(uintE vtx_id, size_t i) {
    return vertex::get_ith_neighbor<symmetricVertex, W>(getInNeighbors(), i);
  }

  template <class F>
  inline size_t countOutNgh(uintE vtx_id, F &f, bool parallel=true) {
    return vertex::countNghs<symmetricVertex, W, F>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), f, parallel);
  }

  template <class F>
  inline size_t countInNgh(uintE vtx_id, F &f, bool parallel=true) {
    return vertex::countNghs<symmetricVertex, W, F>(this, vtx_id,
        getInNeighbors(), getInDegree(), f, parallel);
  }

  template <class E, class M, class R>
  inline E reduceOutNgh(uintE vtx_id, E id, M& m, R& r) {
    return vertex::reduceNghs<symmetricVertex, W, E, M, R>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), id, m, r);
  }

  template <class E, class M, class R>
  inline E reduceInNgh(uintE vtx_id, E id, M& m, R& r) {
    return vertex::reduceNghs<symmetricVertex, W, E, M, R>(this, vtx_id,
        getInNeighbors(), getInDegree(), id, m, r);
  }

  template <class F>
  inline void mapOutNgh(uintE vtx_id, F &f, bool parallel=true) {
    vertex::mapNghs<symmetricVertex, W, F>(this, vtx_id, getOutNeighbors(),
        getOutDegree(), f, parallel);
  }

  template <class F>
  inline void mapInNgh(uintE vtx_id, F &f, bool parallel=true) {
    vertex::mapNghs<symmetricVertex, W, F>(this, vtx_id, getInNeighbors(),
        getInDegree(), f, parallel);
  }

  template <class P, class O>
  inline void filterOutNgh(uintE vtx_id, P &p,
      O& out, tuple<uintE, W>* tmp) {
    vertex::filterNghs<symmetricVertex, W, P, O>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), p, out, tmp);
  }

  template <class P, class O>
  inline void filterInNgh(uintE vtx_id, P &p,
      O& out, tuple<uintE, W>* tmp) {
    vertex::filterNghs<symmetricVertex, W, P, O>(this, vtx_id,
        getInNeighbors(), getInDegree(), p, out, tmp);
  }

  template <class P>
  inline size_t packOutNgh(uintE vtx_id, P &p, tuple<uintE, W>* tmp) {
    uintE k = vertex::packNghs<symmetricVertex, W, P>(this, vtx_id, p,
        getOutNeighbors(), getOutDegree(), tmp);
    setOutDegree(k);
    return k;
  }

  template <class P>
  inline size_t packInNgh(uintE vtx_id, P &p, tuple<uintE, W>* tmp) {
    uintE k = vertex::packNghs<symmetricVertex, W, P>(this, vtx_id, p,
        getInNeighbors(), getInDegree(), tmp);
    setInDegree(k);
    return k;
  }

  inline size_t calculateOutTemporarySpace() {
    return vertex::calculateTemporarySpace(getOutDegree());
  }

  inline size_t calculateInTemporarySpace() {
    return vertex::calculateTemporarySpace(getInDegree());
  }

};

template <class W>
struct asymmetricVertex {
  tuple<uintE, W>* inNeighbors, *outNeighbors;
  uintE outDegree;
  uintE inDegree;
  void del() { free(inNeighbors); free(outNeighbors); }
asymmetricVertex(tuple<uintE, W>* iN, tuple<uintE, W>* oN, uintE id, uintE od)
: inNeighbors(iN), outNeighbors(oN), inDegree(id), outDegree(od) { }

  tuple<uintE, W>* getInNeighbors () { return inNeighbors; }
  tuple<uintE, W>* getOutNeighbors () { return outNeighbors; }
  uintE getInNeighbor(uintE j) { return get<0>(inNeighbors[j]); }
  uintE getOutNeighbor(uintE j) { return get<0>(outNeighbors[j]); }
  W getInWeight(uintE j) { return get<1>(inNeighbors[j]); }
  W getOutWeight(uintE j) { return get<1>(outNeighbors[j]); }
  void setInNeighbor(uintE j, uintE ngh) { get<0>(inNeighbors[j]) = ngh; }
  void setOutNeighbor(uintE j, uintE ngh) { get<0>(outNeighbors[j]) = ngh; }
  void setInWeight(uintE j, W wgh) { get<1>(inNeighbors[j]) = wgh; }
  void setOutWeight(uintE j, W wgh) { get<1>(outNeighbors[j]) = wgh; }
  void setInNeighbors(tuple<uintE, W>* _i) { inNeighbors = _i; }
  void setOutNeighbors(tuple<uintE, W>* _i) { outNeighbors = _i; }

  uintE getInDegree() { return inDegree; }
  uintE getOutDegree() { return outDegree; }
  uintE getInVirtualDegree() { return inDegree; }
  uintE getOutVirtualDegree() { return outDegree; }

  void setInDegree(uintE _d) { inDegree = _d; }
  void setOutDegree(uintE _d) { outDegree = _d; }
  void flipEdges() { swap(inNeighbors,outNeighbors); swap(inDegree,outDegree); }

  auto getInIter(uintE id) { return vertex::get_iter(getInNeighbors(), getInDegree()); }
  auto getOutIter(uintE id) { return vertex::get_iter(getOutNeighbors(), getOutDegree()); }

  inline size_t intersect(asymmetricVertex<W>* other, long our_id, long other_id) {
    return intersection::intersect(this, other, our_id, other_id);
  }

  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(uintE vtx_id, VS& vertexSubset, F &f, G &g,
      bool parallel = 0) {
    vertex::decodeNghsBreakEarly<asymmetricVertex, W, F, G, VS>(this, vtx_id,
        getInNeighbors(), getInDegree(), vertexSubset, f, g, parallel);
  }

  template <class VS, class F, class G>
  inline void decodeOutNghBreakEarly(uintE vtx_id, VS& vertexSubset, F &f, G &g,
      bool parallel = 0) {
    vertex::decodeNghsBreakEarly<asymmetricVertex, W, F, G, VS>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), vertexSubset, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeOutNgh(uintE vtx_id, F &f, G& g) {
     vertex::decodeNghs<asymmetricVertex, W, F, G>(this, vtx_id,
         getOutNeighbors(), getOutDegree(), f, g);
  }

  template <class F, class G>
  inline void decodeInNgh(uintE vtx_id, F &f, G& g) {
     vertex::decodeNghs<asymmetricVertex, W, F, G>(this, vtx_id,
         getInNeighbors(), getInDegree(), f, g);
  }

  template <class F, class G>
  inline void decodeOutNghSparse(uintE vtx_id, uintT o, F &f, G &g) {
    vertex::decodeNghsSparse<asymmetricVertex, W, F>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline void decodeInNghSparse(uintE vtx_id, uintT o, F &f, G &g) {
    vertex::decodeNghsSparse<asymmetricVertex, W, F>(this, vtx_id,
        getInNeighbors(), getInDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseSeq(uintE vtx_id, uintT o, F &f, G &g) {
    return vertex::decodeNghsSparseSeq<asymmetricVertex, W, F>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseSeq(uintE vtx_id, uintT o, F &f, G &g) {
    return vertex::decodeNghsSparseSeq<asymmetricVertex, W, F>(this, vtx_id,
        getInNeighbors(), getInDegree(), o, f, g);
  }

  template <class F, class G>
  inline size_t decodeOutNghSparseBlock(uintE vtx_id, uintT o, uintE block_size, uintE block_num, F &f, G &g) {
    return vertex::decodeNghsSparseBlock<asymmetricVertex, W, F>(
       this, vtx_id, getOutNeighbors(), getOutDegree(), o, block_size, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInNghSparseBlock(uintE vtx_id, uintT o, uintE block_size, uintE block_num, F &f, G &g) {
    return vertex::decodeNghsSparseBlock<asymmetricVertex, W, F>(
       this, vtx_id, getInNeighbors(), getInDegree(), o, block_size, block_num, f, g);
  }

  template <class F, class G>
  inline void copyOutNgh(uintE vtx_id, uintT o, F& f, G& g) {
    vertex::copyNghs<asymmetricVertex, W>(this, vtx_id, getOutNeighbors(),
        getOutDegree(), o, f, g);
  }

  template <class F, class G>
  inline void copyInNgh(uintE vtx_id, uintT o, F& f, G& g) {
    vertex::copyNghs<asymmetricVertex, W>(this, vtx_id, getInNeighbors(),
        getInDegree(), o, f, g);
  }

  inline tuple<uintE, W> get_ith_out_neighbor(uintE vtx_id, size_t i) {
    return vertex::get_ith_neighbor<asymmetricVertex, W>(getOutNeighbors(), i);
  }

  inline tuple<uintE, W> get_ith_in_neighbor(uintE vtx_id, size_t i) {
    return vertex::get_ith_neighbor<asymmetricVertex, W>(getInNeighbors(), i);
  }

  template <class F>
  inline size_t countOutNgh(uintE vtx_id, F &f, bool parallel=true) {
    return vertex::countNghs<asymmetricVertex, W, F>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), f, parallel);
  }

  template <class F>
  inline size_t countInNgh(uintE vtx_id, F &f, bool parallel=true) {
    return vertex::countNghs<asymmetricVertex, W, F>(this, vtx_id,
        getInNeighbors(), getInDegree(), f, parallel);
  }

  template <class E, class M, class R>
  inline E reduceOutNgh(uintE vtx_id, E id, M& m, R& r) {
    return vertex::reduceNghs<asymmetricVertex, W, E, M, R>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), id, m, r);
  }

  template <class E, class M, class R>
  inline E reduceInNgh(uintE vtx_id, E id, M& m, R& r) {
    return vertex::reduceNghs<asymmetricVertex, W, E, M, R>(this, vtx_id,
        getInNeighbors(), getInDegree(), id, m, r);
  }

  template <class F>
  inline void mapOutNgh(uintE vtx_id, F &f, bool parallel=true) {
    vertex::mapNghs<asymmetricVertex, W, F>(this, vtx_id, getOutNeighbors(),
        getOutDegree(), f, parallel);
  }

  template <class F>
  inline void mapInNgh(uintE vtx_id, F &f, bool parallel=true) {
    vertex::mapNghs<asymmetricVertex, W, F>(this, vtx_id, getInNeighbors(),
        getInDegree(), f, parallel);
  }

  template <class P, class O>
  inline void filterOutNgh(uintE vtx_id, P &p,
      O& out, tuple<uintE, W>* tmp) {
    vertex::filterNghs<asymmetricVertex, W, P, O>(this, vtx_id,
        getOutNeighbors(), getOutDegree(), p, out, tmp);
  }

  template <class P, class O>
  inline void filterInNgh(uintE vtx_id, P &p,
      O& out, tuple<uintE, W>* tmp) {
    vertex::filterNghs<asymmetricVertex, W, P, O>(this, vtx_id,
        getInNeighbors(), getInDegree(), p, out, tmp);
  }

  template <class P>
  inline size_t packOutNgh(uintE vtx_id, P &p, tuple<uintE, W>* tmp) {
    uintE k = vertex::packNghs<asymmetricVertex, W, P>(this, vtx_id, p,
        getOutNeighbors(), getOutDegree(), tmp);
    setOutDegree(k);
    return k;
  }

  template <class P>
  inline size_t packInNgh(uintE vtx_id, P &p, tuple<uintE, W>* tmp) {
    uintE k = vertex::packNghs<asymmetricVertex, W, P>(this, vtx_id, p,
        getInNeighbors(), getInDegree(), tmp);
    setInDegree(k);
    return k;
  }

  inline size_t calculateOutTemporarySpace() {
    return vertex::calculateTemporarySpace(getOutDegree());
  }

  inline size_t calculateInTemporarySpace() {
    return vertex::calculateTemporarySpace(getInDegree());
  }
};
