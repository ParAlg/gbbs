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

#include "bitset_managers.h"
#include "gbbs/macros.h"

namespace gbbs {

namespace block_vertex_ops {

template <class It>
size_t intersect(It& a, It& b) {
  size_t i = 0;
  size_t j = 0;
  size_t nA = a.degree();
  size_t nB = b.degree();
  size_t ans = 0;
  bool advance_a = false;
  bool advance_b = false;
  while (i < nA && j < nB) {
    if (advance_a) {
      advance_a = false;
      a.next();
    }
    if (advance_b) {
      advance_b = false;
      b.next();
    }
    if (a.cur() == b.cur()) {
      advance_a = true;
      advance_b = true;
      i++;
      j++;
      ans++;
    } else if (a.cur() < b.cur()) {
      advance_a = true;
      i++;
    } else {
      advance_b = true;
      j++;
    }
  }
  return ans;
}

template <class S, class It>
size_t intersect_seq(S& a, It& b) {
  size_t i = 0;
  size_t j = 0;
  size_t nA = a.size();
  size_t nB = b.degree();
  size_t ans = 0;
  uintE b_cur;
  if (j < nB) b_cur = b.cur();
  while (i < nA && j < nB) {
    if (a[i] == b_cur) {
      i++;
      j++;
      ans++;
      if (j < nB) {
        b_cur = b.next();
      }
    } else if (a[i] < b_cur) {
      i++;
    } else {
      j++;
      if (j < nB) {
        b_cur = b.next();
      }
    }
  }
  return ans;
}

template <class SeqA, class SeqB>
size_t seq_merge(SeqA& A, SeqB& B) {
  using T = typename SeqA::value_type;
  size_t nA = A.size(), nB = B.size();
  size_t i = 0, j = 0;
  size_t ct = 0;
  while (i < nA && j < nB) {
    T& a = A[i];
    T& b = B[j];
    if (a == b) {
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

// template <class S, class It>
// size_t intersect_batch_seq(S& a, It& b) {
//  size_t b_size = 128;
//  uintE block[b_size];
//  size_t n_blocks = parlay::num_blocks(b.degree(), b_size);
//  size_t ans = 0;
//  size_t i=0;
//  size_t block_id=0;
//
//  while (i < nA && block_id < n_blocks) {
//    size_t block_start =
//
//  }
//
//  size_t j=0;
//  size_t nA = a.size(); size_t nB = b.degree();
//
//  while (i < nA && j < nB) {
//    if (a[i] == b_cur) {
//      i++; j++; ans++;
//      if (b.has_next()) b.next();
//      b_cur = std::get<0>(b.cur());
//    } else if (a[i] < std::get<0>(b.cur())) {
//      i++;
//    } else {
//      j++;
//      if (b.has_next()) b.next();
//      b_cur = std::get<0>(b.cur());
//    }
//  }
//  return ans;
//}

/* Used to map over the edges incident to v */
template <class BM /* block_manager */, class W /* weight */,
          class F /* user-specified mapping function */>
inline void map_nghs(uintE vtx_id, BM& block_manager, F& f, bool parallel) {
  parallel_for(0, block_manager.num_blocks(), 1, [&](size_t block_num) {
    block_manager.decode_block(
        block_num, [&](const uintE& ngh, const W& wgh, uintE edge_num) {
          f(vtx_id, ngh, wgh);
        });
  });
}

/* Map over edges incident to v using M, reduce using Monoid */
template <class BM /* block_manager */, class W /* weight */,
          class M /* mapping function */, class Monoid /* reduction monoid */>
inline auto map_reduce(uintE vtx_id, BM& block_manager, M& m, Monoid& reduce,
                       bool parallel = true) -> typename Monoid::T {
  using T = typename Monoid::T;
  size_t num_blocks = block_manager.num_blocks();
  if (num_blocks >= 1) {
    T stk[100];
    T* block_outputs = stk;
    parlay::sequence<T> alloc;
    if (num_blocks > 100) {
      // TODO: should the interface for reduce_nghs expect a tmp memory
      // allocation for sizes larger than 100 (reduce_alloc_thresh)? Concern is
      // that in-line memory allocation for large vertices is far slower than
      // doing a bulk-memory allocation up front and receiving offsets.
      alloc = parlay::sequence<T>(num_blocks);
      block_outputs = alloc.begin();
    } else {
      block_outputs = (T*)stk;
    }

    parallel_for(0, num_blocks, 1, [&](size_t block_num) {
      T cur = reduce.identity;
      block_manager.decode_block(
          block_num, [&](const uintE& ngh, const W& wgh, uintE edge_num) {
            cur = reduce.f(cur, m(vtx_id, ngh, wgh));
          });
      block_outputs[block_num] = cur;
    });

    auto im = gbbs::make_slice(block_outputs, num_blocks);
    T res = parlay::reduce(im, reduce);
    if (num_blocks > 100) {
      gbbs::free_array(block_outputs);
    }
    return res;
  } else {
    return reduce.identity;
  }
}

/* Used to copy edges incident to a vertex */
template <class BM /* block_manager */, class W /* weight */,
          class F /* mapping function */, class G /* output function */>
inline void copyNghs(uintE vtx_id, BM& block_manager, uintT o, F& f, G& g,
                     bool parallel) {
  parallel_for(0, block_manager.num_blocks(), 1, [&](size_t block_num) {
    block_manager.decode_block(
        block_num, [&](const uintE& ngh, const W& wgh, uintE edge_num) {
          auto val = f(vtx_id, ngh, wgh);
          g(ngh, o + edge_num, val);
        });
  });
}

/* For each out-neighbor satisfying cond, call updateAtomic */
template <class BM /* block_manager */, class W /* weight */,
          class F /* mapping function */, class G /* full output function */,
          class H /* empty output function */>
inline void decodeNghsSparse(uintE vtx_id, BM& block_manager, uintT o, F& f,
                             G& g, H& h, bool parallel) {
  parallel_for(0, block_manager.num_blocks(), 1, [&](size_t block_num) {
    block_manager.decode_block(
        block_num, [&](const uintE& ngh, const W& wgh, uintE edge_num) {
          if (f.cond(ngh)) {
            auto m = f.updateAtomic(vtx_id, ngh, wgh);
            g(ngh, o + edge_num, m);
          } else {
            h(ngh, o + edge_num);
          }
        });
  });
}

/* For each out-neighbor satisfying cond, call updateAtomic */
template <class BM /* block_manager */, class W /* weight */,
          class F /* mapping function */, class G /* output function */>
inline void decodeNghs(uintE vtx_id, BM& block_manager, F& f, G& g,
                       bool parallel) {
  parallel_for(0, block_manager.num_blocks(), 1, [&](size_t block_num) {
    block_manager.decode_block(
        block_num, [&](const uintE& ngh, const W& wgh, uintE edge_num) {
          auto m = f.updateAtomic(vtx_id, ngh, wgh);
          g(ngh, m);
        });
  });
}

/* Sequentially process incident edges and quit if cond on self fails. */
template <class BM /* block_manager */, class W /* weight */,
          class F /* mapping function */, class G /* output function */,
          class VS /* vertex_subset type */>
inline void decodeNghsBreakEarly(uintE vtx_id, BM& block_manager,
                                 VS& vertexSubset, F& f, G& g, bool parallel) {
  if (block_manager.get_degree() > 0) {
    if (!parallel) {
      for (size_t i = 0; i < block_manager.num_blocks(); i++) {
        block_manager.decode_block_cond(
            i, [&](const uintE& ngh, const W& wgh, const size_t& edge_num) {
              if (vertexSubset.isIn(ngh)) {
                auto m = f.update(ngh, vtx_id, wgh);
                g(vtx_id, m);
                return f.cond(vtx_id);
              }
              return true;
            });
      }
    } else {
      size_t num_blocks = block_manager.num_blocks();
      parallel_for(0, num_blocks, 1, [&](size_t block_num) {
        block_manager.decode_block_cond(
            block_num,
            [&](const uintE& ngh, const W& wgh, const size_t& edge_num) {
              if (vertexSubset.isIn(ngh)) {
                auto m = f.updateAtomic(ngh, vtx_id, wgh);
                g(vtx_id, m);
                return f.cond(vtx_id);
              }
              return true;
            });
      });
    }
  }
}

template <class BM /* block manager */, class W /* weight */,
          class F /* edgemap struct */, class G /* output function */>
inline size_t decode_block(uintE vtx_id, BM& block_manager, uintT o,
                           uintE block_num, F& f, G& g) {
  size_t k = 0;
  block_manager.decode_block(
      block_num, [&](const uintE& ngh, const W& wgh, uintE edge_num) {
        if (f.cond(ngh)) {
          auto m = f.updateAtomic(vtx_id, ngh, wgh);
          bool wrote = g(ngh, o + k, m);
          if (wrote) {
            k++;
          }
        }
      });
  return k;
}

}  // namespace block_vertex_ops
}  // namespace gbbs
