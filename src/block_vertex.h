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
#include "noop_managers.h"
#include "bitset_managers.h"

namespace block_vertex_ops {

/* Used to map over the edges incident to v */
template <class BM /* block_manager */,
          class W  /* weight */,
          class F  /* user-specified mapping function */>
inline void map_nghs(uintE vtx_id, BM& block_manager, F& f, bool parallel) {
  par_for(0, block_manager.num_blocks(), 1, [&] (size_t block_num) {
    block_manager.decode_block(block_num,
      [&] (const uintE& ngh, const W& wgh, uintE edge_num) {
        f(vtx_id, ngh, wgh);
      }
    );
  }, parallel);
}

/* Map over edges incident to v using M, reduce using Monoid */
template <class BM     /* block_manager */,
          class W      /* weight */,
          class M      /* mapping function */,
          class Monoid /* reduction monoid */>
inline auto map_reduce(uintE vtx_id, BM& block_manager, M& m,
           Monoid& reduce, bool parallel=true) -> typename Monoid::T {
  using T = typename Monoid::T;
  size_t num_blocks = block_manager.num_blocks();
  if (num_blocks >= 1) {
    T stk[100];
    T* block_outputs;
    if (num_blocks > 100) {
      // TODO: should the interface for reduce_nghs expect a tmp memory
      // allocation for sizes larger than 100 (reduce_alloc_thresh)? Concern is
      // that in-line memory allocation for large vertices is far slower than
      // doing a bulk-memory allocation up front and receiving offsets.
      block_outputs = pbbslib::new_array_no_init<T>(num_blocks);
    } else {
      block_outputs = (T*)stk;
    }

    par_for(0, num_blocks, 1, [&] (size_t block_num) {
      T cur = reduce.identity;
      block_manager.decode_block(block_num,
        [&] (const uintE& ngh, const W& wgh, uintE edge_num) {
          cur = reduce.f(cur, m(vtx_id, ngh, wgh));
        }
      );
      block_outputs[block_num] = cur;
    }, parallel);

    auto im = pbbslib::make_sequence(block_outputs, num_blocks);
    T res = pbbslib::reduce(im, reduce);
    if (num_blocks > 100) {
      pbbslib::free_array(block_outputs);
    }
    return res;
  } else {
    return reduce.identity;
  }
}

/* Used to copy edges incident to a vertex */
template <class BM     /* block_manager */,
          class W      /* weight */,
          class F      /* mapping function */,
          class G      /* output function */>
inline void copyNghs(uintE vtx_id, BM& block_manager,
                     uintT o, F& f, G& g, bool parallel) {
  par_for(0, block_manager.num_blocks(), 1, [&] (size_t block_num) {
    block_manager.decode_block(block_num,
      [&] (const uintE& ngh, const W& wgh, uintE edge_num) {
        auto val = f(vtx_id, ngh, wgh);
        g(ngh, o + edge_num, val);
      }
    );
  }, parallel);
}

/* For each out-neighbor satisfying cond, call updateAtomic */
template <class BM     /* block_manager */,
          class W      /* weight */,
          class F      /* mapping function */,
          class G      /* full output function */,
          class H      /* empty output function */>
inline void decodeNghsSparse(uintE vtx_id, BM& block_manager, uintT o, F& f,
                             G& g, H& h, bool parallel) {
  par_for(0, block_manager.num_blocks(), 1, [&] (size_t block_num) {
    block_manager.decode_block(block_num,
      [&] (const uintE& ngh, const W& wgh, uintE edge_num) {
        if (f.cond(ngh)) {
          auto m = f.updateAtomic(vtx_id, ngh, wgh);
          g(ngh, o + edge_num, m);
        } else {
          h(ngh, o + edge_num);
        }
      }
    );
  }, parallel);
}

/* For each out-neighbor satisfying cond, call updateAtomic */
template <class BM     /* block_manager */,
          class W      /* weight */,
          class F      /* mapping function */,
          class G      /* output function */>
inline void decodeNghs(uintE vtx_id, BM& block_manager, F& f, G& g, bool parallel) {
  par_for(0, block_manager.num_blocks(), 1, [&] (size_t block_num) {
    block_manager.decode_block(block_num,
      [&] (const uintE& ngh, const W& wgh, uintE edge_num) {
        auto m = f.updateAtomic(vtx_id, ngh, wgh);
        g(ngh, m);
      }
    );
  }, parallel);
}

/* Sequentially process incident edges and quit if cond on self fails. */
template <class BM     /* block_manager */,
          class W      /* weight */,
          class F      /* mapping function */,
          class G      /* output function */,
          class VS     /* vertex_subset type */>
inline void decodeNghsBreakEarly(uintE vtx_id, BM& block_manager,
                                 VS& vertexSubset, F& f, G& g,
                                 bool parallel) {
  if (block_manager.get_degree() > 0) {
    if (!parallel) {
      for (size_t i=0; i<block_manager.num_blocks(); i++) {
        block_manager.decode_block_cond(i,
          [&] (const uintE& ngh, const W& wgh, const size_t& edge_num) {
            if (vertexSubset.isIn(ngh)) {
              auto m = f.update(ngh, vtx_id, wgh);
              g(vtx_id, m);
              return f.cond(vtx_id);
            }
            return true;
          }
        );
      }
    } else {
      size_t num_blocks = block_manager.num_blocks();
      par_for(0, num_blocks, 1, [&] (size_t block_num) {
        block_manager.decode_block_cond(block_num,
          [&] (const uintE& ngh, const W& wgh, const size_t& edge_num) {
            if (vertexSubset.isIn(ngh)) {
              auto m = f.updateAtomic(ngh, vtx_id, wgh);
              g(vtx_id, m);
              return f.cond(vtx_id);
            }
            return true;
          }
        );
      });
    }
  }
}

template <class BM /* block manager */,
          class W  /* weight */,
          class F  /* edgemap struct */,
          class G  /* output function */>
inline size_t decode_block(uintE vtx_id, BM& block_manager,
                           uintT o, uintE block_num, F& f, G& g) {
  size_t k = 0;
  block_manager.decode_block(block_num,
    [&] (const uintE& ngh, const W& wgh, uintE edge_num) {
      if (f.cond(ngh)) {
        auto m = f.updateAtomic(vtx_id, ngh, wgh);
        bool wrote = g(ngh, o + k, m);
        if (wrote) {
          k++;
        }
      }
    }
  );
  return k;
}

//// functions to pack
//// Caller is responsible for setting the degree on v.
//template <class BM, class W, class P, class E>
//inline size_t packNghs(uintE vtx_id, BM& block_manager, P& p,
//                       E* tmp) {
//  // 1. Pack out all live blocks
//  par_for(0, block_manager.num_blocks(), 1, [&] (size_t block_num) {
//    block_manager.pack_block(block_num, p);
//  }, parallel);
//
//  block_manager.pack_blocks(p);
//
//  if (d < vertex_ops::kAllocThreshold) {
//    uintE k = 0;
//    for (size_t i = 0; i < d; i++) {
//      auto nw = nghs[i];
//      uintE ngh = std::get<0>(nw);
//      W wgh = std::get<1>(nw);
//      if (p(vtx_id, ngh, wgh)) {
//        nghs[k++] = std::make_tuple(ngh, wgh);
//      }
//    }
//    return k;
//  } else {
//    // copy to tmp
//    par_for(0, d, pbbslib::kSequentialForThreshold, [&] (size_t i) { tmp[i] = nghs[i]; });
//    auto pc = [&](const std::tuple<uintE, W>& nw) {
//      return p(vtx_id, std::get<0>(nw), std::get<1>(nw));
//    };
//    size_t k = pbbslib::filterf(tmp, nghs, d, pc);
//    return k;
//  }
//}





}  // namespace block_vertex_ops
