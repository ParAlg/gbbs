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
#include "block_managers.h"

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

}  // namespace block_vertex_ops


/*
 * block_manager supplies:
 *  - num_blocks
 *  - decode_block, decode_block_cond
 */
template <class W, /* weight type */
          class BM /* block manager type */>
struct block_symmetric_vertex {

  BM block_manager; // copy; not a reference.

  block_symmetric_vertex(BM&& block_manager) :
    block_manager(std::move(block_manager)) {}

  uintE getOutDegree() {
    return block_manager.get_degree();
  }
  uintE getInDegree() {
    return getOutDegree();
  }

  uintE getNumInBlocks() { return block_manager.num_blocks(); }
  uintE getNumOutBlocks() { return getNumInBlocks(); }
  inline uintE in_block_degree(uintE block_num) {
    return block_manager.block_degree(block_num);
  }
  inline uintE out_block_degree(uintE block_num) {
    return in_block_degree(block_num);
  }

  template <class F>
  inline void mapOutNgh(uintE vtx_id, F& f, bool parallel = true) {
    block_vertex_ops::map_nghs<BM, W, F>(vtx_id, block_manager, f, parallel);
  }

  template <class F>
  inline void mapInNgh(uintE vtx_id, F& f, bool parallel = true) {
    return mapOutNgh(vtx_id, f, parallel);
  }

  template <class M, class Monoid>
  inline auto reduceOutNgh(uintE vtx_id, M& m, Monoid& reduce,
                           bool parallel=true) -> typename Monoid::T {
    return block_vertex_ops::map_reduce<BM, W, M, Monoid>(
        vtx_id, block_manager, m, reduce, parallel);
  }

  template <class M, class Monoid>
  inline auto reduceInNgh(uintE vtx_id, M& m, Monoid& reduce,
                          bool parallel=true) -> typename Monoid::T {
    return reduceOutNgh(vtx_id, m, reduce, parallel);
  }

  template <class F>
  inline size_t countOutNgh(uintE vtx_id, F& f, bool parallel = true) {
    auto reduce = pbbs::addm<size_t>();
    return reduceOutNgh(vtx_id, f, reduce, parallel);
  }

  template <class F>
  inline size_t countInNgh(uintE vtx_id, F& f, bool parallel = true) {
    return countOutNgh(vtx_id, f, parallel);
  }

  inline std::tuple<uintE, W> get_ith_out_neighbor(uintE vtx_id, size_t i) {
    return block_manager.ith_neighbor(i);
  }

  inline std::tuple<uintE, W> get_ith_in_neighbor(uintE vtx_id, size_t i) {
    return block_manager.ith_neighbor(i);
  }

  /* copying primitives */
  template <class F, class G>
  inline void copyOutNgh(uintE vtx_id, uintT o, F& f, G& g, bool parallel=true) {
    block_vertex_ops::copyNghs<BM, W>(vtx_id, block_manager, o, f, g, parallel);
  }

  template <class F, class G>
  inline void copyInNgh(uintE vtx_id, uintT o, F& f, G& g, bool parallel=true) {
    copyOutNgh(vtx_id, o, f, g, parallel);
  }

  /* edgemap primitives */
  template <class F, class G, class H>
  inline void decodeOutNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    block_vertex_ops::decodeNghsSparse<BM, W, F>(vtx_id, block_manager, o, f, g, h, parallel);
  }

  template <class F, class G, class H>
  inline void decodeInNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h, bool parallel=true) {
    decodeOutNghSparse(vtx_id, o, f, g, h, parallel);
  }

  template <class F, class G>
  inline size_t decodeOutBlock(uintE vtx_id, uintT o, uintE block_num,
                               F& f, G& g) {
    return block_vertex_ops::decode_block<BM, W, F, G>(vtx_id,
        block_manager, o, block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInBlock(uintE vtx_id, uintT o, uintE block_num,
                              F& f, G& g) {
    return decodeOutBlock(vtx_id, o, block_num, f, g);
  }

  template <class F, class G>
  inline void decodeOutNgh(uintE vtx_id, F& f, G& g, bool parallel=true) {
    block_vertex_ops::decodeNghs<BM, W, F, G>(vtx_id,
        block_manager, f, g, parallel);
  }

  template <class F, class G>
  inline void decodeInNgh(uintE vtx_id, F& f, G& g, bool parallel=true) {
    decodeOutNgh(vtx_id, f, g, parallel);
  }

  template <class VS, class F, class G>
  inline void decodeOutNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g,
                                     bool parallel = false) {
    block_vertex_ops::decodeNghsBreakEarly<BM, W, F, G, VS>(
        vtx_id, block_manager, vertexSubset, f, g, parallel);
  }

  template <class VS, class F, class G>
  inline void decodeInNghBreakEarly(uintE vtx_id, VS& vertexSubset, F& f, G& g,
                                    bool parallel = false) {
    decodeOutNghBreakEarly(vtx_id, vertexSubset, f, g, parallel);
  }

};

