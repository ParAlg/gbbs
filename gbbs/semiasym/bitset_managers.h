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

#include "bitset.h"
#include "utils.h"

#include "gbbs/encodings/byte_pd_amortized.h"
#include "gbbs/flags.h"
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

/* Used to map over the edges incident to v */
template <class W /* weight */, class F /* user-specified mapping function */,
          class BM /* block_manager */>
inline void map_nghs(uintE vtx_id, BM& block_manager, F& f, bool parallel) {
  parallel_for(0, block_manager.num_blocks(),
               [&](size_t block_num) {
                 block_manager.decode_block(
                     block_num, [&](const uintE& ngh, const W& wgh,
                                    uintE edge_num) { f(vtx_id, ngh, wgh); });
               },
               1);
}

/* Map over edges incident to v using M, reduce using Monoid */
template <class W /* weight */, class M /* mapping function */,
          class Monoid /* reduction monoid */, class BM /* block_manager */>
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
      alloc = parlay::sequence<T>::uninitialized(num_blocks);
      block_outputs = alloc.begin();
    }

    parallel_for(0, num_blocks,
                 [&](size_t block_num) {
                   T cur = reduce.identity;
                   block_manager.decode_block(
                       block_num,
                       [&](const uintE& ngh, const W& wgh, uintE edge_num) {
                         cur = reduce.f(cur, m(vtx_id, ngh, wgh));
                       });
                   block_outputs[block_num] = cur;
                 },
                 1);

    auto im = gbbs::make_slice(block_outputs, num_blocks);
    T res = parlay::reduce(im, reduce);
    return res;
  } else {
    return reduce.identity;
  }
}

/* Used to copy edges incident to a vertex */
template <class W /* weight */, class F /* mapping function */,
          class G /* output function */, class BM /* block_manager */>
inline void copyNghs(uintE vtx_id, BM& block_manager, uintT o, F& f, G& g,
                     bool parallel) {
  parallel_for(0, block_manager.num_blocks(),
               [&](size_t block_num) {
                 block_manager.decode_block(
                     block_num,
                     [&](const uintE& ngh, const W& wgh, uintE edge_num) {
                       auto val = f(vtx_id, ngh, wgh);
                       g(ngh, o + edge_num, val);
                     });
               },
               1);
}

/* For each out-neighbor satisfying cond, call updateAtomic */
template <class W /* weight */, class F /* mapping function */,
          class G /* full output function */,
          class H /* empty output function */, class BM /* block_manager */>
inline void decodeNghsSparse(uintE vtx_id, BM& block_manager, uintT o, F& f,
                             G& g, H& h, bool parallel) {
  parallel_for(0, block_manager.num_blocks(),
               [&](size_t block_num) {
                 block_manager.decode_block(
                     block_num,
                     [&](const uintE& ngh, const W& wgh, uintE edge_num) {
                       if (f.cond(ngh)) {
                         auto m = f.updateAtomic(vtx_id, ngh, wgh);
                         g(ngh, o + edge_num, m);
                       } else {
                         h(ngh, o + edge_num);
                       }
                     });
               },
               1);
}

/* For each out-neighbor satisfying cond, call updateAtomic */
template <class W /* weight */, class F /* mapping function */,
          class G /* output function */, class BM /* block_manager */>
inline void decodeNghs(uintE vtx_id, BM& block_manager, F& f, G& g,
                       bool parallel) {
  parallel_for(0, block_manager.num_blocks(),
               [&](size_t block_num) {
                 block_manager.decode_block(
                     block_num,
                     [&](const uintE& ngh, const W& wgh, uintE edge_num) {
                       auto m = f.updateAtomic(vtx_id, ngh, wgh);
                       g(ngh, m);
                     });
               },
               1);
}

/* Sequentially process incident edges and quit if cond on self fails. */
template <class W /* weight */, class F /* mapping function */,
          class G /* output function */, class VS /* vertex_subset type */,
          class BM /* block_manager */>
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
      parallel_for(0, num_blocks,
                   [&](size_t block_num) {
                     block_manager.decode_block_cond(
                         block_num, [&](const uintE& ngh, const W& wgh,
                                        const size_t& edge_num) {
                           if (vertexSubset.isIn(ngh)) {
                             auto m = f.updateAtomic(ngh, vtx_id, wgh);
                             g(vtx_id, m);
                             return f.cond(vtx_id);
                           }
                           return true;
                         });
                   },
                   1);
    }
  }
}

template <class W /* weight */, class F /* edgemap struct */,
          class G /* output function */, class BM /* block manager */>
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

template <template <class W> class vertex, class W>
struct uncompressed_bitset_neighbors {
  using WV = vertex<W>;
  using E = typename WV::edge_type;
  using metadata = bitsets::metadata;

  uintE vtx_id;
  uintE vtx_degree;
  uintE vtx_original_degree;
  uintE vtx_num_blocks;
  uint8_t* blocks_start;
  uint8_t* block_data_start;

  vtx_info* v_infos;

  static constexpr uintE edges_per_block = 256;
  static constexpr uintE bytes_per_block =
      edges_per_block / 8 + sizeof(metadata);
  static constexpr uintE bitset_bytes_per_block =
      bytes_per_block - sizeof(metadata);
  static constexpr uintE kFullBlockPackThreshold = 2;
  static constexpr uintE kBlockAllocThreshold = 20;

  E* e0;

  uncompressed_bitset_neighbors(const uintE vtx_id, uint8_t* blocks,
                                uintE vtx_original_degree, vtx_info* v_infos,
                                E* e0)
      : vtx_id(vtx_id),
        vtx_original_degree(vtx_original_degree),
        v_infos(v_infos),
        e0(e0) {
    auto& v_info = v_infos[vtx_id];
    vtx_degree = v_info.vtx_degree;
    vtx_num_blocks = v_info.vtx_num_blocks;

    // TODO: ensure that blocks is 4-byte aligned.
    blocks_start = blocks + v_info.vtx_block_offset;
    block_data_start = blocks_start + (vtx_num_blocks * sizeof(metadata));
  }

  __attribute__((always_inline)) inline E* get_edges() { return e0; }

  __attribute__((always_inline)) inline uintE get_degree() {
    return vtx_degree;
  }

  __attribute__((always_inline)) inline size_t num_blocks() {
    return vtx_num_blocks;
  }

  __attribute__((always_inline)) inline size_t get_num_blocks() {
    return vtx_num_blocks;
  }

  //  uintE vtx_num_blocks;     // number of blocks associated with v
  //  size_t vtx_block_offset;  // pointer into the block structure
  //
  //  E* vtx_edges;  // pointer to the original edges (prevents one random read)

  inline void clear_vertex() {
    vtx_degree = 0;
    vtx_num_blocks = 0;
    v_infos[vtx_id].vtx_num_blocks = 0;
    v_infos[vtx_id].vtx_degree = 0;
  }

  __attribute__((always_inline)) inline uintE block_degree(uintE block_id) {
    return bitsets::block_degree(blocks_start, block_id, vtx_num_blocks,
                                 vtx_degree);
  }

  template <class P>
  inline size_t pack(P& p, uint8_t* tmp, bool parallel = true,
                     const flags fl = 0) {
    return pack_blocks(vtx_id, p, tmp, parallel, fl);
  }

  inline size_t calculateTemporarySpaceBytes() {
    if (vtx_degree > 0) {
      size_t nblocks = vtx_num_blocks;
      if (nblocks > kBlockAllocThreshold) {
        return (sizeof(uintE) * nblocks) + (bytes_per_block * nblocks);
      }
    }
    return 0;
  }

  /* mapping/reducing/counting primitives */
  template <class F>
  inline void map(F& f, bool parallel = true) {
    block_vertex_ops::map_nghs<W, F>(vtx_id, *this, f, parallel);
  }

  template <class M, class Monoid>
  inline auto reduce(M& m, Monoid& reduce, bool parallel = true) ->
      typename Monoid::T {
    return block_vertex_ops::map_reduce<W, M, Monoid>(vtx_id, *this, m, reduce,
                                                      parallel);
  }

  template <class F>
  inline size_t count(F& f, bool parallel = true) {
    auto reduce_f = parlay::addm<size_t>();
    return reduce(f, reduce_f, parallel);
  }

  /* edgemap primitives */
  template <class F, class G, class H>
  inline void decodeSparse(uintT o, F& f, G& g, H& h, bool parallel = true) {
    block_vertex_ops::decodeNghsSparse<W, F>(vtx_id, *this, o, f, g, h,
                                             parallel);
  }

  template <class F, class G>
  inline size_t decode_block(uintT o, uintE block_num, F& f, G& g) {
    return block_vertex_ops::decode_block<W, F, G>(vtx_id, *this, o, block_num,
                                                   f, g);
  }

  template <class F, class G>
  inline void decode(F& f, G& g, bool parallel = true) {
    block_vertex_ops::decodeNghs<W, F, G>(vtx_id, *this, f, g, parallel);
  }

  template <class VS, class F, class G>
  inline void decodeBreakEarly(VS& vertexSubset, F& f, G& g,
                               bool parallel = false) {
    block_vertex_ops::decodeNghsBreakEarly<W, F, G, VS>(
        vtx_id, *this, vertexSubset, f, g, parallel);
  }

  /* ===================== Implementations ====================== */

  template <class F>
  __attribute__((always_inline)) inline void decode_block(uintE block_id, F f) {
    metadata* block_metadata = (metadata*)blocks_start;
    uintE offset = block_metadata[block_id].offset;
    uintE next_block_offset = (block_id == vtx_num_blocks - 1)
                                  ? vtx_degree
                                  : block_metadata[block_id + 1].offset;
    uintE block_degree = next_block_offset - offset;
    if (block_degree == 0) return;

    uintE orig_block_num = block_metadata[block_id].block_num;

    uint8_t* block_bits = block_data_start + bitset_bytes_per_block * block_id;

    uintE block_start = orig_block_num * edges_per_block;
    uintE block_end =
        std::min(block_start + edges_per_block, vtx_original_degree);
    assert(block_start < block_end);

    E* e = get_edges();

    size_t block_size = block_end - block_start;
    size_t block_size_num_longs = (block_size + 64 - 1) / 64;

    // "select" on a block
    uint64_t* long_block_bits = (uint64_t*)block_bits;
    size_t cur_offset = block_start;
    for (size_t idx = 0; idx < block_size_num_longs; idx++) {
      uint64_t cur_long = long_block_bits[idx];
      while (cur_long > 0) {
        unsigned select_idx =
            _tzcnt_u64(cur_long);  // #trailing zeros in cur_long
        auto& ee = e[cur_offset + select_idx];
        f(std::get<0>(ee), std::get<1>(ee), offset++);
        assert((cur_long & (1UL << select_idx)) > 0);
        cur_long = _blsr_u64(cur_long);  // clears lowest set bit
      }
      cur_offset += 64;  // next long
    }
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block_cond(uintE block_id,
                                                               F f) {
    metadata* block_metadata = (metadata*)blocks_start;
    uintE offset = block_metadata[block_id].offset;
    uintE next_block_offset = (block_id == vtx_num_blocks - 1)
                                  ? vtx_degree
                                  : block_metadata[block_id + 1].offset;
    uintE block_degree = next_block_offset - offset;
    if (block_degree == 0) {
      return;
    }

    uintE orig_block_num = block_metadata[block_id].block_num;

    uint8_t* block_bits = block_data_start + bitset_bytes_per_block * block_id;

    uintE block_start = orig_block_num * edges_per_block;
    uintE block_end =
        std::min(block_start + edges_per_block, vtx_original_degree);

    E* e = get_edges();

    size_t block_size = block_end - block_start;
    size_t block_size_num_longs = (block_size + 64 - 1) / 64;

    // "select" on a block
    uint64_t* long_block_bits = (uint64_t*)block_bits;
    size_t cur_offset = block_start;
    for (size_t idx = 0; idx < block_size_num_longs; idx++) {
      uint64_t cur_long = long_block_bits[idx];
      while (cur_long > 0) {
        unsigned select_idx = _tzcnt_u64(cur_long);  // index of first nz bit
        auto& ee = e[cur_offset + select_idx];
        if (!f(std::get<0>(ee), std::get<1>(ee), offset++)) {
          return;
        }
        assert((cur_long & (1UL << select_idx)) > 0);
        cur_long = _blsr_u64(cur_long);  // reset lowest bit
      }
      cur_offset += 64;  // next long
    }
  }

  /* Only called when the discrepency between full and total blocks is large.
   * Specifically, when #full_blocks*kFullBlockPackThreshold >= vtx_num_blocks
   * Note: Defers to pack(..) to finish updating the degree/scan info */
  inline void repack_blocks_par(uint8_t* tmp, bool parallel) {
    uint8_t stk[bytes_per_block * kBlockAllocThreshold];  // temporary space
    uintE int_stk[kBlockAllocThreshold];
    uint8_t* tmp_space = (uint8_t*)stk;
    parlay::sequence<uint8_t> tmp_alloc;
    uintE* tmp_ints = int_stk;
    parlay::sequence<uintE> int_alloc;
    size_t total_bytes = vtx_num_blocks * bytes_per_block;
    if ((tmp == nullptr) && (vtx_num_blocks > kBlockAllocThreshold)) {
      tmp_alloc = parlay::sequence<uint8_t>::uninitialized(total_bytes);
      tmp_space = tmp_alloc.begin();
      int_alloc = parlay::sequence<uintE>::uninitialized(vtx_num_blocks);
      tmp_ints = int_alloc.begin();
    }

    // caller supplies:
    // vtx_num_blocks*sizeof(uintE) +
    // vtx_num_blocks*bytes_per_block

    if (tmp) {
      tmp_space = tmp;
    }

    metadata* block_metadata = (metadata*)blocks_start;

    size_t n_full_blocks = vtx_num_blocks - 1;  // full blocks
    size_t bytes_to_copy = n_full_blocks * bytes_per_block;
    {
      // fetch original block
      size_t last_block_num = vtx_num_blocks - 1;
      size_t orig_block_num = block_metadata[last_block_num].block_num;

      // get block size
      size_t block_start = orig_block_num * edges_per_block;
      size_t block_end = std::min(block_start + edges_per_block,
                                  static_cast<size_t>(vtx_original_degree));

      // #bytes for this block
      size_t last_block_size = block_end - block_start;
      size_t last_block_bytes =
          sizeof(metadata) +
          bitsets::get_bitset_block_size_in_bytes(last_block_size);

      bytes_to_copy += last_block_bytes;
    }
    assert(bytes_to_copy <= ((n_full_blocks + 1) * bytes_per_block));

    // Is a blocked memcpy faster here?
    parallel_for(0, bytes_to_copy,
                 [&](size_t i) { tmp_space[i] = blocks_start[i]; }, 512);

    // Tmp space for integers starts consecutively after tmp block space
    if (tmp) {
      tmp_ints = (uintE*)(tmp + bytes_to_copy);
    }

    // 2. Write 1 to tmp_ints (new_locs) if full, 0 if empty
    auto new_locs = gbbs::make_slice(tmp_ints, vtx_num_blocks);
    auto tmp_metadata = (metadata*)tmp_space;
    parallel_for(0, vtx_num_blocks, [&](size_t block_id) {
      new_locs[block_id] =
          static_cast<uintE>((tmp_metadata[block_id].offset > 0));
    });

    // 3. Scan new_locs to get new block indices for full blocks
    size_t new_num_blocks = parlay::scan_inplace(new_locs);

    // 4. Copy saved blocks to new positions.
    auto real_metadata = (metadata*)blocks_start;
    auto real_block_data = blocks_start + (new_num_blocks * sizeof(metadata));
    auto tmp_block_data = tmp_space + (vtx_num_blocks * sizeof(metadata));
    parallel_for(0, vtx_num_blocks, [&](size_t block_id) {
      uintE block_entries = tmp_metadata[block_id].offset;
      if (block_entries > 0) {  // live
        uintE new_block_id = new_locs[block_id];
        // (a) copy metadata
        real_metadata[new_block_id] = tmp_metadata[block_id];
        uintE orig_block_num = real_metadata[new_block_id].block_num;

        uint8_t* tmp_block_bits =
            tmp_block_data + bitset_bytes_per_block * block_id;
        uint8_t* real_block_bits =
            real_block_data + bitset_bytes_per_block * new_block_id;

        uintE block_start = orig_block_num * edges_per_block;
        uintE block_end =
            std::min(block_start + edges_per_block, vtx_original_degree);

        size_t this_block_size = block_end - block_start;
        size_t num_bytes_to_copy =
            bitsets::get_bitset_block_size_in_bytes(this_block_size);
        assert(num_bytes_to_copy <= bytes_per_block);

        // (b) copy bitset data
        for (size_t i = 0; i < num_bytes_to_copy; i++) {
          real_block_bits[i] = tmp_block_bits[i];
        }
      }
    });

    // 5. Update num_blocks info both locally and in v_infos
    vtx_num_blocks = new_num_blocks;
    v_infos[vtx_id].vtx_num_blocks = new_num_blocks;
  }

  // P : (uintE, uintE, wgh) -> bool
  // only keep edges s.t. P(...) = true.
  template <class P>
  inline size_t pack_blocks(uintE vtx_id, P& p, uint8_t* tmp, bool parallel,
                            const flags fl) {
    if (vtx_degree == 0) {
      return 0;
    }
    metadata* block_metadata = (metadata*)blocks_start;

    // 1. pack each block
    parallel_for(0, vtx_num_blocks,
                 [&](size_t block_id) {
                   uintE orig_block_num = block_metadata[block_id].block_num;
                   uint8_t* block_bits =
                       block_data_start + bitset_bytes_per_block * block_id;

                   uintE block_start = orig_block_num * edges_per_block;
                   uintE block_end = std::min(block_start + edges_per_block,
                                              vtx_original_degree);

                   E* e = get_edges();

                   size_t block_size = block_end - block_start;
                   size_t block_size_num_longs = (block_size + 64 - 1) / 64;

                   // "select" on a block
                   uint64_t* long_block_bits = (uint64_t*)block_bits;
                   size_t cur_offset = block_start;
                   size_t live_edges = 0;
                   for (size_t idx = 0; idx < block_size_num_longs; idx++) {
                     uint64_t cur_long = long_block_bits[idx];
                     // size_t cnt = _mm_popcnt_u64(cur_long); // #bits set to
                     // one
                     uint64_t long_to_write = cur_long;
                     while (cur_long > 0) {
                       unsigned select_idx = _tzcnt_u64(cur_long);
                       auto& ee = e[cur_offset + select_idx];
                       if (!p(vtx_id, std::get<0>(ee), std::get<1>(ee))) {
                         long_to_write ^= (1UL << select_idx);
                       } else {
                         live_edges++;
                       }
                       assert((cur_long & (1UL << select_idx)) > 0);
                       cur_long = _blsr_u64(cur_long);
                     }
                     long_block_bits[idx] = long_to_write;
                     cur_offset += 64;  // next long
                   }

                   // Temporarily store #live_edges in offset positions.
                   block_metadata[block_id].offset = live_edges;
                 },
                 1);

    // 2. Reduce to get the #empty_blocks
    auto full_block_seq =
        parlay::delayed_seq<size_t>(vtx_num_blocks, [&](size_t i) {
          return static_cast<size_t>(block_metadata[i].offset > 0);
        });
    size_t full_blocks = parlay::reduce(full_block_seq);

    if ((full_blocks * kFullBlockPackThreshold <= vtx_num_blocks) ||
        ((full_blocks < vtx_num_blocks) && (fl & compact_blocks))) {
      repack_blocks_par(tmp, parallel);
    }
    uintE sum = 0;
    for (size_t i = 0; i < vtx_num_blocks; i++) {
      uintE cur = block_metadata[i].offset;
      block_metadata[i].offset = sum;
      sum += cur;
    }
    vtx_degree = sum;

    // Update the degree in vtx_info.
    v_infos[vtx_id].vtx_degree = sum;
    assert(sum <= vtx_original_degree);
    return sum;  // return the new degree
  }

  std::tuple<uintE, W> ith_neighbor(size_t i) {
    metadata* block_metadata = (metadata*)blocks_start;
    auto offsets_imap = parlay::delayed_seq<size_t>(
        vtx_num_blocks, [&](size_t ind) { return block_metadata[ind].offset; });

    auto lte = [&](const size_t& l, const size_t& r) { return l <= r; };
    size_t block = parlay::binary_search(offsets_imap, i, lte);
    assert(block > 0);
    block = block - 1;

    std::tuple<uintE, W> out;
    auto decode_f = [&](const uintE& v, const W& wgh, const uintE& offset) {
      if (offset == i) {
        out = std::make_tuple(v, wgh);
        return false;
      }
      return true;
    };
    decode_block_cond(block, decode_f);
    return out;
  }

  struct iter {
    E* edges;
    uintE vtx_degree;
    uintE vtx_original_degree;
    uintE vtx_num_blocks;
    uint8_t* blocks_start;
    uint8_t* block_data_start;

    uintE cur_block;            // block id
    uintE cur_block_degree;     // #live edges
    uintE cur_block_start;      // start offset (in edges)
    uintE cur_block_num_longs;  // block size (in terms of #longs)

    uint64_t* cur_block_longs;
    uint64_t cur_long_value;
    uintE cur_long_idx;

    uintE last_ngh;
    uintE proc;
    uintE proc_cur_block;

    iter(E* edges, uintE vtx_degree, uintE vtx_original_degree,
         uintE vtx_num_blocks, uint8_t* blocks_start, uint8_t* block_data_start)
        : edges(edges),
          vtx_degree(vtx_degree),
          vtx_original_degree(vtx_original_degree),
          vtx_num_blocks(vtx_num_blocks),
          blocks_start(blocks_start),
          block_data_start(block_data_start),
          proc(0) {
      proc = 0;
      if (vtx_degree > 0) {
        next_nonempty_block(/* on_initialization = */ true);
        next();  // sets last_ngh
      }
    }

    // precondition: there is a subsequent non-empty block
    __attribute__((always_inline)) inline void next_nonempty_block(
        bool on_initialization = false) {
      if (on_initialization) {
        cur_block = 0;
      } else {
        cur_block++;
      }
      metadata* block_metadata = (metadata*)blocks_start;

      // set cur_block_degree
      uintE offset = block_metadata[cur_block].offset;
      uintE next_block_offset = (cur_block == vtx_num_blocks - 1)
                                    ? vtx_degree
                                    : block_metadata[cur_block + 1].offset;
      cur_block_degree = next_block_offset - offset;
      assert(cur_block_degree > 0);

      // reset proc_cur_block
      proc_cur_block = 0;

      // set cur_block_start and cur_block_size
      uintE orig_block_num = block_metadata[cur_block].block_num;
      cur_block_start = orig_block_num * edges_per_block;
      uintE cur_block_end =
          std::min(cur_block_start + edges_per_block, vtx_original_degree);
      uintE cur_block_size = cur_block_end - cur_block_start;
      cur_block_num_longs = (cur_block_size + 64 - 1) / 64;

      // set cur_block_longs
      cur_block_longs =
          (uint64_t*)(block_data_start + bitset_bytes_per_block * cur_block);

      // initialize value to read, and the idx.
      cur_long_value = cur_block_longs[0];
      cur_long_idx = 0;
    }

    __attribute__((always_inline)) inline uintE degree() { return vtx_degree; }
    __attribute__((always_inline)) inline uintE cur() { return last_ngh; }

    // updates last_ngh
    __attribute__((always_inline)) inline uintE next() {
      while (cur_long_value == 0) {
        // done with block?
        cur_long_idx++;
        proc_cur_block += 64;
        if (cur_long_idx == cur_block_num_longs) {
          next_nonempty_block();
        }
        cur_long_value = cur_block_longs[cur_long_idx];
      }
      // cur_long_value > 0

      unsigned select_idx =
          _tzcnt_u64(cur_long_value);  // #trailing zeros in cur_long
      cur_long_value = _blsr_u64(cur_long_value);  // clears lowest set bit
      last_ngh =
          std::get<0>(edges[cur_block_start + proc_cur_block + select_idx]);
      proc++;
      return last_ngh;
    }

    __attribute__((always_inline)) inline bool has_next() {
      return proc < vtx_degree;
    }
  };

  auto get_iter() {
    return iter(get_edges(), vtx_degree, vtx_original_degree, vtx_num_blocks,
                blocks_start, block_data_start);
  }
};

template <template <class W> class vertex, class W>
struct compressed_bitset_neighbors {
  using WV = vertex<W>;
  using E = typename WV::edge_type;
  using metadata = bitsets::metadata;

  uintE vtx_id;
  uintE vtx_degree;
  uintE vtx_original_degree;
  uintE vtx_num_blocks;
  uint8_t* blocks_start;
  uint8_t* block_data_start;
  static constexpr uintE kFullBlockPackThreshold = 4;
  static constexpr uintE kBlockAllocThreshold = 20;

  vtx_info* v_infos;

  static constexpr uintE edges_per_block = PARALLEL_DEGREE;
  static constexpr uintE bytes_per_block =
      edges_per_block / 8 + sizeof(metadata);
  static constexpr uintE bitset_bytes_per_block =
      bytes_per_block - sizeof(metadata);

  E* e0;

  compressed_bitset_neighbors(const uintE vtx_id, uint8_t* blocks,
                              uintE vtx_original_degree, vtx_info* v_infos,
                              E* e0)
      : vtx_id(vtx_id),
        vtx_original_degree(vtx_original_degree),
        v_infos(v_infos),
        e0(e0) {
    auto& v_info = v_infos[vtx_id];
    vtx_degree = v_info.vtx_degree;
    vtx_num_blocks = v_info.vtx_num_blocks;

    blocks_start = blocks + v_info.vtx_block_offset;
    block_data_start = blocks_start + (vtx_num_blocks * sizeof(metadata));
  }

  __attribute__((always_inline)) inline E* get_edges() { return e0; }

  __attribute__((always_inline)) inline uintE get_degree() {
    return vtx_degree;
  }

  inline void clear_vertex() {
    vtx_degree = 0;
    vtx_num_blocks = 0;
    v_infos[vtx_id].vtx_num_blocks = 0;
    v_infos[vtx_id].vtx_degree = 0;
  }

  __attribute__((always_inline)) inline size_t num_blocks() {
    return vtx_num_blocks;
  }

  __attribute__((always_inline)) inline size_t get_num_blocks() {
    return vtx_num_blocks;
  }

  __attribute__((always_inline)) inline uintE block_degree(uintE block_id) {
    return bitsets::block_degree(blocks_start, block_id, vtx_num_blocks,
                                 vtx_degree);
  }

  template <class P>
  inline size_t pack(P& p, uint8_t* tmp, bool parallel = true,
                     const flags fl = 0) {
    return pack_blocks(vtx_id, p, tmp, parallel, fl);
  }

  inline size_t calculateTemporarySpaceBytes() {
    if (vtx_degree > 0) {
      size_t nblocks = vtx_num_blocks;
      if (nblocks > kBlockAllocThreshold) {
        return (sizeof(uintE) * nblocks) + (bytes_per_block * nblocks);
      }
    }
    return 0;
  }

  /* mapping/reducing/counting primitives */
  template <class F>
  inline void map(F& f, bool parallel = true) {
    block_vertex_ops::map_nghs<W, F>(vtx_id, *this, f, parallel);
  }

  template <class M, class Monoid>
  inline auto reduce(M& m, Monoid& reduce, bool parallel = true) ->
      typename Monoid::T {
    return block_vertex_ops::map_reduce<W, M, Monoid>(vtx_id, *this, m, reduce,
                                                      parallel);
  }

  template <class F>
  inline size_t count(F& f, bool parallel = true) {
    auto reduce_f = parlay::addm<size_t>();
    return reduce(f, reduce_f, parallel);
  }

  /* edgemap primitives */
  template <class F, class G, class H>
  inline void decodeSparse(uintT o, F& f, G& g, H& h, bool parallel = true) {
    block_vertex_ops::decodeNghsSparse<W, F>(vtx_id, *this, o, f, g, h,
                                             parallel);
  }

  template <class F, class G>
  inline size_t decode_block(uintT o, uintE block_num, F& f, G& g) {
    return block_vertex_ops::decode_block<W, F, G>(vtx_id, *this, o, block_num,
                                                   f, g);
  }

  template <class F, class G>
  inline void decode(F& f, G& g, bool parallel = true) {
    block_vertex_ops::decodeNghs<W, F, G>(vtx_id, *this, f, g, parallel);
  }

  template <class VS, class F, class G>
  inline void decodeBreakEarly(VS& vertexSubset, F& f, G& g,
                               bool parallel = false) {
    block_vertex_ops::decodeNghsBreakEarly<W, F, G, VS>(
        vtx_id, *this, vertexSubset, f, g, parallel);
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block(uintE block_id, F f) {
    metadata* block_metadata = (metadata*)blocks_start;
    uintE offset = block_metadata[block_id].offset;
    uintE orig_block_num = block_metadata[block_id].block_num;

    E* e = get_edges();

    // 1. decode the compressed block
    size_t k = 0;
    std::tuple<uintE, W> block_decode[edges_per_block];
    auto t_f = [&](const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
      assert(k < edges_per_block);
      block_decode[k++] = std::make_tuple(ngh, wgh);
    };
    bytepd_amortized::template decode_block<W>(
        t_f, e, vtx_id, vtx_original_degree, orig_block_num);

    uint64_t* long_block_bits =
        (uint64_t*)(block_data_start + bitset_bytes_per_block * block_id);
    uintE block_start = orig_block_num * edges_per_block;
    size_t block_end =
        std::min(block_start + edges_per_block, vtx_original_degree);
    size_t block_size = block_end - block_start;
    size_t block_size_num_longs = (block_size + 64 - 1) / 64;

    // "select" on a block
    size_t cur_offset = 0;
    for (size_t idx = 0; idx < block_size_num_longs; idx++) {
      uint64_t cur_long = long_block_bits[idx];
      while (cur_long > 0) {
        unsigned select_idx =
            _tzcnt_u64(cur_long);  // #trailing zeros in cur_long
        auto& ee = block_decode[cur_offset + select_idx];
        f(std::get<0>(ee), std::get<1>(ee), offset++);
        assert((cur_long & (1UL << select_idx)) > 0);
        cur_long = _blsr_u64(cur_long);  // clears lowest set bit
      }
      cur_offset += 64;  // next long
    }
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block_cond(uintE block_id,
                                                               F f) {
    metadata* block_metadata = (metadata*)blocks_start;
    uintE offset = block_metadata[block_id].offset;
    uintE orig_block_num = block_metadata[block_id].block_num;

    E* e = get_edges();

    size_t k = 0;
    uint8_t* block_bits = block_data_start + bitset_bytes_per_block * block_id;
    auto t_f = [&](const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
      if (bitsets::is_bit_set(block_bits,  // about 7% overhead
                              k++)) {      // check if the k-th bit is set
        return f(ngh, wgh, offset++);  // and apply f with the correct offset
      }
      return true;
    };
    bytepd_amortized::template decode_block_cond<W>(
        t_f, e, vtx_id, vtx_original_degree, orig_block_num);
    //    debug(uintE block_start = orig_block_num * edges_per_block;
    //    uintE block_end =
    //        std::min(block_start + edges_per_block, vtx_original_degree);
    //    assert(k == (block_end - block_start));
    //    );
  }

  /* Only called when the discrepency between full and total blocks is large.
   * Specifically, when #full_blocks*kFullBlockPackThreshold >= vtx_num_blocks
   * Note: Defers to pack(..) to finish updating the degree/scan info */
  inline void repack_blocks_par(uint8_t* tmp, bool parallel) {
    uint8_t stk[bytes_per_block * kBlockAllocThreshold];  // temporary space
    uintE int_stk[kBlockAllocThreshold];

    uint8_t* tmp_space = (uint8_t*)stk;
    parlay::sequence<uint8_t> tmp_alloc;
    uintE* tmp_ints = (uintE*)int_stk;
    parlay::sequence<uintE> int_alloc;
    size_t total_bytes = vtx_num_blocks * bytes_per_block;
    if ((tmp == nullptr) && (vtx_num_blocks > kBlockAllocThreshold)) {
      tmp_alloc = parlay::sequence<uint8_t>::uninitialized(total_bytes);
      tmp_space = tmp_alloc.begin();
      int_alloc = parlay::sequence<uintE>::uninitialized(vtx_num_blocks);
      tmp_ints = int_alloc.begin();
    }

    // caller supplies:
    // vtx_num_blocks*sizeof(uintE) +
    // vtx_num_blocks*bytes_per_block

    if (tmp) {
      tmp_space = tmp;
    }

    metadata* block_metadata = (metadata*)blocks_start;

    size_t n_full_blocks = vtx_num_blocks - 1;  // full blocks
    size_t bytes_to_copy = n_full_blocks * bytes_per_block;
    {
      // fetch original block
      size_t last_block_num = vtx_num_blocks - 1;
      size_t orig_block_num = block_metadata[last_block_num].block_num;

      // get block size
      size_t block_start = orig_block_num * edges_per_block;
      size_t block_end = std::min(block_start + edges_per_block,
                                  static_cast<size_t>(vtx_original_degree));

      // #bytes for this block
      size_t last_block_size = block_end - block_start;
      size_t last_block_bytes =
          sizeof(metadata) +
          bitsets::get_bitset_block_size_in_bytes(last_block_size);

      bytes_to_copy += last_block_bytes;
    }

    // Is a blocked memcpy faster here?
    parallel_for(0, bytes_to_copy,
                 [&](size_t i) { tmp_space[i] = blocks_start[i]; }, 512);

    // Tmp space for integers starts consecutively after tmp block space
    if (tmp) {
      tmp_ints = (uintE*)(tmp + bytes_to_copy);
    }

    // 2. Write 1 to tmp_ints (new_locs) if full, 0 if empty
    auto new_locs = gbbs::make_slice(tmp_ints, vtx_num_blocks);
    auto tmp_metadata = (metadata*)tmp_space;
    parallel_for(0, vtx_num_blocks, [&](size_t block_id) {
      new_locs[block_id] =
          static_cast<uintE>((tmp_metadata[block_id].offset > 0));
    });

    // 3. Scan new_locs to get new block indices for full blocks
    size_t new_num_blocks = parlay::scan_inplace(new_locs);

    // 4. Copy saved blocks to new positions.
    auto real_metadata = (metadata*)blocks_start;
    auto real_block_data = blocks_start + (new_num_blocks * sizeof(metadata));
    auto tmp_block_data = tmp_space + (vtx_num_blocks * sizeof(metadata));
    parallel_for(0, vtx_num_blocks, [&](size_t block_id) {
      uintE block_entries = tmp_metadata[block_id].offset;
      if (block_entries > 0) {  // live
        uintE new_block_id = new_locs[block_id];
        // (a) copy metadata
        real_metadata[new_block_id] = tmp_metadata[block_id];
        uintE orig_block_num = real_metadata[new_block_id].block_num;

        uint8_t* tmp_block_bits =
            tmp_block_data + bitset_bytes_per_block * block_id;
        uint8_t* real_block_bits =
            real_block_data + bitset_bytes_per_block * new_block_id;

        uintE block_start = orig_block_num * edges_per_block;
        uintE block_end =
            std::min(block_start + edges_per_block, vtx_original_degree);

        size_t this_block_size = block_end - block_start;
        size_t block_bytes_to_copy =
            bitsets::get_bitset_block_size_in_bytes(this_block_size);

        // (b) copy bitset data
        for (size_t i = 0; i < block_bytes_to_copy; i++) {
          real_block_bits[i] = tmp_block_bits[i];
        }
      }
    });

    // 5. Update num_blocks info both locally and in v_infos
    vtx_num_blocks = new_num_blocks;
    v_infos[vtx_id].vtx_num_blocks = new_num_blocks;
  }

  // P : (uintE, uintE, wgh) -> bool
  // only keep edges s.t. P(...) = true.
  template <class P>
  inline size_t pack_blocks(uintE vtx_id, P& p, uint8_t* tmp, bool parallel,
                            const flags fl) {
    if (vtx_degree == 0) {
      return 0;
    }
    metadata* block_metadata = (metadata*)blocks_start;

    // 1. pack each block
    parallel_for(
        0, vtx_num_blocks,
        [&](size_t block_id) {
          uintE orig_block_num = block_metadata[block_id].block_num;
          E* e = get_edges();

          // (i) decode the block

          size_t k = 0;
          std::tuple<uintE, W> block_decode[edges_per_block];
          auto t_f = [&](const uintE& ngh, const W& wgh,
                         const uintE& orig_edge_id) {
            assert(k < edges_per_block);
            block_decode[k++] = std::make_tuple(ngh, wgh);
          };
          bytepd_amortized::template decode_block<W>(
              t_f, e, vtx_id, vtx_original_degree, orig_block_num);

          uint64_t* long_block_bits =
              (uint64_t*)(block_data_start + bitset_bytes_per_block * block_id);
          uintE block_start = orig_block_num * edges_per_block;
          size_t block_end =
              std::min(block_start + edges_per_block, vtx_original_degree);
          size_t block_size = block_end - block_start;
          size_t block_size_num_longs = (block_size + 64 - 1) / 64;

          // "select" on a block
          size_t cur_offset = 0;
          size_t live_edges = 0;
          for (size_t idx = 0; idx < block_size_num_longs; idx++) {
            uint64_t cur_long = long_block_bits[idx];
            uint64_t long_to_write = cur_long;
            while (cur_long > 0) {
              unsigned select_idx =
                  _tzcnt_u64(cur_long);  // #trailing zeros in cur_long
              auto& ee = block_decode[cur_offset + select_idx];
              if (!p(vtx_id, std::get<0>(ee), std::get<1>(ee))) {
                long_to_write ^= (1UL << select_idx);
              } else {
                live_edges++;
              }
              assert((cur_long & (1UL << select_idx)) > 0);
              cur_long = _blsr_u64(cur_long);  // clears lowest set bit
            }
            long_block_bits[idx] = long_to_write;
            cur_offset += 64;  // next long
          }
          // Temporarily store #live_edges in offset positions.
          block_metadata[block_id].offset = live_edges;
        },
        1);

    // 2. Reduce to get the #empty_blocks
    auto full_block_seq =
        parlay::delayed_seq<size_t>(vtx_num_blocks, [&](size_t i) {
          return static_cast<size_t>(block_metadata[i].offset > 0);
        });
    size_t full_blocks = parlay::reduce(full_block_seq);

    if ((full_blocks * kFullBlockPackThreshold <= vtx_num_blocks) ||
        ((full_blocks < vtx_num_blocks) && (fl & compact_blocks))) {
      repack_blocks_par(tmp, parallel);
    }

    uintE sum = 0;
    for (size_t i = 0; i < vtx_num_blocks; i++) {
      uintE cur = block_metadata[i].offset;
      block_metadata[i].offset = sum;
      sum += cur;
    }
    vtx_degree = sum;

    // Update the degree in vtx_info.
    v_infos[vtx_id].vtx_degree = sum;
    return sum;  // return the new degree
  }

  std::tuple<uintE, W> ith_neighbor(size_t i) {
    metadata* block_metadata = (metadata*)blocks_start;
    auto offsets_imap = parlay::delayed_seq<size_t>(
        vtx_num_blocks, [&](size_t ind) { return block_metadata[ind].offset; });

    auto lte = [&](const size_t& l, const size_t& r) { return l <= r; };
    size_t block = parlay::binary_search(offsets_imap, i, lte);
    assert(block > 0);
    block = block - 1;

    std::tuple<uintE, W> out;
    auto decode_f = [&](const uintE& v, const W& wgh, const uintE& offset) {
      if (offset == i) {
        out = std::make_tuple(v, wgh);
        return false;
      }
      return true;
    };
    decode_block_cond(block, decode_f);
    return out;
  }

  struct iter {
    E* edges;
    uintE vtx_id;
    uintE vtx_degree;
    uintE vtx_original_degree;
    uintE vtx_num_blocks;
    uint8_t* blocks_start;
    uint8_t* block_data_start;

    uintE cur_block;            // block id
    uintE cur_block_degree;     // #live edges
    uintE cur_block_start;      // start offset (in edges)
    uintE cur_block_num_longs;  // block size (in terms of #longs)

    uint64_t* cur_block_longs;
    uint64_t cur_long_value;
    uintE cur_long_idx;

    uintE block_decode[PARALLEL_DEGREE];

    uintE last_ngh;
    uintE proc;
    uintE proc_cur_block;

    iter(E* edges, uintE vtx_id, uintE vtx_degree, uintE vtx_original_degree,
         uintE vtx_num_blocks, uint8_t* blocks_start, uint8_t* block_data_start)
        : edges(edges),
          vtx_id(vtx_id),
          vtx_degree(vtx_degree),
          vtx_original_degree(vtx_original_degree),
          vtx_num_blocks(vtx_num_blocks),
          blocks_start(blocks_start),
          block_data_start(block_data_start),
          proc(0) {
      proc = 0;
      if (vtx_degree > 0) {
        next_nonempty_block(/* on_initialization = */ true);
        next();  // sets last_ngh
      }
    }

    // precondition: there is a subsequent non-empty block
    __attribute__((always_inline)) inline void next_nonempty_block(
        bool on_initialization = false) {
      if (on_initialization) {
        cur_block = 0;
      } else {
        cur_block++;
      }
      metadata* block_metadata = (metadata*)blocks_start;

      // set cur_block_degree
      uintE offset = block_metadata[cur_block].offset;
      uintE next_block_offset = (cur_block == vtx_num_blocks - 1)
                                    ? vtx_degree
                                    : block_metadata[cur_block + 1].offset;
      cur_block_degree = next_block_offset - offset;
      assert(cur_block_degree > 0);

      // reset proc_cur_block
      proc_cur_block = 0;

      // set cur_block_start and cur_block_size
      uintE orig_block_num = block_metadata[cur_block].block_num;
      cur_block_start = orig_block_num * edges_per_block;
      uintE cur_block_end =
          std::min(cur_block_start + edges_per_block, vtx_original_degree);
      uintE cur_block_size = cur_block_end - cur_block_start;
      cur_block_num_longs = (cur_block_size + 64 - 1) / 64;

      // set cur_block_longs
      cur_block_longs =
          (uint64_t*)(block_data_start + bitset_bytes_per_block * cur_block);

      // initialize value to read, and the idx.
      cur_long_value = cur_block_longs[0];
      cur_long_idx = 0;

      // decode current compressed block into block_decode
      size_t k = 0;
      auto map_f = [&](const uintE& v, const W& wgh, size_t off) {
        block_decode[k++] = v;
      };
      bytepd_amortized::template decode_block<W>(
          map_f, edges, vtx_id, vtx_original_degree, orig_block_num);
    }

    __attribute__((always_inline)) inline uintE degree() { return vtx_degree; }
    __attribute__((always_inline)) inline uintE cur() { return last_ngh; }

    // updates last_ngh
    __attribute__((always_inline)) inline uintE next() {
      while (cur_long_value == 0) {
        // done with block?
        cur_long_idx++;
        proc_cur_block += 64;
        if (cur_long_idx == cur_block_num_longs) {
          next_nonempty_block();
        }
        cur_long_value = cur_block_longs[cur_long_idx];
      }
      // cur_long_value > 0

      unsigned select_idx =
          _tzcnt_u64(cur_long_value);  // #trailing zeros in cur_long
      cur_long_value = _blsr_u64(cur_long_value);  // clears lowest set bit

      last_ngh = block_decode[proc_cur_block + select_idx];
      proc++;
      return last_ngh;
    }

    __attribute__((always_inline)) inline bool has_next() {
      return proc < vtx_degree;
    }
  };

  auto get_iter() {
    return iter(get_edges(), vtx_id, vtx_degree, vtx_original_degree,
                vtx_num_blocks, blocks_start, block_data_start);
  }
};

}  // namespace gbbs
