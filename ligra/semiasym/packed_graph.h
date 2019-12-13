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

/*
 * This file contains an implementation of a bit-packed graph, which is a
 * decremental graph structure that supports edge-deletions using a bit-packed
 * vector.
 */

#pragma once

#include <functional>
#include <tuple>
#include <type_traits>

#include "bitset_managers.h"
#include "block_vertex.h"
#include "graph.h"

/*
 * block_manager supplies:
 *  - num_blocks
 *  - decode_block, decode_block_cond
 */
template <class W, /* weight type */
          class BM /* block manager type */>
struct packed_symmetric_vertex {
  using iter_type = typename BM::iter;
  BM block_manager;  // copy; not a reference.

  packed_symmetric_vertex(BM&& block_manager)
      : block_manager(std::move(block_manager)) {}

  __attribute__((always_inline)) inline uintE getOutDegree() {
    return block_manager.get_degree();
  }
  __attribute__((always_inline)) inline uintE getInDegree() {
    return getOutDegree();
  }

  __attribute__((always_inline)) inline uintE getNumInBlocks() {
    return block_manager.num_blocks();
  }
  __attribute__((always_inline)) inline uintE getNumOutBlocks() {
    return getNumInBlocks();
  }
  __attribute__((always_inline)) inline uintE in_block_degree(uintE block_num) {
    return block_manager.block_degree(block_num);
  }
  __attribute__((always_inline)) inline uintE out_block_degree(
      uintE block_num) {
    return in_block_degree(block_num);
  }

  /* mapping/reducing/counting primitives */
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
                           bool parallel = true) -> typename Monoid::T {
    return block_vertex_ops::map_reduce<BM, W, M, Monoid>(vtx_id, block_manager,
                                                          m, reduce, parallel);
  }

  template <class M, class Monoid>
  inline auto reduceInNgh(uintE vtx_id, M& m, Monoid& reduce,
                          bool parallel = true) -> typename Monoid::T {
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
  inline void copyOutNgh(uintE vtx_id, uintT o, F& f, G& g,
                         bool parallel = true) {
    block_vertex_ops::copyNghs<BM, W>(vtx_id, block_manager, o, f, g, parallel);
  }

  template <class F, class G>
  inline void copyInNgh(uintE vtx_id, uintT o, F& f, G& g,
                        bool parallel = true) {
    copyOutNgh(vtx_id, o, f, g, parallel);
  }

  /* edgemap primitives */
  template <class F, class G, class H>
  inline void decodeOutNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h,
                                 bool parallel = true) {
    block_vertex_ops::decodeNghsSparse<BM, W, F>(vtx_id, block_manager, o, f, g,
                                                 h, parallel);
  }

  template <class F, class G, class H>
  inline void decodeInNghSparse(uintE vtx_id, uintT o, F& f, G& g, H& h,
                                bool parallel = true) {
    decodeOutNghSparse(vtx_id, o, f, g, h, parallel);
  }

  template <class F, class G>
  inline size_t decodeOutBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                               G& g) {
    return block_vertex_ops::decode_block<BM, W, F, G>(vtx_id, block_manager, o,
                                                       block_num, f, g);
  }

  template <class F, class G>
  inline size_t decodeInBlock(uintE vtx_id, uintT o, uintE block_num, F& f,
                              G& g) {
    return decodeOutBlock(vtx_id, o, block_num, f, g);
  }

  template <class F, class G>
  inline void decodeOutNgh(uintE vtx_id, F& f, G& g, bool parallel = true) {
    block_vertex_ops::decodeNghs<BM, W, F, G>(vtx_id, block_manager, f, g,
                                              parallel);
  }

  template <class F, class G>
  inline void decodeInNgh(uintE vtx_id, F& f, G& g, bool parallel = true) {
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

  /* packing primitives */
  template <class P>
  inline size_t packOutNgh(uintE vtx_id, P& p, uint8_t* tmp, bool parallel = true, const flags fl=0) {
    return block_manager.pack_blocks(vtx_id, p, tmp, parallel, fl);
  }

  inline size_t calculateOutTemporarySpaceBytes() {
    if (block_manager.vtx_degree > 0) {
      size_t nblocks =  block_manager.vtx_num_blocks;
      if (nblocks > block_manager.kBlockAllocThreshold) {
        return (sizeof(uintE) * nblocks) + (block_manager.bytes_per_block * nblocks);
      }
    }
    return 0;
  }

  inline void clear_vertex() {
    block_manager.clear_vertex();
  }

  /* packing primitives */
  template <class P>
  inline size_t packInNghs(uintE vtx_id, P& p, uint8_t* tmp, bool parallel = true, const flags fl=0) {
    return packOutNgh(vtx_id, p, tmp, parallel, fl);
  }

  auto getOutIter() -> iter_type {
    return block_manager.get_iter();
  }

  auto getInIter() -> iter_type {
    return getOutIter();
  }

  size_t intersect(packed_symmetric_vertex<W, BM>& other) {
    auto it = getOutIter();
    auto other_it = other.getOutIter();
    return block_vertex_ops::intersect(it, other_it);
  }

  template <class S>
  size_t intersect(S& seq, packed_symmetric_vertex<W, BM>& other) {
    auto other_it = other.getOutIter();
    return block_vertex_ops::intersect_seq(seq, other_it);
  }
};

// Augments an ordinary (immutable) graph with the ability to filter/pack out
// edges.
template <template <class W> class vertex_type, class W>
struct packed_graph {
  size_t n; /* number of vertices */
  size_t m; /* number of edges; updated by decremental updates  */
  symmetric_graph<vertex, W>& GA;
  using vertex = vertex_type<W>;
  using E = typename vertex::edge_type;
  using weight_type = W;

  vtx_info* VI;
  uint8_t* blocks;

  size_t bs;
  size_t bs_in_bytes;
  size_t metadata_size;

  // Initializes the block memory for each vertex.
  void init_block_memory() {
    // 1. Calculate the #bytes corresponding to each vertex

    auto block_bytes_offs = pbbs::sequence<size_t>(n + 1);
    parallel_for(0, n, [&](size_t i) {
      uintE degree = GA.get_vertex(i).getOutDegree();
      block_bytes_offs[i] =
          bitsets::bytes_for_degree_and_bs(degree, bs, bs_in_bytes);
    });
    block_bytes_offs[n] = 0;

    size_t block_mem_to_alloc =
        pbbslib::scan_add_inplace(block_bytes_offs.slice());
    std::cout << "# total memory for packed_graph = " << block_mem_to_alloc << std::endl;

    auto blocks_seq = pbbs::delayed_seq<size_t>(n, [&] (size_t i) {
      uintE degree = GA.get_vertex(i).getOutDegree();
      if (degree == 0) { return static_cast<size_t>(0); }
      size_t nb = pbbs::num_blocks(degree, bs);
      return nb;
    });
    size_t total_blocks = pbbslib::reduce_add(blocks_seq);

    // allocate blocks
    blocks = pbbs::new_array_no_init<uint8_t>(block_mem_to_alloc);
    VI = pbbs::new_array_no_init<vtx_info>(n);
    std::cout << "# packed graph: block_size = " << bs << std::endl;
    std::cout << "# total blocks = " << total_blocks << " sizeof(vtx_info) = " << (sizeof(vtx_info)) << " total vtx_info bytes = " << (n*sizeof(vtx_info)) << std::endl;
    std::cout << "# total memory usage = " << (block_mem_to_alloc + (n*sizeof(vtx_info))) << " bytes " << std::endl;

    // initialize blocks and vtx_info
    parallel_for(
        0, n,
        [&](size_t v) {
          uintE degree = GA.get_vertex(v).getOutDegree();
          size_t block_byte_offset = block_bytes_offs[v];
          size_t vtx_bytes = block_bytes_offs[v + 1] - block_byte_offset;

          size_t num_blocks = pbbs::num_blocks(degree, bs);
          // set vertex_info for v
          VI[v] = vtx_info(degree, num_blocks, block_byte_offset);

          // initialize blocks corresponding to v's neighbors
          uint8_t* our_block_start = blocks + block_byte_offset;
          bitsets::bitset_init_blocks(our_block_start, degree, num_blocks, bs, bs_in_bytes,
                                      vtx_bytes);
        },
        1);
  }

  template <
      bool bool_enable = true,
      typename std::enable_if<
          std::is_same<vertex<W>, symmetric_vertex<W>>::value && bool_enable,
          int>::type = 0>
  void init_block_size_and_metadata() {
    using block_manager = sym_bitset_manager<vertex, W>;
    using metadata = typename block_manager::metadata;
    bs = block_manager::edges_per_block;
    bs_in_bytes = bitsets::bitset_bs_bytes(bs);
    metadata_size = sizeof(metadata);
  }

  template <bool bool_enable = true,
            typename std::enable_if<
                std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value &&
                    bool_enable,
                int>::type = 0>
  void init_block_size_and_metadata() {
    using block_manager = compressed_sym_bitset_manager<vertex, W>;
    using metadata = typename block_manager::metadata;
    bs = block_manager::edges_per_block;
    bs_in_bytes = bitsets::bitset_bs_bytes(bs);
    metadata_size = sizeof(metadata);
  }

  packed_graph(symmetric_graph<vertex, W>& GA) : n(GA.n), m(GA.m), GA(GA) {
    init_block_size_and_metadata();  // conditioned on vertex type
    init_block_memory();
  }

  inline size_t getOutDegree(uintE v) {
    return VI[v].vtx_degree;
  }

  inline size_t getInDegree(uintE v) {
    return getOutDegree();
  }

  // Symmetric: get_vertex
  template <
      bool bool_enable = true,
      typename std::enable_if<
          std::is_same<vertex<W>, symmetric_vertex<W>>::value && bool_enable,
          int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using block_manager = sym_bitset_manager<vertex, W>;
    auto vtx_data = GA.V[v];
    uintE original_degree = vtx_data.degree;
    uintE offset = vtx_data.offset;
#ifndef NVM
    auto sym_blocks = block_manager(v, blocks, original_degree, VI, GA.e0 + offset);
#else
    auto sym_blocks = block_manager(v, blocks, original_degree, VI, GA.e0 + offset, GA.e1 + offset);
#endif
    return packed_symmetric_vertex<W, block_manager>(std::move(sym_blocks));
  }

  // Compressed Symmetric: get_vertex
  template <bool bool_enable = true,
            typename std::enable_if<
                std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value &&
                    bool_enable,
                int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using block_manager = compressed_sym_bitset_manager<vertex, W>;
    auto vtx_data = GA.V[v];
    uintE original_degree = vtx_data.degree;
    size_t offset = vtx_data.offset;
#ifndef NVM
    auto sym_blocks = block_manager(v, blocks, original_degree, VI, GA.e0 + offset);
#else
    auto sym_blocks = block_manager(v, blocks, original_degree, VI, GA.e0 + offset, GA.e1 + offset);
#endif
    return packed_symmetric_vertex<W, block_manager>(std::move(sym_blocks));
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    par_for(0, n, 1, [&](size_t v) {
      auto vert_v = get_vertex(v);
      vert_v.mapOutNgh(v, f, parallel_inner_map);
    });
  }

  void del() {
    std::cout << "# deleting packed_graph" << std::endl;
    pbbs::free_array(VI);
    pbbs::free_array(blocks);
  }
};

// Used to infer template arguments
template <template <class W> class vertex, class W>
auto build_packed_graph(symmetric_graph<vertex, W>& GA) {
  return packed_graph<vertex, W>(GA);
}

template <template <class W> class vertex, class W, class P>
packed_graph<vertex, W> filter_graph(symmetric_graph<vertex, W>& G, P& pred_f) {
  // TODO: do allocations, but in a (medium) constant number of allocations.
  auto GA = packed_graph<vertex, W>(G);
  {
    parallel_for(0, G.n, [&] (size_t v) {
      auto vtx = GA.get_vertex(v);
      if (vtx.getOutDegree() > 0) {
        vtx.packOutNgh(v, pred_f, /* tmp = */ nullptr, /* parallel = */true, /* flags = */ compact_blocks);
      }
    }, 1);
  }
  auto degree_seq = pbbs::delayed_seq<size_t>(GA.n, [&] (size_t i) {
    return GA.getOutDegree(i);
  });
  auto new_m = pbbslib::reduce_add(degree_seq);
  GA.m = new_m;
  std::cout << "# Returning new packed graph, new m = " << new_m << std::endl;
  return GA;
}

template <template <class W> class vertex, class W, class P>
void filter_graph(packed_graph<vertex, W>& GA, P& pred_f) {
  // TODO: do allocations, but in a (medium) constant number of allocations.
  {
    parallel_for(0, GA.n, [&] (size_t v) {
      auto vtx = GA.get_vertex(v);
      if (vtx.getOutDegree() > 0) {
        vtx.packOutNgh(v, pred_f, /* tmp = */ nullptr, /* parallel = */ true, /* flags = */ compact_blocks);
      }
    }, 1);
  }
  auto degree_seq = pbbs::delayed_seq<size_t>(GA.n, [&] (size_t i) {
    return GA.getOutDegree(i);
  });
  auto new_m = pbbslib::reduce_add(degree_seq);
  GA.m = new_m;
  std::cout << "# Packing packed graph: new m = " << new_m << std::endl;
}

