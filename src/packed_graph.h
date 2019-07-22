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

#include "graph.h"
#include "block_vertex.h"
#include "bitset_managers.h"

// TODO: define concept for packable graph?
//#ifdef CONCEPTS
//template<typename W>
//concept bool PackableGraph =
//  requires(T t, size_t u) {
//  typename T::value_type;
//  { t.size() } -> size_t;
//  { t.slice() };
//  { t[u] };
//};
//
//template<typename T>
//concept bool Range =
//  Seq<T> && requires(T t, size_t u) {
//  { t[u] } -> typename T::value_type&;
//  typename T::iterator;
//};
//#define SEQ Seq
//#define RANGE Range
//#else
//#define SEQ typename
//#define RANGE typename
//#endif

/*
 * block_manager supplies:
 *  - num_blocks
 *  - decode_block, decode_block_cond
 */
template <class W, /* weight type */
          class BM /* block manager type */>
struct packed_symmetric_vertex {

  BM block_manager; // copy; not a reference.

  packed_symmetric_vertex(BM&& block_manager) :
    block_manager(std::move(block_manager)) {}

  __attribute__((always_inline)) inline uintE getOutDegree() {
    return block_manager.get_degree();
  }
  __attribute__((always_inline)) inline uintE getInDegree() {
    return getOutDegree();
  }

  __attribute__((always_inline)) inline uintE getNumInBlocks() { return block_manager.num_blocks(); }
  __attribute__((always_inline)) inline uintE getNumOutBlocks() { return getNumInBlocks(); }
  __attribute__((always_inline)) inline uintE in_block_degree(uintE block_num) {
    return block_manager.block_degree(block_num);
  }
  __attribute__((always_inline)) inline uintE out_block_degree(uintE block_num) {
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

  /* packing primitives */
  template <class P, class E>
  inline size_t packOutNgh(uintE vtx_id, P& p, E* tmp, bool parallel = true) {
    block_manager.pack_blocks(vtx_id, p, tmp, parallel);
  }

  /* packing primitives */
  template <class P, class E>
  inline size_t packInNghs(uintE vtx_id, P& p, E* tmp, bool parallel = true) {
    return packOutNgh(vtx_id, p, tmp, parallel);
  }
};

// Augments an ordinary (immutable) graph with the ability to filter/pack out
// edges.
template <template <class W> class vertex, class W>
struct packed_graph {
  size_t n; /* number of vertices */
  size_t m; /* number of edges; updated by decremental updates  */
  graph<vertex, W>& GA;
  using E = typename vertex<W>::E;
  using weight_type = W;

  vtx_info<E>* VI;
  uint8_t* blocks;

  size_t bs;
  size_t bs_in_bytes;
  size_t metadata_size;

  // Initializes the block memory for each vertex.
  void init_block_memory() {
    // 1. Calculate the #bytes corresponding to each vertex
    auto block_bytes_offs = pbbs::sequence<size_t>(n+1);
    parallel_for(0, n, [&] (size_t i) {
      uintE degree = GA.get_vertex(i).getOutDegree();
      block_bytes_offs[i] = bitsets::bytes_for_degree_and_bs(degree, bs, bs_in_bytes);
    });
    block_bytes_offs[n] = 0;

    size_t block_mem_to_alloc = pbbslib::scan_add_inplace(block_bytes_offs.slice());

    cout << "total memory for packed_graph = " << block_mem_to_alloc << endl;

    // allocate blocks
    blocks = pbbs::new_array_no_init<uint8_t>(block_mem_to_alloc);

    VI = pbbs::new_array_no_init<vtx_info<E>>(n);

    // initialize blocks and vtx_info
    parallel_for(0, n, [&] (size_t v) {
      uintE degree = GA.get_vertex(v).getOutDegree();
      E* edges = GA.get_vertex(v).getOutNeighbors();
      size_t block_byte_offset = block_bytes_offs[v];
      size_t vtx_bytes = block_bytes_offs[v+1] - block_byte_offset;

      size_t num_blocks = pbbs::num_blocks(degree, bs);
      // set vertex_info for v
      VI[v] = vtx_info<E>(degree, num_blocks, block_byte_offset, edges);

      // initialize blocks corresponding to v's neighbors
      uint8_t* our_block_start = blocks + block_byte_offset;
      bitsets::bitset_init_blocks(our_block_start, degree, num_blocks, bs, vtx_bytes);
    }, 1);
  }

  template <bool bool_enable=true, typename
  std::enable_if<std::is_same<vertex<W>, symmetricVertex<W>>::value && bool_enable, int>::type = 0>
  void init_block_size_and_metadata() {
    using block_manager = sym_bitset_manager<vertex, W>;
    using metadata = typename block_manager::metadata;
    bs = block_manager::edges_per_block;
    bs_in_bytes = bitsets::bitset_bs_bytes(bs);
    metadata_size = sizeof(metadata);
  }

  template <bool bool_enable=true, typename
  std::enable_if<std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value && bool_enable, int>::type = 0>
  void init_block_size_and_metadata() {
    using block_manager = compressed_sym_bitset_manager<vertex, W>;
    using metadata = typename block_manager::metadata;
    bs = block_manager::edges_per_block;
    bs_in_bytes = bitsets::bitset_bs_bytes(bs);
    metadata_size = sizeof(metadata);
  }

  packed_graph(graph<vertex, W>& GA) : n(GA.n), m(GA.m), GA(GA) {
    init_block_size_and_metadata(); // conditioned on vertex type
    init_block_memory();
  }

  // Symmetric: get_vertex
  template <bool bool_enable=true, typename
    std::enable_if<std::is_same<vertex<W>, symmetricVertex<W>>::value && bool_enable, int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using block_manager = sym_bitset_manager<vertex, W>;
    auto sym_blocks = block_manager(v, blocks, VI);
    return packed_symmetric_vertex<W, block_manager>(std::move(sym_blocks));
  }

  // Compressed Symmetric: get_vertex
  template <bool bool_enable=true, typename
    std::enable_if<std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value && bool_enable, int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using block_manager = compressed_sym_bitset_manager<vertex, W>;
    auto sym_blocks = block_manager(v, blocks, VI);
    return packed_symmetric_vertex<W, block_manager>(std::move(sym_blocks));
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map=true) {
    par_for(0, n, 1, [&] (size_t v) {
      auto vert_v = get_vertex(v);
      vert_v.mapOutNgh(v, f, parallel_inner_map);
    });
  }

  void del() {
  }
};
