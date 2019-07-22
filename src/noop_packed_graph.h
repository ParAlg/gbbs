#pragma once

#pragma once

#include <functional>
#include <tuple>
#include <type_traits>

#include "graph.h"
#include "block_vertex.h"
#include "noop_managers.h"

static constexpr size_t kUncompressedBlockSize = 1024;

/*
 * block_manager supplies:
 *  - num_blocks
 *  - decode_block, decode_block_cond
 */
template <class W, /* weight type */
          class BM /* block manager type */>
struct noop_symmetric_vertex {

  BM block_manager; // copy; not a reference.

  noop_symmetric_vertex(BM&& block_manager) :
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


// An allocation-free wrapper around an ordinary (immutable) graph.
// * The bit-packing structure is a no-op (returns all blocks in the
//   original graph for a vertex, and does not support packing).
// This is mainly for measuring the overhead due to the new interfaces (which
// seems to be pretty minimal---around 10% compared to using the immutable
// interfaces)
template <template <class W> class vertex, class W>
struct noop_packed_graph {
  size_t n; /* number of vertices */
  size_t m; /* number of edges  */
  graph<vertex<W>>& GA;

  using E = typename vertex<W>::E;
  using WW = W;

  noop_packed_graph(graph<vertex<W>>& GA) : n(GA.n), m(GA.m), GA(GA) { }

  template <bool bool_enable=true, typename
    std::enable_if<std::is_same<vertex<W>, symmetricVertex<W>>::value && bool_enable, int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using block_manager = sym_noop_manager<vertex, W>;
#ifndef NVM
    auto sym_blocks = block_manager(GA.V[v], v, kUncompressedBlockSize);
#else
    auto sym_blocks = block_manager(GA.V0[v], GA.V1[v], v, kUncompressedBlockSize);
#endif
    return noop_symmetric_vertex<W, block_manager>(std::move(sym_blocks));
  }

  template <bool bool_enable=true, typename
    std::enable_if<std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value && bool_enable, int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using block_manager = compressed_sym_noop_manager<vertex, W>;
#ifndef NVM
    auto sym_blocks = block_manager(GA.V[v], v);
#else
    auto sym_blocks = block_manager(GA.V0[v], GA.V1[v], v);
#endif
    return noop_symmetric_vertex<W, block_manager>(std::move(sym_blocks));
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
