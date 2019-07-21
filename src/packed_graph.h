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

//// A wrapper around an ordinary (immutable) graph.
//// * Extends vertices with a bit-packing structure that admits fast decremental
////   updates.
//// * constructor takes an ordinary graph, and generates a packed_graph.
//template <class vertex, class W>
//struct packed_graph {
//  size_t n; /* number of vertices */
//  size_t m; /* number of edges; updated by decremental updates  */
//
//  using D = typename vertex::block_decode;
//
//  vertex* V;
//  block_vertex<W, D>* BV;
//
//  // block_vertex holds:
//  // * bitset (blocks)
//  // * degree
//
//  packed_graph(graph<vertex>& GA) {
//
//  }
//
//  // Would be good if we can write the NVM code here s.t. we can have internal
//  // parallelism over out-edges. (i.e. we fetch the correct inner-vertex's block
//  // based on numa-node, so a get_block function defined on the vertex)
//  auto get_vertex(uintE v) {
//    return BV[v];
//  }
//
//  template <class F>
//  void map_edges(F f, bool parallel_inner_map=true) {
//    par_for(0, n, 1, [&] (size_t v) {
//      BV[v].mapOutNgh(v, f, parallel_inner_map);
//    });
//  }
//}

  // block manager knows about (1) block size and (2) num blocks.


static constexpr size_t kUncompressedBlockSize = 1024;

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
    return block_symmetric_vertex<W, block_manager>(std::move(sym_blocks));
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
    return block_symmetric_vertex<W, block_manager>(std::move(sym_blocks));
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map=true) {
    par_for(0, n, 1, [&] (size_t v) {
      auto vert_v = get_vertex(v);
      vert_v.mapOutNgh(v, f, parallel_inner_map);
    });
  }
};

template <template <class W> class vertex, class W>
struct packed_graph {
  size_t n; /* number of vertices */
  size_t m; /* number of edges; updated by decremental updates  */
  graph<vertex<W>>& GA;

  // block degree computable by differencing two starts.
  struct vtx_info = {
    uintE degree; // vertex's (packed) degree.
    uintE num_blocks; // number of blocks associated with v
    size_t offset; // pointer into the block structure
    uintE* edges; // pointer to the original edges (prevents one random read)
    vtx_info(uintE degree, uintE num_blocks, size_t offset, uintE* edges) :
      degree(degree), num_blocks(num_blocks), offsets(offset), edges(edges) {}
  };

  vtx_info* VI;
  uint8_t* blocks;

  template <bool bool_enable=true, typename
  std::enable_if<std::is_same<vertex<W>, symmetricVertex<W>>::value && bool_enable, int>::type = 0>
  packed_graph(graph<vertex<W>>& GA) : n(GA.n), m(GA.m), GA(GA) {
    constexpr size_t bs = kUncompressedBitsetBlockSize;
    auto block_offs = pbbs::sequence<size_t>(n+1);
    parallel_for(0, n, [&] (size_t i) {
      size_t degree = GA.get_vertex(i).getOutDegree();
      block_offs[i] = pbbs::num_blocks(degree, bs);
    });
    block_offs[n] = 0;
    size_t total_blocks = pbbslib::scan_add_inplace(block_offs.slice());
    size_t bs_in_bytes = bitset_bs_bytes(kUncompressedBitsetBlockSize);
    size_t block_mem_to_alloc = total_blocks*bs_in_bytes;

    cout << "total blocks to allocate = " << total_blocks
      << " bs_in_bytes = " << bs_in_bytes << endl;
    cout << "total memory = " << block_mem_to_alloc << endl;

    // allocate blocks
    blocks = pbbs::new_array_no_init<uint8_t>(block_mem_to_alloc);

    // initialize blocks and vtx_info
    parallel_for(0, n, [&] (size_t v) {
      uintE degree = GA.get_vertex(i).getOutDegree();
      uintE* edges = GA.get_vertex(i).getOutNeighbors();
      size_t offset = block_offs[v];
      size_t num_blocks = block_offs[v+1] - offset;
      VI[v] = vtx_info(degree, num_blocks, offset, edges);

      uint8_t* our_block_start = blocks + (offset*bs_in_bytes);
      bitset_init_all_blocks(our_block_start, num_blocks, bs, bs_in_bytes);
    }, 1);
  }

  template <bool bool_enable=true, typename
    std::enable_if<std::is_same<vertex<W>, symmetricVertex<W>>::value && bool_enable, int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    auto& v_info = VI[v];
    using block_manager = sym_bitset_manager<vertex, W>;
    uint8_t* v_blocks_start = blocks + v_info.offset;
    auto sym_blocks = block_manager(v, v_info.degree, v_info.num_blocks, v_blocks_start);
    return block_symmetric_vertex<W, block_manager>(std::move(sym_blocks));
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map=true) {
    par_for(0, n, 1, [&] (size_t v) {
      auto vert_v = get_vertex(v);
      vert_v.mapOutNgh(v, f, parallel_inner_map);
    });
  }
};






// Need two separate bit-packing implementations based on the underlying vertex
// type (one for uncompressed vertices, and one for compressed vertices).

// Implement filter_graph, etc internally on the packed_graph data structure.
// The hope is that this will be strictly faster than the previous
// implementations we used.

// Constructor method for a packed_graph takes an ordinary graph and converts it
// to a packed_graph. Note that graphs in the system are now immutable. The only
// way to modify a graph is by generating a packed_graph

// templatize on the vertex type <Vertex>. Supply this to the bit-packing
// structure, which uses type_traits to generate one of two implementations. We
// then hold a packed_structure_name<V> as a class member.

// Question now is how to implement the vertex interfaces using the packed
// structure.

// One option is to design a packed_vertex and packed_compressed_vertex
// structure.
//
// This requires a large amount of code rewriting/duplication, though.
//
//
//
//
// Another way is to change vertex and compressed_vertex to operate on a given
// block size. For vertex, we can configure the block size, and for
// compressed_vertex it is just the block_size of the underlying compression
// format.
//
// Now, the decoding methods are structured to be specified
// (1) which blocks to decode using the block-structure supplied by the graph class
// (2) the specification of what to do for each block which is of the form of
// checking the i'th bit of the block to see whether the i'th neighbor needs
// decoding.
//
// For ordinary (non-packed) graphs, the decoding structure simply generates all
// of the blocks for the vertex, and returns true for all i-th bit queries.
//
// After implementing we can measure the overhead this extra code/inlining adds
// to the regular (unpacked) structure. Hopefully it is minimal.

// Issue is for intersection. Need to write a (parallel) block intersection.
// Note that the same algorithm could be used for both compressed and
// uncompressed since both work on blocks now.

// Is there a higher-level way of synthesizing a vertex implementation for the
// blocked_format now? i.e. a blocked_vertex.h class

// A blocked_vertex takes an implementation to decode blocks (trivial for
// uncompressed, the decompression method for the compressed case).




