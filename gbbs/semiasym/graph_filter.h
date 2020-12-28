// This code is part of the project "Sage: Parallel Semi-Asymmetric Graph
// Algorithms for NVRAMs" presented at VLDB 2020.
// Copyright (c) 2020 Laxman Dhulipala and the GBBS Team
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

// This file contains an implementation of a bit-packed graph, which is a
// decremental graph structure that supports edge-deletions using a bit-packed
// vector.

#pragma once

#include <functional>
#include <tuple>
#include <type_traits>

#include "bitset_managers.h"
#include "gbbs/graph.h"
#include "gbbs/flags.h"

namespace gbbs {
namespace sage {

/*
 * BitsetNeighbors supplies:
 *  - num_blocks TODO(update here?)
 *  - decode_block, decode_block_cond
 */
template <class W, /* weight type */
          class BitsetNeighbors /* block manager type */>
struct packed_symmetric_vertex {
  using iter_type = typename BitsetNeighbors::iter;
  BitsetNeighbors neighbors;  // copy; not a reference.

  packed_symmetric_vertex(BitsetNeighbors&& neighbors)
      : neighbors(std::move(neighbors)) {}

  BitsetNeighbors in_neighbors() { return neighbors; }
  BitsetNeighbors out_neighbors() { return neighbors; }

  __attribute__((always_inline)) inline uintE out_degree() {
    return neighbors.get_degree();
  }
  __attribute__((always_inline)) inline uintE in_degree() {
    return out_degree();
  }
};


// Initializes the block memory for each vertex.
// Returns a pair of (vtx_info*, blocks*)
template <class Graph>
std::pair<vtx_info*, uint8_t*> init_block_memory(Graph& GA, size_t bs, size_t bs_in_bytes, gbbs::flags fl= 0) {
  size_t n = GA.n;

  // 1. Calculate the #bytes corresponding to each vertex
  auto block_bytes_offs = pbbs::sequence<size_t>(n + 1);
  parallel_for(0, n, [&](size_t i) {
    uintE degree = (fl & in_edges) ? GA.get_vertex(i).in_degree() : GA.get_vertex(i).out_degree();
    block_bytes_offs[i] =
        bitsets::bytes_for_degree_and_bs(degree, bs, bs_in_bytes);
  });
  block_bytes_offs[n] = 0;

  size_t block_mem_to_alloc =
      pbbslib::scan_add_inplace(block_bytes_offs.slice());
  std::cout << "# total memory for symmetric_packed_graph = " << block_mem_to_alloc << std::endl;

  auto blocks_seq = pbbs::delayed_seq<size_t>(n, [&] (size_t i) {
    uintE degree = (fl & in_edges) ? GA.get_vertex(i).in_degree() : GA.get_vertex(i).out_degree();
    if (degree == 0) { return static_cast<size_t>(0); }
    size_t nb = pbbs::num_blocks(degree, bs);
    return nb;
  });
  size_t total_blocks = pbbslib::reduce_add(blocks_seq);

  // allocate blocks
  auto blocks = pbbs::new_array_no_init<uint8_t>(block_mem_to_alloc);
  auto VI = pbbs::new_array_no_init<vtx_info>(n);
  std::cout << "# packed graph: block_size = " << bs << std::endl;
  std::cout << "# total blocks = " << total_blocks << " sizeof(vtx_info) = " << (sizeof(vtx_info)) << " total vtx_info bytes = " << (n*sizeof(vtx_info)) << std::endl;
  std::cout << "# total memory usage = " << (block_mem_to_alloc + (n*sizeof(vtx_info))) << " bytes; in_edges = " << (fl & in_edges) << std::endl;

  // initialize blocks and vtx_info
  parallel_for(
      0, n,
      [&](size_t v) {
        uintE degree = (fl & in_edges) ? GA.get_vertex(v).in_degree() : GA.get_vertex(v).out_degree();
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

  return {VI, blocks};
}


// Augments an ordinary (immutable) graph with the ability to filter/pack out
// edges.
template <template <class W> class vertex_type, class W>
struct symmetric_packed_graph {
  size_t n;  // number of vertices
  size_t m;  // number of edges; updated by decremental updates

  symmetric_graph<vertex_type, W>& GA;  // the underlying (immutable) graph
  using vertex = vertex_type<W>;
  using weight_type = W;
  using E = typename vertex::edge_type;

  vtx_info* VI;  // metadata for each vertex into the blocks structure
  uint8_t* blocks;  // stores the block information for each vertex's neighbors

  size_t bs;  // number of vertices in each block
  size_t bs_in_bytes;
  size_t metadata_size;

  template <
      bool bool_enable = true,
      typename std::enable_if<
          std::is_same<vertex, symmetric_vertex<W>>::value && bool_enable,
          int>::type = 0>
  void init_block_size_and_metadata() {
    using neighborhood_type = uncompressed_bitset_neighbors<vertex_type, W>;
    using metadata = typename neighborhood_type::metadata;
    bs = neighborhood_type::edges_per_block;
    bs_in_bytes = bitsets::bitset_bs_bytes(bs);
    metadata_size = sizeof(metadata);
  }

  template <bool bool_enable = true,
            typename std::enable_if<
                std::is_same<vertex, csv_bytepd_amortized<W>>::value &&
                    bool_enable,
                int>::type = 0>
  void init_block_size_and_metadata() {
    using neighborhood_type = compressed_bitset_neighbors<vertex_type, W>;
    using metadata = typename neighborhood_type::metadata;
    bs = neighborhood_type::edges_per_block;
    bs_in_bytes = bitsets::bitset_bs_bytes(bs);
    metadata_size = sizeof(metadata);
  }

  symmetric_packed_graph(symmetric_graph<vertex_type, W>& GA) : n(GA.n), m(GA.m), GA(GA) {
    init_block_size_and_metadata();  // conditioned on vertex type; initializes bs, bs_in_bytes, metadata_size
    std::tie(VI, blocks) = init_block_memory(GA, bs, bs_in_bytes);  // initializes VI and blocks based on bs and bs_in_bytes
  }

  inline size_t out_degree(uintE v) {
    return VI[v].vtx_degree;
  }

  inline size_t in_degree(uintE v) {
    return out_degree();
  }

  // Symmetric: get_vertex
  template <
      bool bool_enable = true,
      typename std::enable_if<
          std::is_same<vertex, symmetric_vertex<W>>::value && bool_enable,
          int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using neighborhood_type = uncompressed_bitset_neighbors<vertex_type, W>;
    auto vtx_data = GA.v_data[v];
    uintE original_degree = vtx_data.degree;
    uintE offset = vtx_data.offset;
#ifndef SAGE
    auto sym_blocks = neighborhood_type(v, blocks, original_degree, VI, GA.e0 + offset);
#else
    auto sym_blocks = neighborhood_type(v, blocks, original_degree, VI, GA.e0 + offset, GA.e1 + offset);
#endif
    return packed_symmetric_vertex<W, neighborhood_type>(std::move(sym_blocks));
  }

  // Compressed Symmetric: get_vertex
  template <bool bool_enable = true,
            typename std::enable_if<
                std::is_same<vertex, csv_bytepd_amortized<W>>::value &&
                    bool_enable,
                int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using neighborhood_type = compressed_bitset_neighbors<vertex_type, W>;
    auto vtx_data = GA.v_data[v];
    uintE original_degree = vtx_data.degree;
    size_t offset = vtx_data.offset;
#ifndef SAGE
    auto sym_blocks = neighborhood_type(v, blocks, original_degree, VI, GA.e0 + offset);
#else
    auto sym_blocks = neighborhood_type(v, blocks, original_degree, VI, GA.e0 + offset, GA.e1 + offset);
#endif
    return packed_symmetric_vertex<W, neighborhood_type>(std::move(sym_blocks));
  }

  template <class P>
  uintE packNeighbors(uintE id, P& p, uint8_t* tmp) {
    auto vtx = get_vertex(id);
    return vtx.out_neighbors().pack(p, tmp, false);
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    par_for(0, n, 1, [&](size_t v) {
      auto vert_v = get_vertex(v);
      vert_v.out_neighbors().map(f, parallel_inner_map);
    });
  }

  void del() {
    std::cout << "# deleting packed_graph" << std::endl;
    pbbs::free_array(VI);
    pbbs::free_array(blocks);
  }
};

// Used to infer template arguments
template <template <class W> class vertex_type, class W>
auto build_symmetric_packed_graph(symmetric_graph<vertex_type, W>& GA) {
  return symmetric_packed_graph<vertex_type, W>(GA);
}

template <template <class W> class vertex_type, class W, class P>
symmetric_packed_graph<vertex_type, W> filter_graph(symmetric_graph<vertex_type, W>& G, P& pred_f) {
  // TODO: do allocations, but in a (medium) constant number of allocations.
  auto GA = symmetric_packed_graph<vertex_type, W>(G);
  {
    parallel_for(0, G.n, [&] (size_t v) {
      auto vtx = GA.get_vertex(v);
      if (vtx.out_degree() > 0) {
        vtx.out_neighbors().pack(pred_f, /* tmp = */ nullptr, /* parallel = */true, /* flags = */ compact_blocks);
      }
    }, 1);
  }
  auto degree_seq = pbbs::delayed_seq<size_t>(GA.n, [&] (size_t i) {
    return GA.out_degree(i);
  });
  auto new_m = pbbslib::reduce_add(degree_seq);
  GA.m = new_m;
  std::cout << "# Returning new packed graph, new m = " << new_m << std::endl;
  return GA;
}

template <template <class W> class vertex_type, class W, class P>
void filter_graph(symmetric_packed_graph<vertex_type, W>& GA, P& pred_f) {
  // TODO: do allocations, but in a (medium) constant number of allocations.
  {
    parallel_for(0, GA.n, [&] (size_t v) {
      auto vtx = GA.get_vertex(v);
      if (vtx.out_degree() > 0) {
        vtx.out_neighbors().pack(pred_f, /* tmp = */ nullptr, /* parallel = */ true, /* flags = */ compact_blocks);
      }
    }, 1);
  }
  auto degree_seq = pbbs::delayed_seq<size_t>(GA.n, [&] (size_t i) {
    return GA.out_degree(i);
  });
  auto new_m = pbbslib::reduce_add(degree_seq);
  GA.m = new_m;
  std::cout << "# Packing packed graph: new m = " << new_m << std::endl;
}



/*
 * BitsetNeighbors supplies:
 *  - num_blocks TODO(update here?)
 *  - decode_block, decode_block_cond
 */
template <class W, /* weight type */
          class BitsetNeighbors /* block manager type */>
struct packed_asymmetric_vertex {
  using iter_type = typename BitsetNeighbors::iter;
  BitsetNeighbors in_nghs;  // copy; not a reference.
  BitsetNeighbors out_nghs;  // copy; not a reference.

  packed_asymmetric_vertex(BitsetNeighbors&& in_nghs, BitsetNeighbors&& out_nghs)
      : in_nghs(std::move(in_nghs)), out_nghs(std::move(out_nghs)) {}

  BitsetNeighbors in_neighbors() { return in_nghs; }
  BitsetNeighbors out_neighbors() { return out_nghs; }

  __attribute__((always_inline)) inline uintE out_degree() {
    return out_nghs.get_degree();
  }
  __attribute__((always_inline)) inline uintE in_degree() {
    return in_nghs.get_degree();
  }
};


// Augments an ordinary (immutable) graph with the ability to filter/pack out
// edges.
template <template <class W> class vertex_type, class W>
struct asymmetric_packed_graph {
  size_t n;  // number of vertices
  size_t m;  // number of (directed) edges; updated by decremental updates

  asymmetric_graph<vertex_type, W>& GA;  // the underlying (immutable) graph
  using vertex = vertex_type<W>;
  using weight_type = W;
  using E = typename vertex::edge_type;

  vtx_info* in_VI;  // metadata for each vertex into the in_blocks structure
  uint8_t* in_blocks;  // stores the block information for each vertex's in_neighbors

  vtx_info* out_VI;  // metadata for each vertex into the out_blocks structure
  uint8_t* out_blocks;  // stores the block information for each vertex's out_neighbors

  size_t bs;  // number of vertices in each block
  size_t bs_in_bytes;
  size_t metadata_size;

  template <
      bool bool_enable = true,
      typename std::enable_if<
          std::is_same<vertex, asymmetric_vertex<W>>::value && bool_enable,
          int>::type = 0>
  void init_block_size_and_metadata() {
    using neighborhood_type = uncompressed_bitset_neighbors<vertex_type, W>;
    using metadata = typename neighborhood_type::metadata;
    bs = neighborhood_type::edges_per_block;
    bs_in_bytes = bitsets::bitset_bs_bytes(bs);
    metadata_size = sizeof(metadata);
  }

  template <bool bool_enable = true,
            typename std::enable_if<
                std::is_same<vertex, cav_bytepd_amortized<W>>::value &&
                    bool_enable,
                int>::type = 0>
  void init_block_size_and_metadata() {
    using neighborhood_type = compressed_bitset_neighbors<vertex_type, W>;
    using metadata = typename neighborhood_type::metadata;
    bs = neighborhood_type::edges_per_block;
    bs_in_bytes = bitsets::bitset_bs_bytes(bs);
    metadata_size = sizeof(metadata);
  }

  asymmetric_packed_graph(asymmetric_graph<vertex_type, W>& GA) : n(GA.n), m(GA.m), GA(GA) {
    init_block_size_and_metadata();  // conditioned on vertex type; initializes bs, bs_in_bytes, metadata_size
    std::tie(in_VI, in_blocks) = init_block_memory(GA, bs, bs_in_bytes, in_edges);
    std::tie(out_VI, out_blocks) = init_block_memory(GA, bs, bs_in_bytes);
  }

  inline size_t out_degree(uintE v) {
    return out_VI[v].vtx_degree;
  }

  inline size_t in_degree(uintE v) {
    return in_VI[v].vtx_degree;
  }

  // Symmetric: get_vertex
  template <
      bool bool_enable = true,
      typename std::enable_if<
          std::is_same<vertex, asymmetric_vertex<W>>::value && bool_enable,
          int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using neighborhood_type = uncompressed_bitset_neighbors<vertex_type, W>;

    auto in_vtx_data = GA.v_in_data[v];
    uintE in_original_degree = in_vtx_data.degree;
    uintE in_offset = in_vtx_data.offset;

    auto out_vtx_data = GA.v_out_data[v];
    uintE out_original_degree = out_vtx_data.degree;
    uintE out_offset = out_vtx_data.offset;

#ifndef SAGE
    auto in_neighborhood = neighborhood_type(v, in_blocks, in_original_degree, in_VI, GA.in_edges_0 + in_offset);
    auto out_neighborhood = neighborhood_type(v, out_blocks, out_original_degree, out_VI, GA.out_edges_0 + out_offset);
#else
    auto in_neighborhood = neighborhood_type(v, in_blocks, in_original_degree, in_VI, GA.in_edges_0 + in_offset, GA.in_eddes_1 + in_offset);
    auto out_neighborhood = neighborhood_type(v, out_blocks, out_original_degree, out_VI, GA.out_edges_0 + out_offset, GA.out_edges_1 + out_offset);
#endif
    return packed_asymmetric_vertex<W, neighborhood_type>(std::move(in_neighborhood), std::move(out_neighborhood));
  }

  // Compressed Symmetric: get_vertex
  template <bool bool_enable = true,
            typename std::enable_if<
                std::is_same<vertex, cav_bytepd_amortized<W>>::value &&
                    bool_enable,
                int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using neighborhood_type = compressed_bitset_neighbors<vertex_type, W>;

    auto in_vtx_data = GA.v_in_data[v];
    uintE in_original_degree = in_vtx_data.degree;
    uintE in_offset = in_vtx_data.offset;

    auto out_vtx_data = GA.v_out_data[v];
    uintE out_original_degree = out_vtx_data.degree;
    uintE out_offset = out_vtx_data.offset;

#ifndef SAGE
    auto in_neighborhood = neighborhood_type(v, in_blocks, in_original_degree, in_VI, GA.in_edges_0 + in_offset);
    auto out_neighborhood = neighborhood_type(v, out_blocks, out_original_degree, out_VI, GA.out_edges_0 + out_offset);
#else
    auto in_neighborhood = neighborhood_type(v, in_blocks, in_original_degree, in_VI, GA.in_edges_0 + in_offset, GA.in_edges_1 + in_offset);
    auto out_neighborhood = neighborhood_type(v, out_blocks, out_original_degree, out_VI, GA.out_edges_0 + out_offset, GA.out_edges_1 + out_offset);
#endif
    return packed_asymmetric_vertex<W, neighborhood_type>(std::move(in_neighborhood), std::move(out_neighborhood));
  }

  template <class P>
  uintE packNeighbors(uintE id, P& p, uint8_t* tmp) {
    auto vtx = get_vertex(id);
    return vtx.out_neighbors().pack(p, tmp, false);
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    par_for(0, n, 1, [&](size_t v) {
      auto vert_v = get_vertex(v);
      vert_v.out_neighbors().map(f, parallel_inner_map);
    });
  }

  void del() {
    std::cout << "# deleting packed_graph" << std::endl;
    pbbslib::free_arrays(out_VI, in_VI, out_blocks, in_blocks);
  }
};

// Used to infer template arguments
template <template <class W> class vertex_type, class W>
auto build_asymmetric_packed_graph(asymmetric_graph<vertex_type, W>& GA) {
  return asymmetric_packed_graph<vertex_type, W>(GA);
}


}  // namespace sage
}  // namespace gbbs
