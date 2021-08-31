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
#include "gbbs/flags.h"
#include "gbbs/graph.h"

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
template <class Graph, class Pred>
std::tuple<parlay::sequence<vtx_info>, parlay::sequence<uint8_t>, size_t>
init_block_memory(Graph& GA, size_t bs, size_t bs_in_bytes, Pred& vtx_pred,
                  gbbs::flags fl = 0) {
  timer ibm;
  ibm.start();
  size_t n = GA.n;

  // Calculate the #bytes corresponding to each vertex.
  auto block_bytes_offs = sequence<size_t>(n + 1);
  parallel_for(0, n, [&](size_t i) {
    uintE degree = 0;
    if (vtx_pred(i)) {
      degree = (fl & in_edges) ? GA.get_vertex(i).in_degree()
                               : GA.get_vertex(i).out_degree();
    }
    block_bytes_offs[i] =
        bitsets::bytes_for_degree_and_bs(degree, bs, bs_in_bytes);
  });
  block_bytes_offs[n] = 0;

  size_t block_mem_to_alloc =
      parlay::scan_inplace(make_slice(block_bytes_offs));
  std::cout << "# total memory for block memory = " << block_mem_to_alloc
            << std::endl;

  //  auto blocks_seq = parlay::delayed_seq<size_t>(n, [&] (size_t i) {
  //    uintE degree = (fl & in_edges) ? GA.get_vertex(i).in_degree() :
  //    GA.get_vertex(i).out_degree();
  //    if (degree == 0) { return static_cast<size_t>(0); }
  //    size_t nb = parlay::num_blocks(degree, bs);
  //    return nb;
  //  });
  //  size_t total_blocks = parlay::reduce(blocks_seq);
  //  std::cout << "# total blocks = " << total_blocks << std::endl;

  // allocate blocks
  auto blocks = parlay::sequence<uint8_t>::uninitialized(block_mem_to_alloc);
  auto VI = parlay::sequence<vtx_info>::uninitialized(n);
  std::cout << "# packed graph: block_size = " << bs << std::endl;
  std::cout << "# sizeof(vtx_info) = " << (sizeof(vtx_info))
            << " total vtx_info bytes = " << (n * sizeof(vtx_info))
            << std::endl;
  std::cout << "# total memory usage = "
            << (block_mem_to_alloc + (n * sizeof(vtx_info)))
            << " bytes; in_edges = " << (fl & in_edges) << std::endl;

  // initialize blocks and vtx_info
  parallel_for(0, n, [&](size_t v) {
    uintE degree = 0;
    if (vtx_pred(v)) {
      degree = (fl & in_edges) ? GA.get_vertex(v).in_degree()
                               : GA.get_vertex(v).out_degree();
    }
    size_t block_byte_offset = block_bytes_offs[v];
    size_t vtx_bytes = block_bytes_offs[v + 1] - block_byte_offset;

    size_t num_blocks = parlay::num_blocks(degree, bs);
    // set vertex_info for v
    VI[v] = vtx_info(degree, num_blocks, block_byte_offset);

    // initialize blocks corresponding to v's neighbors
    uint8_t* our_block_start = blocks.begin() + block_byte_offset;
    bitsets::bitset_init_blocks(our_block_start, degree, num_blocks, bs,
                                bs_in_bytes, vtx_bytes);
  });

  ibm.stop();
  ibm.next("init block memory time");
  return {std::move(VI), std::move(blocks), block_mem_to_alloc};
}

template <class Graph>
std::tuple<parlay::sequence<vtx_info>, parlay::sequence<uint8_t>, size_t>
init_block_memory(Graph& GA, size_t bs, size_t bs_in_bytes,
                  gbbs::flags fl = 0) {
  auto pred_true = [&](const uintE& v) { return true; };
  return init_block_memory(GA, bs, bs_in_bytes, pred_true, fl);
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

  parlay::sequence<vtx_info>
      VI;  // metadata for each vertex into the blocks structure
  parlay::sequence<uint8_t>
      blocks;  // stores the block information for each vertex's neighbors
  size_t blocks_mem;

  size_t bs;  // number of vertices in each block
  size_t bs_in_bytes;
  size_t metadata_size;

  size_t block_mem_to_alloc;

  template <bool bool_enable = true,
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

  template <
      bool bool_enable = true,
      typename std::enable_if<
          std::is_same<vertex, csv_bytepd_amortized<W>>::value && bool_enable,
          int>::type = 0>
  void init_block_size_and_metadata() {
    using neighborhood_type = compressed_bitset_neighbors<vertex_type, W>;
    using metadata = typename neighborhood_type::metadata;
    bs = neighborhood_type::edges_per_block;
    bs_in_bytes = bitsets::bitset_bs_bytes(bs);
    metadata_size = sizeof(metadata);
  }

  symmetric_packed_graph(symmetric_graph<vertex_type, W>& GA)
      : n(GA.n), m(GA.m), GA(GA) {
    init_block_size_and_metadata();  // conditioned on vertex type; initializes
                                     // bs, bs_in_bytes, metadata_size
    std::tie(VI, blocks, blocks_mem) = init_block_memory(
        GA, bs,
        bs_in_bytes);  // initializes VI and blocks based on bs and bs_in_bytes
  }

  inline size_t out_degree(uintE v) { return VI[v].vtx_degree; }

  inline size_t in_degree(uintE v) { return out_degree(); }

  // Symmetric: get_vertex
  template <bool bool_enable = true,
            typename std::enable_if<
                std::is_same<vertex, symmetric_vertex<W>>::value && bool_enable,
                int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using neighborhood_type = uncompressed_bitset_neighbors<vertex_type, W>;
    auto vtx_data = GA.v_data[v];
    uintE original_degree = vtx_data.degree;
    uintT offset = vtx_data.offset;
    auto sym_blocks = neighborhood_type(v, blocks.begin(), original_degree,
                                        VI.begin(), GA.e0 + offset);
    return packed_symmetric_vertex<W, neighborhood_type>(std::move(sym_blocks));
  }

  // Compressed Symmetric: get_vertex
  template <
      bool bool_enable = true,
      typename std::enable_if<
          std::is_same<vertex, csv_bytepd_amortized<W>>::value && bool_enable,
          int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using neighborhood_type = compressed_bitset_neighbors<vertex_type, W>;
    auto vtx_data = GA.v_data[v];
    uintE original_degree = vtx_data.degree;
    size_t offset = vtx_data.offset;
    auto sym_blocks = neighborhood_type(v, blocks.begin(), original_degree,
                                        VI.begin(), GA.e0 + offset);
    return packed_symmetric_vertex<W, neighborhood_type>(std::move(sym_blocks));
  }

  template <class P>
  uintE packNeighbors(uintE id, P& p, uint8_t* tmp) {
    auto vtx = get_vertex(id);
    return vtx.out_neighbors().pack(p, tmp, false);
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    parallel_for(0, n,
                 [&](size_t v) {
                   auto vert_v = get_vertex(v);
                   vert_v.out_neighbors().map(f, parallel_inner_map);
                 },
                 1);
  }
};

// Used to infer template arguments
template <template <class W> class vertex_type, class W>
auto build_symmetric_packed_graph(symmetric_graph<vertex_type, W>& GA) {
  return symmetric_packed_graph<vertex_type, W>(GA);
}

template <template <class W> class vertex_type, class W, class P>
symmetric_packed_graph<vertex_type, W> filter_graph(
    symmetric_graph<vertex_type, W>& G, P& pred_f) {
  // TODO: do allocations, but in a (medium) constant number of allocations.
  auto GA = symmetric_packed_graph<vertex_type, W>(G);
  {
    parallel_for(0, G.n,
                 [&](size_t v) {
                   auto vtx = GA.get_vertex(v);
                   if (vtx.out_degree() > 0) {
                     vtx.out_neighbors().pack(pred_f, /* tmp = */ nullptr,
                                              /* parallel = */ true,
                                              /* flags = */ compact_blocks);
                   }
                 },
                 1);
  }
  auto degree_seq = parlay::delayed_seq<size_t>(
      GA.n, [&](size_t i) { return GA.out_degree(i); });
  auto new_m = parlay::reduce(degree_seq);
  GA.m = new_m;
  std::cout << "# Returning new packed graph, new m = " << new_m << std::endl;
  return GA;
}

template <template <class W> class vertex_type, class W, class P>
void filter_graph(symmetric_packed_graph<vertex_type, W>& GA, P& pred_f) {
  // TODO: do allocations, but in a (medium) constant number of allocations.
  {
    parallel_for(0, GA.n,
                 [&](size_t v) {
                   auto vtx = GA.get_vertex(v);
                   if (vtx.out_degree() > 0) {
                     vtx.out_neighbors().pack(pred_f, /* tmp = */ nullptr,
                                              /* parallel = */ true,
                                              /* flags = */ compact_blocks);
                   }
                 },
                 1);
  }
  auto degree_seq = parlay::delayed_seq<size_t>(
      GA.n, [&](size_t i) { return GA.out_degree(i); });
  auto new_m = parlay::reduce(degree_seq);
  GA.m = new_m;
  std::cout << "# Packing packed graph: new m = " << new_m << std::endl;
}

/* ========================== Asymmetric Utils ========================= */

/*
 * BitsetNeighbors supplies:
 *  - num_blocks TODO(update here?)
 *  - decode_block, decode_block_cond
 */
template <class W, /* weight type */
          class BitsetNeighbors /* block manager type */>
struct packed_asymmetric_vertex {
  using iter_type = typename BitsetNeighbors::iter;
  BitsetNeighbors in_nghs;   // copy; not a reference.
  BitsetNeighbors out_nghs;  // copy; not a reference.

  packed_asymmetric_vertex(BitsetNeighbors&& in_nghs,
                           BitsetNeighbors&& out_nghs)
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

  asymmetric_graph<vertex_type, W>* GA;  // the underlying (immutable) graph
  using vertex = vertex_type<W>;
  using weight_type = W;
  using E = typename vertex::edge_type;

  parlay::sequence<vtx_info>
      in_VI;  // metadata for each vertex into the in_blocks structure
  parlay::sequence<uint8_t>
      in_blocks;  // stores the block information for each vertex's in_neighbors
  size_t in_blocks_mem;

  parlay::sequence<vtx_info>
      out_VI;  // metadata for each vertex into the out_blocks structure
  parlay::sequence<uint8_t> out_blocks;  // stores the block information for
                                         // each vertex's out_neighbors
  size_t out_blocks_mem;

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

  template <
      bool bool_enable = true,
      typename std::enable_if<
          std::is_same<vertex, cav_bytepd_amortized<W>>::value && bool_enable,
          int>::type = 0>
  void init_block_size_and_metadata() {
    using neighborhood_type = compressed_bitset_neighbors<vertex_type, W>;
    using metadata = typename neighborhood_type::metadata;
    bs = neighborhood_type::edges_per_block;
    bs_in_bytes = bitsets::bitset_bs_bytes(bs);
    metadata_size = sizeof(metadata);
  }

  asymmetric_packed_graph() : n(0), m(0) {}

  asymmetric_packed_graph(asymmetric_graph<vertex_type, W>& GA)
      : n(GA.n), m(GA.m), GA(&GA) {
    init_block_size_and_metadata();  // conditioned on vertex type; initializes
                                     // bs, bs_in_bytes, metadata_size
    std::tie(in_VI, in_blocks, in_blocks_mem) =
        init_block_memory(GA, bs, bs_in_bytes, in_edges);
    std::tie(out_VI, out_blocks, out_blocks_mem) =
        init_block_memory(GA, bs, bs_in_bytes);
  }

  template <class P>
  asymmetric_packed_graph(asymmetric_graph<vertex_type, W>& GA, P vtx_pred)
      : n(GA.n), m(GA.m), GA(&GA) {
    init_block_size_and_metadata();  // conditioned on vertex type; initializes
                                     // bs, bs_in_bytes, metadata_size
    std::tie(in_VI, in_blocks, in_blocks_mem) =
        init_block_memory(GA, bs, bs_in_bytes, vtx_pred, in_edges);
    std::tie(out_VI, out_blocks, out_blocks_mem) =
        init_block_memory(GA, bs, bs_in_bytes, vtx_pred);
    auto degree_seq = parlay::delayed_seq<size_t>(
        n, [&](size_t i) { return out_VI[i].vtx_degree; });
    m = parlay::reduce(degree_seq);
  }

  inline size_t out_degree(uintE v) { return out_VI[v].vtx_degree; }

  inline size_t in_degree(uintE v) { return in_VI[v].vtx_degree; }

  // Uncompressed: get_vertex
  template <
      bool bool_enable = true,
      typename std::enable_if<
          std::is_same<vertex, asymmetric_vertex<W>>::value && bool_enable,
          int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using neighborhood_type = uncompressed_bitset_neighbors<vertex_type, W>;

    auto in_vtx_data = GA->v_in_data[v];
    uintE in_original_degree = in_vtx_data.degree;
    uintT in_offset = in_vtx_data.offset;

    auto out_vtx_data = GA->v_out_data[v];
    uintE out_original_degree = out_vtx_data.degree;
    uintT out_offset = out_vtx_data.offset;

    auto in_neighborhood =
        neighborhood_type(v, in_blocks.begin(), in_original_degree,
                          in_VI.begin(), GA->in_edges + in_offset);
    auto out_neighborhood =
        neighborhood_type(v, out_blocks.begin(), out_original_degree,
                          out_VI.begin(), GA->out_edges + out_offset);
    return packed_asymmetric_vertex<W, neighborhood_type>(
        std::move(in_neighborhood), std::move(out_neighborhood));
  }

  // Compressed : get_vertex
  template <
      bool bool_enable = true,
      typename std::enable_if<
          std::is_same<vertex, cav_bytepd_amortized<W>>::value && bool_enable,
          int>::type = 0>
  __attribute__((always_inline)) inline auto get_vertex(uintE v) {
    using neighborhood_type = compressed_bitset_neighbors<vertex_type, W>;

    auto in_vtx_data = GA->v_in_data[v];
    uintE in_original_degree = in_vtx_data.degree;
    uintT in_offset = in_vtx_data.offset;

    auto out_vtx_data = GA->v_out_data[v];
    uintE out_original_degree = out_vtx_data.degree;
    uintT out_offset = out_vtx_data.offset;

    auto in_neighborhood =
        neighborhood_type(v, in_blocks.begin(), in_original_degree,
                          in_VI.begin(), GA->in_edges + in_offset);
    auto out_neighborhood =
        neighborhood_type(v, out_blocks.begin(), out_original_degree,
                          out_VI.begin(), GA->out_edges + out_offset);
    return packed_asymmetric_vertex<W, neighborhood_type>(
        std::move(in_neighborhood), std::move(out_neighborhood));
  }

  template <class P>
  uintE packNeighbors(uintE id, P& p, uint8_t* tmp) {
    auto vtx = get_vertex(id);
    return vtx.out_neighbors().pack(p, tmp, false);
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    parallel_for(0, n,
                 [&](size_t v) {
                   auto vert_v = get_vertex(v);
                   vert_v.out_neighbors().map(f, parallel_inner_map);
                 },
                 1);
  }

  // Clears all vertices (zeros out their in-/out-degree) that do not satisfy
  // the predicate P.
  template <class P>
  void clear_vertices(P pred_v) {
    auto deg_imap = parlay::delayed_seq<uintT>(n, [&](size_t i) {
      if (!pred_v(i)) {  // delete
        clear_vertex(i);
      }
      return out_VI[i].vtx_degree;
    });
    std::cout << "Before clearing, m = " << m << std::endl;
    m = parlay::reduce(deg_imap);
    std::cout << "After clearing, m = " << m << std::endl;
  }

  // Clears all vertices (zeros out their in-/out-degree) that do not satisfy
  // the predicate P.
  template <class Seq, class P>
  void clear_vertices(Seq& S, P pred_v) {
    auto deg_imap = parlay::delayed_seq<uintT>(S.size(), [&](size_t i) {
      uintE v = S[i];
      uintE ret = 0;
      if (!pred_v(v)) {  // delete
        ret = out_VI[v].vtx_degree;
        clear_vertex(v);
      }
      return ret;
    });
    std::cout << "Before clearing, m = " << m << std::endl;
    m -= parlay::reduce(deg_imap);
    std::cout << "After clearing, m = " << m << std::endl;
  }

  template <class P, class Seq>
  void filter_graph(P& pred_f, Seq& S) {
    // TODO: do allocations, but in a (medium) constant number of allocations.

    auto reverse_pred = [&](const uintE& u, const uintE& v, const W& wgh) {
      return pred_f(v, u, wgh);
    };

    auto degree_seq = parlay::delayed_seq<size_t>(S.size(), [&](size_t i) {
      uintE v = S[i];
      return out_degree(v);
    });
    size_t pre_degrees = parlay::reduce(degree_seq);

    parallel_for(0, S.size(),
                 [&](size_t i) {
                   uintE v = S[i];
                   auto vtx = get_vertex(v);
                   // Pack-out both the in- and out-edges.
                   if (vtx.out_degree() > 0) {
                     vtx.out_neighbors().pack(pred_f, /* tmp = */ nullptr,
                                              /* parallel = */ true,
                                              /* flags = */ compact_blocks);
                   }
                   if (vtx.in_degree() > 0) {
                     vtx.in_neighbors().pack(reverse_pred, /* tmp = */ nullptr,
                                             /* parallel = */ true,
                                             /* flags = */ compact_blocks);
                   }
                 },
                 128);

    size_t post_degrees = parlay::reduce(degree_seq);
    size_t delta = pre_degrees - post_degrees;
    m -= delta;
    std::cout << "# Filtered graph, new m = " << m << std::endl;
  }

 private:
  inline void clear_vertex(uintE v) {
    out_VI[v].vtx_degree = 0;
    in_VI[v].vtx_degree = 0;
  }
};

// Used to infer template arguments
template <template <class W> class vertex_type, class W>
auto build_asymmetric_packed_graph(asymmetric_graph<vertex_type, W>& GA) {
  return asymmetric_packed_graph<vertex_type, W>(GA);
}

// Used to infer template arguments
template <template <class W> class vertex_type, class W, class Pred>
auto build_asymmetric_packed_graph(asymmetric_graph<vertex_type, W>& GA,
                                   Pred vtx_pred) {
  return asymmetric_packed_graph<vertex_type, W>(GA, vtx_pred);
}

// The predicate pred_f is defined over the _out-edges_ of the graph.
template <template <class W> class vertex_type, class W, class P>
asymmetric_packed_graph<vertex_type, W> filter_graph(
    asymmetric_graph<vertex_type, W>& G, P& pred_f) {
  // TODO: do allocations, but in a (medium) constant number of allocations.
  auto GA = asymmetric_packed_graph<vertex_type, W>(G);

  auto reverse_pred = [&](const uintE& u, const uintE& v, const W& wgh) {
    return pred_f(v, u, wgh);
  };
  {
    parallel_for(0, G.n,
                 [&](size_t v) {
                   auto vtx = GA.get_vertex(v);
                   // Pack-out both the in- and out-edges.
                   if (vtx.out_degree() > 0) {
                     vtx.out_neighbors().pack(pred_f, /* tmp = */ nullptr,
                                              /* parallel = */ true,
                                              /* flags = */ compact_blocks);
                   }
                   if (vtx.in_degree() > 0) {
                     vtx.in_neighbors().pack(reverse_pred, /* tmp = */ nullptr,
                                             /* parallel = */ true,
                                             /* flags = */ compact_blocks);
                   }
                 },
                 1);
  }
  auto degree_seq = parlay::delayed_seq<size_t>(
      GA.n, [&](size_t i) { return GA.out_degree(i); });
  auto new_m = parlay::reduce(degree_seq);
  GA.m = new_m;
  std::cout << "# Returning new asymmetric packed graph, new m = " << new_m
            << std::endl;
  return GA;
}

template <template <class W> class vertex_type, class W, class P>
void filter_graph(asymmetric_packed_graph<vertex_type, W>& GA, P& pred_f) {
  // TODO: do allocations, but in a (medium) constant number of allocations.

  auto reverse_pred = [&](const uintE& u, const uintE& v, const W& wgh) {
    return pred_f(v, u, wgh);
  };
  {
    parallel_for(0, GA.n,
                 [&](size_t v) {
                   auto vtx = GA.get_vertex(v);
                   // Pack-out both the in- and out-edges.
                   if (vtx.out_degree() > 0) {
                     vtx.out_neighbors().pack(pred_f, /* tmp = */ nullptr,
                                              /* parallel = */ true,
                                              /* flags = */ compact_blocks);
                   }
                   if (vtx.in_degree() > 0) {
                     vtx.in_neighbors().pack(reverse_pred, /* tmp = */ nullptr,
                                             /* parallel = */ true,
                                             /* flags = */ compact_blocks);
                   }
                 },
                 1);
  }
  auto degree_seq = parlay::delayed_seq<size_t>(
      GA.n, [&](size_t i) { return GA.out_degree(i); });
  auto new_m = parlay::reduce(degree_seq);
  GA.m = new_m;
  std::cout << "# Packing asymmetric packed graph: new m = " << new_m
            << std::endl;
}

template <template <class W> class vertex_type, class W, class Seq, class P>
void filter_graph(asymmetric_packed_graph<vertex_type, W>& GA, Seq& S,
                  P& pred_f) {
  // TODO: do allocations, but in a (medium) constant number of allocations.

  auto reverse_pred = [&](const uintE& u, const uintE& v, const W& wgh) {
    return pred_f(v, u, wgh);
  };
  {
    parallel_for(0, GA.n,
                 [&](size_t v) {
                   auto vtx = GA.get_vertex(v);
                   // Pack-out both the in- and out-edges.
                   if (vtx.out_degree() > 0) {
                     vtx.out_neighbors().pack(pred_f, /* tmp = */ nullptr,
                                              /* parallel = */ true,
                                              /* flags = */ compact_blocks);
                   }
                   if (vtx.in_degree() > 0) {
                     vtx.in_neighbors().pack(reverse_pred, /* tmp = */ nullptr,
                                             /* parallel = */ true,
                                             /* flags = */ compact_blocks);
                   }
                 },
                 1);
  }
  auto degree_seq = parlay::delayed_seq<size_t>(
      GA.n, [&](size_t i) { return GA.out_degree(i); });
  auto new_m = parlay::reduce(degree_seq);
  GA.m = new_m;
  std::cout << "# Packing asymmetric packed graph: new m = " << new_m
            << std::endl;
}

}  // namespace sage
}  // namespace gbbs
