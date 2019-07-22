#pragma once

#include "macros.h"
#include "bitset.h"
#include "encodings/byte_pd_amortized.h"


// block degree computable by differencing two starts.
template <class E>
struct vtx_info {
  uintE vtx_degree; // vertex's (packed) degree.
  uintE vtx_num_blocks; // number of blocks associated with v
  size_t vtx_block_offset; // pointer into the block structure

  E* vtx_edges; // pointer to the original edges (prevents one random read)

  vtx_info(uintE degree, uintE num_blocks, size_t block_offset, E* edges) :
    vtx_degree(degree), vtx_num_blocks(num_blocks),
    vtx_block_offset(block_offset), vtx_edges(edges) {}
};

template <template <class W> class vertex, class W>
struct sym_bitset_manager {
  using WV = vertex<W>;
  using E = typename WV::E;
  using metadata = bitsets::metadata;

  uintE vtx_id;
  uintE vtx_degree;
  uintE vtx_num_blocks;
  uint8_t* blocks_start;
  uint8_t* block_data_start;

  vtx_info<E>* v_infos;

  static constexpr uintE edges_per_block = 2048;
  static constexpr uintE bytes_per_block = edges_per_block/8 + sizeof(metadata);
  static constexpr uintE bitset_bytes_per_block =
      bytes_per_block - sizeof(metadata);

  E* e0;
  sym_bitset_manager(const uintE vtx_id, uint8_t* blocks, vtx_info<E>* v_infos) :
        vtx_id(vtx_id), v_infos(v_infos) {
    auto& v_info = v_infos[vtx_id];
    vtx_degree = v_info.vtx_degree;
    vtx_num_blocks = v_info.vtx_num_blocks;
    e0 = v_info.vtx_edges;

    blocks_start = blocks + v_info.vtx_block_offset;
    block_data_start = blocks_start + (vtx_num_blocks*sizeof(metadata));
  }

  __attribute__((always_inline)) inline E* get_edges() {
    return e0; // todo: NVM
  }

  __attribute__((always_inline)) inline uintE get_degree() {
    return vtx_degree;
  }

  __attribute__((always_inline)) inline size_t num_blocks() {
    return vtx_num_blocks;
  }

  __attribute__((always_inline)) inline uintE block_degree(uintE block_id) {
    return bitsets::block_degree(blocks_start, block_id, vtx_num_blocks, vtx_degree);
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block(uintE block_id, F f) {
    metadata* block_metadata = (metadata*)blocks_start;
    uintE offset = block_metadata[block_id].offset;
    uintE orig_block_num = block_metadata[block_id].block_num;

    uint8_t* block_bits = block_data_start + bitset_bytes_per_block*block_id;

    uintE block_start = orig_block_num*edges_per_block;
    uintE block_end = std::min(block_start + edges_per_block, vtx_degree);

    // the following matches the perf of noop
    // uintE block_start = block_id*edges_per_block;
    // uintE block_end = std::min(block_start + edges_per_block, vtx_degree);
    // uintE offset = block_start;

    E* e = get_edges();

    // This is one way of decoding (check bit at a time). The other way is to
    // use a fetch_next_bit function.
    // Probably also faster to have a look-up table on the byte
    for (size_t k=0; k<(block_end - block_start); k++) {
      if (bitsets::is_bit_set(block_bits, k)) { // check if the k-th bit is set
        auto& ee = e[block_start + k]; // if so, fetch the k-th edge
        f(std::get<0>(ee), std::get<1>(ee), offset++); // and apply f with the correct offset
      }
    }
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block_cond(uintE block_id, F f) {
    metadata* block_metadata = (metadata*)blocks_start;
    uintE offset = block_metadata[block_id].offset;
    uintE orig_block_num = block_metadata[block_id].block_num;

    uint8_t* block_bits = block_data_start + bitset_bytes_per_block*block_id;

    uintE block_start = orig_block_num*edges_per_block;
    uintE block_end = std::min(block_start + edges_per_block, vtx_degree);

    // the following matches the perf of noop
    // uintE block_start = block_id*edges_per_block;
    // uintE block_end = std::min(block_start + edges_per_block, vtx_degree);
    // uintE offset = block_start;

    E* e = get_edges();

    // This is one way of decoding (check bit at a time). The other way is to
    // use a fetch_next_bit function.
    // Probably also faster to have a look-up table on the byte
    for (size_t k=0; k<(block_end - block_start); k++) {
      if (bitsets::is_bit_set(block_bits, k)) { // check if the k-th bit is set
        auto& ee = e[block_start + k]; // if so, fetch the k-th edge
        bool ret = f(std::get<0>(ee), std::get<1>(ee), offset++); // and apply f with the correct offset
        if (!ret) break;
      }
    }
  }

  template <class P, class E>
  inline void pack_blocks(uintE vtx_id, P& p, E* tmp, bool parallel) {

//    metadata* block_metadata = (metadata*)blocks_start;
//
//    // 1. pack each block. atomically increment counter to mark newly emptied
//    // blocks.
//    uintE num_empty_blocks = 0;
//    par_for(0, vtx_num_blocks, [&] (size_t block_id) {
//      uintE offset = block_metadata[block_id].offset;
//      uintE orig_block_num = block_metadata[block_id].block_num;
//      uint8_t* block_bits = block_data_start + bitset_bytes_per_block*block_id;
//
//      uintE block_start = orig_block_num*edges_per_block;
//      uintE block_end = std::min(block_start + edges_per_block, vtx_degree);
//
//      E* e = get_edges();
//
//      for (size_t k=0; k<(block_end - block_start); k++) {
//        if (bitsets::is_bit_set(block_bits, k)) { // check if the k-th bit is set
//          auto& ee = e[block_start + k]; // if so, fetch the k-th edge
//          bool ret = f(std::get<0>(ee), std::get<1>(ee), offset++); // and apply f with the correct offset
//          if (!ret) break;
//        }
//      }
//
//    }, parallel);
//
//    // 2. if any blocks became empty, compact the bitset structure. Update
//    // vtx_info for this vertex.

  }

  std::tuple<uintE, W> ith_neighbor(size_t i) {
    E* edges = e0;
#ifdef NVM
    if (numa_node() == 0) edges = e0;
    else edges = e1;
#endif
    return edges[i];
  }
};

template <template <class W> class vertex, class W>
struct compressed_sym_bitset_manager {
  using WV = vertex<W>;
  using E = typename WV::E;
  using metadata = bitsets::metadata;

  uintE vtx_id;
  uintE vtx_degree;
  uintE vtx_num_blocks;
  uint8_t* blocks_start;
  uint8_t* block_data_start;

  vtx_info<E>* v_infos;

  static constexpr uintE edges_per_block = PARALLEL_DEGREE;
  static constexpr uintE bytes_per_block = edges_per_block/8 + sizeof(uintE);
  static constexpr uintE bitset_bytes_per_block =
      bytes_per_block - sizeof(uintE);

  E* e0;
  compressed_sym_bitset_manager(const uintE vtx_id, uint8_t* blocks, vtx_info<E>* v_infos) :
        vtx_id(vtx_id), v_infos(v_infos) {
    auto& v_info = v_infos[vtx_id];
    vtx_degree = v_info.vtx_degree;
    vtx_num_blocks = v_info.vtx_num_blocks;
    e0 = v_info.vtx_edges;

    blocks_start = blocks + v_info.vtx_block_offset;
    block_data_start = blocks_start + (vtx_num_blocks*sizeof(metadata));
  }

  __attribute__((always_inline)) inline E* get_edges() {
    return e0; // todo: NVM
  }

  __attribute__((always_inline)) inline uintE get_degree() {
    return vtx_degree;
  }

  __attribute__((always_inline)) inline size_t num_blocks() {
    return vtx_num_blocks;
  }

  __attribute__((always_inline)) inline uintE block_degree(uintE block_id) {
    return bitsets::block_degree(blocks_start, block_id, vtx_num_blocks, vtx_degree);
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block(uintE block_id, F f) {
    metadata* block_metadata = (metadata*)blocks_start;
    uintE offset = block_metadata[block_id].offset;
    uintE orig_block_num = block_metadata[block_id].block_num;

    uint8_t* block_bits = block_data_start + bitset_bytes_per_block*block_id;

    // the following matches the perf of noop
    // uintE block_start = block_id*edges_per_block;
    // uintE block_end = std::min(block_start + edges_per_block, vtx_degree);
    // uintE offset = block_start;

    E* e = get_edges();

    // This is one way of decoding (check bit at a time). The other way is to
    // use a fetch_next_bit function.
    // Probably also faster to have a look-up table on the byte

    size_t k = 0;
    auto t_f = [&] (const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
      if (bitsets::is_bit_set(block_bits, k++)) { // check if the k-th bit is set
        f(ngh, wgh, offset++); // and apply f with the correct offset
      }
    };
    bytepd_amortized::template decode_block<W>(t_f, e, vtx_id, vtx_degree, orig_block_num);
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block_cond(uintE block_id, F f) {
    metadata* block_metadata = (metadata*)blocks_start;
    uintE offset = block_metadata[block_id].offset;
    uintE orig_block_num = block_metadata[block_id].block_num;

    uint8_t* block_bits = block_data_start + bitset_bytes_per_block*block_id;

    // the following matches the perf of noop
    // uintE block_start = block_id*edges_per_block;
    // uintE block_end = std::min(block_start + edges_per_block, vtx_degree);
    // uintE offset = block_start;
    E* e = get_edges();

    size_t k = 0;
    auto t_f = [&] (const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
      if (bitsets::is_bit_set(block_bits, k++)) { // check if the k-th bit is set
        return f(ngh, wgh, offset++); // and apply f with the correct offset
      }
      return true;
    };
    bytepd_amortized::template decode_block_cond<W>(t_f, e, vtx_id, vtx_degree, orig_block_num);
  }

//  std::tuple<uintE, W> ith_neighbor(size_t i) {
//    E* edges = e0;
//#ifdef NVM
//    if (numa_node() == 0) edges = e0;
//    else edges = e1;
//#endif
//    return edges[i];
//  }
};
