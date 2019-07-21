#pragma once

#include "macros.h"
#include "bitset.h"
#include "encodings/byte_pd_amortized.h"




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

  static constexpr uintE edges_per_block = 128;
  static constexpr uintE bytes_per_block = edges_per_block/8 + sizeof(uintE);
  static constexpr uintE bitset_bytes_per_block =
      bytes_per_block - sizeof(uintE);

  E* e0;
  sym_bitset_manager(const uintE vtx_id, const uintE vtx_degree,
      const uintE vtx_num_blocks, uint8_t* blocks_start, E* e0) :
        vtx_id(vtx_id),
        vtx_degree(vtx_degree),
        vtx_num_blocks(vtx_num_blocks),
        blocks_start(blocks_start),
        e0(e0) {
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
  __attribute__((always_inline)) inline void decode_block_cond(uintE block_num, F f) {


//    uintE block_start = block_num*block_size;
//    uintE block_end = std::min(block_start + block_size, degree);
//    E* e = e0;
//
//    for (size_t i=block_start; i<block_end; i++) {
//      auto& ee = e[i];
//      bool ret = f(std::get<0>(ee), std::get<1>(ee), i);
//      if (!ret) break;
//    }
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
