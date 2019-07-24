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
  uintE vtx_original_degree;
  uintE vtx_num_blocks;
  uint8_t* blocks_start;
  uint8_t* block_data_start;

  vtx_info<E>* v_infos;

  static constexpr uintE edges_per_block = 256;
  static constexpr uintE bytes_per_block = edges_per_block/8 + sizeof(metadata);
  static constexpr uintE bitset_bytes_per_block =
      bytes_per_block - sizeof(metadata);
  static constexpr uintE kFullBlockPackThreshold = 8;
  static constexpr uintE kBlockAllocThreshold = 100;

  E* e0;
  sym_bitset_manager(const uintE vtx_id, uint8_t* blocks, uintE vtx_original_degree, vtx_info<E>* v_infos) :
        vtx_id(vtx_id), vtx_original_degree(vtx_original_degree), v_infos(v_infos) {
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
    uintE next_block_offset = (block_id == vtx_num_blocks-1) ? vtx_degree : block_metadata[block_id+1].offset;
    uintE block_degree = next_block_offset - offset;
    if (block_degree == 0) { return; }

    uintE orig_block_num = block_metadata[block_id].block_num;

    uint8_t* block_bits = block_data_start + bitset_bytes_per_block*block_id;

    uintE block_start = orig_block_num*edges_per_block;
    uintE block_end = std::min(block_start + edges_per_block, vtx_original_degree);
    assert(block_start < block_end);

    // the following matches the perf of noop
//     uintE block_start = block_id*edges_per_block;
//     uintE block_end = std::min(block_start + edges_per_block, vtx_degree);
//     uintE offset = block_start;

    E* e = get_edges();

    // This is one way of decoding (check bit at a time). The other way is to
    // use a fetch_next_bit function.
    // Probably also faster to have a look-up table on the byte

//    uintE n_bytes = (block_end-block_start) >> 3; // div-by-8
//    for (uintE b=0; b<n_bytes; b++) {
//      if (block_bits[b] > 0) {
//        uint8_t* byte_loc = &(block_bits[b]);
//        uintE idx_loc = b << 3; // * 8
//        for (uintE k=0; k<8; k++) {
//          if (bitsets::is_bit_set(byte_loc, k)) {
//            auto& ee = e[block_start + idx_loc + k];
//            f(std::get<0>(ee), std::get<1>(ee), offset++); // and apply f with the correct offset
//          }
//        }
//      }
//    }

    for (size_t k=0; k<(block_end - block_start); k++) {
//      bool isset = bitsets::is_bit_set(block_bits, k);
//      assert(isset);
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
    uintE next_block_offset = (block_id == vtx_num_blocks-1) ? vtx_degree : block_metadata[block_id+1].offset;
    uintE block_degree = next_block_offset - offset;
    if (block_degree == 0) { return; }

    uintE orig_block_num = block_metadata[block_id].block_num;

    uint8_t* block_bits = block_data_start + bitset_bytes_per_block*block_id;

    uintE block_start = orig_block_num*edges_per_block;
    uintE block_end = std::min(block_start + edges_per_block, vtx_original_degree);

//     the following matches the perf of noop
//     uintE block_start = block_id*edges_per_block;
//     uintE block_end = std::min(block_start + edges_per_block, vtx_degree);
//     uintE offset = block_start;

    E* e = get_edges();

    // This is one way of decoding (check bit at a time). The other way is to
    // use a fetch_next_bit function.
    // Probably also faster to have a look-up table on the byte
    for (size_t k=0; k<(block_end - block_start); k++) {
//      bool isset = bitsets::is_bit_set(block_bits, k);
//      assert(isset);
      if (bitsets::is_bit_set(block_bits, k)) { // check if the k-th bit is set
        auto& ee = e[block_start + k]; // if so, fetch the k-th edge
        bool ret = f(std::get<0>(ee), std::get<1>(ee), offset++); // and apply f with the correct offset
        if (!ret) break;
      }
    }
  }


  // Only called when the discrepency between full and total blocks is large.
  // Specifically, when #full_blocks*kFullBlockPackThreshold >= vtx_num_blocks
  inline size_t repack_blocks_par(bool parallel) {
    uint8_t stk[bytes_per_block*kBlockAllocThreshold]; // temporary space
    uintE int_stk[kBlockAllocThreshold];
    uint8_t* tmp_space = (uint8_t*)stk;
    uintE* tmp_ints = (uintE*)int_stk;
    size_t total_bytes = vtx_num_blocks*bytes_per_block;
    if (vtx_num_blocks > kBlockAllocThreshold) {
      tmp_space = pbbs::new_array_no_init<uint8_t>(total_bytes);
      uintE* tmp_ints = pbbs::new_array_no_init<uintE>(vtx_num_blocks);
    }

    // 1. Copy all data to tmp space
    parallel_for(0, total_bytes, [&] (size_t i) {
      tmp_space[i] = blocks_start[i];
    }, 512); // tune threshold

    // 2. Write 1 to tmp_ints (new_locs) if full, 0 if empty
    auto new_locs = pbbslib::make_sequence(tmp_ints, vtx_num_blocks);
    auto tmp_metadata = (metadata*)tmp_space;
    parallel_for(0, vtx_num_blocks, [&] (size_t block_id) {
      new_locs[block_id] = static_cast<uintE>((tmp_metadata[block_id].offset > 0));
    });

    // 3. Scan new_locs to get new block indices for full blocks
    size_t new_num_blocks = pbbslib::scan_add_inplace(new_locs);

    // 4. Copy saved blocks to new positions.
    auto real_metadata = (metadata*)blocks_start;
    auto real_block_data = blocks_start + (new_num_blocks*sizeof(metadata));
    auto tmp_block_data = tmp_space + (vtx_num_blocks*sizeof(metadata));
    parallel_for(0, vtx_num_blocks, [&] (size_t block_id) {
      uintE block_entries = tmp_metadata[block_id].offset;
      if (block_entries > 0) { // live
        uintE new_block_id = new_locs[block_id];
        // (a) copy metadata
        real_metadata[new_block_id] = tmp_metadata[block_id];
        uintE orig_block_num = real_metadata[new_block_id].block_num;

        uint8_t* tmp_block_bits = tmp_block_data + bitset_bytes_per_block*block_id;
        uint8_t* real_block_bits = real_block_data + bitset_bytes_per_block*new_block_id;

        uintE block_start = orig_block_num*edges_per_block;
        uintE block_end = std::min(block_start + edges_per_block, vtx_original_degree);

        size_t this_block_size = block_end - block_start;
        assert(this_block_size % 8 == 0);
        size_t bytes_to_copy = this_block_size / 8;

        // (b) copy bitset data
        for (size_t i=0; i<bytes_to_copy; i++) {
          real_block_bits[i] = tmp_block_bits[i];
        }
      }
    });



    if (vtx_num_blocks > kBlockAllocThreshold) {
      pbbs::free_array(tmp_space);
      pbbs::free_array(tmp_ints);
    }
  }

  // P : (uintE, uintE, wgh) -> bool
  // only keep edges s.t. P(...) = true.
  template <class P>
  inline size_t pack_blocks(uintE vtx_id, P& p, bool parallel) {
    metadata* block_metadata = (metadata*)blocks_start;

    // 1. pack each block
    par_for(0, vtx_num_blocks, 4, [&] (size_t block_id) {
      uintE offset = block_metadata[block_id].offset;
      uintE orig_block_num = block_metadata[block_id].block_num;
      uint8_t* block_bits = block_data_start + bitset_bytes_per_block*block_id;

      uintE block_start = orig_block_num*edges_per_block;
      uintE block_end = std::min(block_start + edges_per_block, vtx_original_degree);

      E* e = get_edges();

      size_t live_edges = 0;
      for (size_t k=0; k<(block_end - block_start); k++) {
        if (bitsets::is_bit_set(block_bits, k)) { // k-th edge is present in G
          auto& ee = e[block_start + k]; // fetch the k-th edge
          if (!p(vtx_id, std::get<0>(ee), std::get<1>(ee))) {
            // unset the k-th bit
            bitsets::flip_bit(block_bits, k);
          } else {
            // otherwise increment the count of live edges
            live_edges++;
          }
        }
      }
      // Temporarily store #live_edges in offset positions.
      block_metadata[block_id].offset = live_edges;
    }, parallel);

    // 2. Reduce to get the #empty_blocks
    auto full_block_seq = pbbslib::make_sequence<size_t>(vtx_num_blocks, [&] (size_t i) {
      return static_cast<size_t>(block_metadata[i].offset > 0);
    });
    size_t full_blocks = pbbslib::reduce_add(full_block_seq);

    if (vtx_num_blocks > 500) {
      cout << "full_blocks = " << full_blocks << " vtx_blocks = " << vtx_num_blocks << endl;
    }

//    if (full_blocks*kFullBlockPackThreshold >= vtx_num_blocks) {
//      repack_blocks_par(parallel);
//      return;
//    }
    // Otherwise, not enough empty blocks to warrant a full re-pack. Just update
    // offset values.

    auto ptr_seq = pbbs::indirect_value_seq<uintE>(vtx_num_blocks, [&] (size_t i) {
      return &(block_metadata[i].offset);
    });

    uintE sum = pbbslib::scan_add_inplace(ptr_seq);
    vtx_degree = sum;

    // Update the degree in vtx_info.
    v_infos[vtx_id].vtx_degree = sum;
    return sum; // return the new degree
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
