#pragma once

#include "bitset.h"
#include "encodings/byte_pd_amortized.h"
#include "macros.h"

// block degree computable by differencing two starts.
template <class E>
struct vtx_info {
  uintE vtx_degree;         // vertex's (packed) degree.
  uintE vtx_num_blocks;     // number of blocks associated with v
  size_t vtx_block_offset;  // pointer into the block structure

  E* vtx_edges;  // pointer to the original edges (prevents one random read)

  vtx_info(uintE degree, uintE num_blocks, size_t block_offset, E* edges)
      : vtx_degree(degree),
        vtx_num_blocks(num_blocks),
        vtx_block_offset(block_offset),
        vtx_edges(edges) {}
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

  static constexpr uintE edges_per_block = 2048;
  static constexpr uintE bytes_per_block =
      edges_per_block / 8 + sizeof(metadata);
  static constexpr uintE bitset_bytes_per_block =
      bytes_per_block - sizeof(metadata);
  static constexpr uintE kFullBlockPackThreshold = 2;
  static constexpr uintE kBlockAllocThreshold = 20;

  E* e0;
  sym_bitset_manager(const uintE vtx_id, uint8_t* blocks,
                     uintE vtx_original_degree, vtx_info<E>* v_infos)
      : vtx_id(vtx_id),
        vtx_original_degree(vtx_original_degree),
        v_infos(v_infos) {
    auto& v_info = v_infos[vtx_id];
    vtx_degree = v_info.vtx_degree;
    vtx_num_blocks = v_info.vtx_num_blocks;
    e0 = v_info.vtx_edges;

    blocks_start = blocks + v_info.vtx_block_offset;
    block_data_start = blocks_start + (vtx_num_blocks * sizeof(metadata));
  }

  __attribute__((always_inline)) inline E* get_edges() {
    return e0;  // todo: NVM
  }

  __attribute__((always_inline)) inline uintE get_degree() {
    return vtx_degree;
  }

  __attribute__((always_inline)) inline size_t num_blocks() {
    return vtx_num_blocks;
  }

  __attribute__((always_inline)) inline uintE block_degree(uintE block_id) {
    return bitsets::block_degree(blocks_start, block_id, vtx_num_blocks,
                                 vtx_degree);
  }

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

    // the following matches the perf of noop
    //     uintE block_start = block_id*edges_per_block;
    //     uintE block_end = std::min(block_start + edges_per_block,
    //     vtx_degree); uintE offset = block_start;

    E* e = get_edges();

    size_t block_size = block_end - block_start;
    size_t block_size_num_longs = (block_size+64-1)/64;

    // "select" on a block
    uint64_t* long_block_bits = (uint64_t*)block_bits;
    size_t cur_offset = block_start;
    for (size_t idx = 0; idx < block_size_num_longs; idx++) {
      uint64_t cur_long = long_block_bits[idx];
      size_t cnt = _mm_popcnt_u64(cur_long); // #bits set to one
      if (cnt > 0) {
        for (size_t i=0; i<cnt; i++) {
          unsigned select_idx = _tzcnt_u64(cur_long); // select64_tz(cur_long, 1);
          auto& ee = e[cur_offset + select_idx];
          f(std::get<0>(ee), std::get<1>(ee), offset++);
          assert((cur_long & (1UL << select_idx)) > 0);
          cur_long = _blsr_u64(cur_long);
        }
      }
      cur_offset += 64; // next long
    }

//    for (size_t k=0; k<(block_end - block_start); k++) {
//      if (bitsets::is_bit_set(block_bits, k)) { // check if the k-th bit is set
//        auto& ee = e[block_start + k]; // if so, fetch the k-th edge
//        f(std::get<0>(ee), std::get<1>(ee), offset++); // and apply f with the correct offset
//      }
//    }
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

    //     the following matches the perf of noop
    //     uintE block_start = block_id*edges_per_block;
    //     uintE block_end = std::min(block_start + edges_per_block,
    //     vtx_degree); uintE offset = block_start;

    E* e = get_edges();

    size_t block_size = block_end - block_start;
    size_t block_size_num_longs = (block_size+64-1)/64;

    // "select" on a block
    uint64_t* long_block_bits = (uint64_t*)block_bits;
    size_t cur_offset = block_start;
    for (size_t idx = 0; idx < block_size_num_longs; idx++) {
      uint64_t cur_long = long_block_bits[idx];
      size_t cnt = _mm_popcnt_u64(cur_long); // #bits set to one
      if (cnt > 0) {
        for (size_t i=0; i<cnt; i++) {
          unsigned select_idx = _tzcnt_u64(cur_long); // index of first nz bit
          auto& ee = e[cur_offset + select_idx];
          if (!f(std::get<0>(ee), std::get<1>(ee), offset++)) {
            return;
          }
          assert((cur_long & (1UL << select_idx)) > 0);
          cur_long = _blsr_u64(cur_long); // reset lowest bit
        }
      }
      cur_offset += 64; // next long
    }

//    // This is one way of decoding (check bit at a time). The other way is to
//    // use a fetch_next_bit function.
//    // Probably also faster to have a look-up table on the byte
//    for (size_t k = 0; k < (block_end - block_start); k++) {
//      //      bool isset = bitsets::is_bit_set(block_bits, k);
//      //      assert(isset);
//      if (bitsets::is_bit_set(block_bits, k)) {  // check if the k-th bit is set
//        auto& ee = e[block_start + k];           // if so, fetch the k-th edge
//        bool ret = f(std::get<0>(ee), std::get<1>(ee),
//                     offset++);  // and apply f with the correct offset
//        if (!ret) break;
//      }
//    }
  }

  /* Only called when the discrepency between full and total blocks is large.
   * Specifically, when #full_blocks*kFullBlockPackThreshold >= vtx_num_blocks
   * Note: Defers to pack(..) to finish updating the degree/scan info */
  inline void repack_blocks_par(bool parallel) {
    uint8_t stk[bytes_per_block * kBlockAllocThreshold];  // temporary space
    uintE int_stk[kBlockAllocThreshold];
    uint8_t* tmp_space = (uint8_t*)stk;
    uintE* tmp_ints = (uintE*)int_stk;
    size_t total_bytes = vtx_num_blocks * bytes_per_block;
    if (vtx_num_blocks > kBlockAllocThreshold) {
      tmp_space = pbbs::new_array_no_init<uint8_t>(total_bytes);
      tmp_ints = pbbs::new_array_no_init<uintE>(vtx_num_blocks);
    }

    // 1. Copy all data to tmp space
    parallel_for(0, total_bytes,
                 [&](size_t i) { tmp_space[i] = blocks_start[i]; },
                 512);  // tune threshold

    // 2. Write 1 to tmp_ints (new_locs) if full, 0 if empty
    auto new_locs = pbbslib::make_sequence(tmp_ints, vtx_num_blocks);
    auto tmp_metadata = (metadata*)tmp_space;
    parallel_for(0, vtx_num_blocks, [&](size_t block_id) {
      new_locs[block_id] =
          static_cast<uintE>((tmp_metadata[block_id].offset > 0));
    });

    // 3. Scan new_locs to get new block indices for full blocks
    size_t new_num_blocks = pbbslib::scan_add_inplace(new_locs);

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
        size_t bytes_to_copy =
            (this_block_size + 8 - 1) / 8;  // ceil(this_block_size/8);

        // (b) copy bitset data
        for (size_t i = 0; i < bytes_to_copy; i++) {
          real_block_bits[i] = tmp_block_bits[i];
        }
      }
    });

    // 5. Update num_blocks info both locally and in v_infos
    uintE old_vtx_num_blocks = vtx_num_blocks;
    vtx_num_blocks = new_num_blocks;
    v_infos[vtx_id].vtx_num_blocks = new_num_blocks;

    if (old_vtx_num_blocks > kBlockAllocThreshold) {
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
    par_for(0, vtx_num_blocks, 1,
            [&](size_t block_id) {
              uintE orig_block_num = block_metadata[block_id].block_num;
              uint8_t* block_bits =
                  block_data_start + bitset_bytes_per_block * block_id;

              uintE block_start = orig_block_num * edges_per_block;
              uintE block_end =
                  std::min(block_start + edges_per_block, vtx_original_degree);

              E* e = get_edges();

              size_t live_edges = 0;
              for (size_t k = 0; k < (block_end - block_start); k++) {
                if (bitsets::is_bit_set(block_bits, k)) {     // k-th edge is present in G
                  auto& ee = e[block_start + k];  // fetch the k-th edge
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
            },
            parallel);

    // 2. Reduce to get the #empty_blocks
    auto full_block_seq =
        pbbslib::make_sequence<size_t>(vtx_num_blocks, [&](size_t i) {
          return static_cast<size_t>(block_metadata[i].offset > 0);
        });
    size_t full_blocks = pbbslib::reduce_add(full_block_seq);

    if (vtx_num_blocks > 500) {
      cout << "full_blocks = " << full_blocks
           << " vtx_blocks = " << vtx_num_blocks << endl;
    }

    if (full_blocks * kFullBlockPackThreshold <= vtx_num_blocks) {
      repack_blocks_par(parallel);
    }

    // Update offset values.
    auto ptr_seq = pbbs::indirect_value_seq<uintE>(
        vtx_num_blocks, [&](size_t i) { return &(block_metadata[i].offset); });
    uintE sum = pbbslib::scan_add_inplace(ptr_seq);
    vtx_degree = sum;

    // Update the degree in vtx_info.
    v_infos[vtx_id].vtx_degree = sum;
    return sum;  // return the new degree
  }

  std::tuple<uintE, W> ith_neighbor(size_t i) {
    E* edges = e0;
#ifdef NVM
    if (numa_node() == 0)
      edges = e0;
    else
      edges = e1;
#endif
    return edges[i];
  }

  struct iter {
    E* edges;
    uintE vtx_id;
    uintE vtx_degree;
    uintE vtx_original_degree;
    uintE vtx_num_blocks;
    uint8_t* blocks_start;
    uint8_t* block_data_start;

    uintE cur_block; // block id
    uintE cur_block_degree; // #live edges
    uintE cur_block_start; // start offset (in edges)
    uintE cur_block_size; // block size (in edges)
    uint8_t* cur_block_bits;

    std::tuple<uintE, W> last_edge;
    uintE proc;
    uintE proc_cur_block;

    iter(E* edges,
        uintE vtx_id,
        uintE vtx_degree,
        uintE vtx_original_degree,
        uintE vtx_num_blocks,
        uint8_t* blocks_start,
        uint8_t* block_data_start) :
        edges(edges),
        vtx_id(vtx_id),
        vtx_degree(vtx_degree),
        vtx_original_degree(vtx_original_degree),
        vtx_num_blocks(vtx_num_blocks),
        blocks_start(blocks_start),
        block_data_start(block_data_start),
        proc(0) {
      proc = 0;
      if (vtx_degree > 0) {
        find_next_nonempty_block(/* on_initialization = */ true);
        next(); // sets last_edge
      }
    }

// unsigned short __builtin_ia32_lzcnt_u16(unsigned short);
// unsigned int __builtin_ia32_lzcnt_u32(unsigned int);
// unsigned long long __builtin_ia32_lzcnt_u64 (unsigned long long);

    // precondition: there is a subsequent non-empty block
    inline void find_next_nonempty_block(bool on_initialization=false) {
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

      if (cur_block_degree == 0) {
        find_next_nonempty_block();
      } else { // >= 1 live edges this block. reset block ctrs
        proc_cur_block = 0;

        // set cur_block_start and cur_block_size
        uintE orig_block_num = block_metadata[cur_block].block_num;
        cur_block_start = orig_block_num * edges_per_block;
        uintE cur_block_end =
          std::min(cur_block_start + edges_per_block, vtx_original_degree);
        cur_block_size = cur_block_end - cur_block_start;

        // set cur_block_bits
        cur_block_bits = block_data_start + bitset_bytes_per_block * cur_block;
      }
    }

    __attribute__((always_inline)) inline uintE degree() { return vtx_degree; }
    __attribute__((always_inline)) inline std::tuple<uintE, W> cur() { return last_edge; }
    __attribute__((always_inline)) inline void next() {
      while (true) {
        while (proc_cur_block < cur_block_size) {
          if ((proc_cur_block & ((1 << 16) - 1)) == 0) { // % 8
            uintE idx = proc_cur_block >> 3; // byte
            if (cur_block_bits[idx] == 0) {
              if (cur_block_bits[idx+1] == 0) {
                proc_cur_block += 16;
                continue;
              }
              proc_cur_block += 8;
              continue;
            }
          } else if ((proc_cur_block & 0x7) == 0) { // % 8
            if (cur_block_bits[proc_cur_block >> 3] == 0) {
              proc_cur_block += 8;
              continue;
            }
          }
          // cout << "proc_cur_block = " << proc_cur_block << endl;
          if (bitsets::is_bit_set(cur_block_bits, proc_cur_block)) {
            // cout << "bit set for proc_cur_block = " << proc_cur_block << endl;
            last_edge = edges[cur_block_start + proc_cur_block];
            proc++;
            proc_cur_block++;
            return;
          }
          proc_cur_block++;
        }
        // otherwise: finished this block, find next non_empty and call next()
        find_next_nonempty_block();
      }
    }

    __attribute__((always_inline)) inline bool has_next() {
      return proc < vtx_degree;
    }
  };

  auto get_iter() {
    return iter(get_edges(),
        vtx_id,
        vtx_degree,
        vtx_original_degree,
        vtx_num_blocks,
        blocks_start,
        block_data_start);
  }

};

template <template <class W> class vertex, class W>
struct compressed_sym_bitset_manager {
  using WV = vertex<W>;
  using E = typename WV::E;
  using metadata = bitsets::metadata;

  uintE vtx_id;
  uintE vtx_degree;
  uintE vtx_original_degree;
  uintE vtx_num_blocks;
  uint8_t* blocks_start;
  uint8_t* block_data_start;
  static constexpr uintE kFullBlockPackThreshold = 4;
  static constexpr uintE kBlockAllocThreshold = 20;

  vtx_info<E>* v_infos;

  static constexpr uintE edges_per_block = PARALLEL_DEGREE;
  static constexpr uintE bytes_per_block = edges_per_block / 8 + sizeof(metadata);
  static constexpr uintE bitset_bytes_per_block =
      bytes_per_block - sizeof(metadata);

  E* e0;
  compressed_sym_bitset_manager(const uintE vtx_id, uint8_t* blocks,
                                uintE vtx_original_degree, vtx_info<E>* v_infos)
      : vtx_id(vtx_id),
      vtx_original_degree(vtx_original_degree),
      v_infos(v_infos) {
    auto& v_info = v_infos[vtx_id];
    vtx_degree = v_info.vtx_degree;
    vtx_num_blocks = v_info.vtx_num_blocks;
    e0 = v_info.vtx_edges;

    blocks_start = blocks + v_info.vtx_block_offset;
    block_data_start = blocks_start + (vtx_num_blocks * sizeof(metadata));
  }

  __attribute__((always_inline)) inline E* get_edges() {
    return e0;  // todo: NVM
  }

  __attribute__((always_inline)) inline uintE get_degree() {
    return vtx_degree;
  }

  __attribute__((always_inline)) inline size_t num_blocks() {
    return vtx_num_blocks;
  }

  __attribute__((always_inline)) inline uintE block_degree(uintE block_id) {
    return bitsets::block_degree(blocks_start, block_id, vtx_num_blocks,
                                 vtx_degree);
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block(uintE block_id, F f) {
    metadata* block_metadata = (metadata*)blocks_start;
    uintE offset = block_metadata[block_id].offset;
    uintE orig_block_num = block_metadata[block_id].block_num;

    uint8_t* block_bits = block_data_start + bitset_bytes_per_block * block_id;

    // the following matches the perf of noop
    // uintE block_start = block_id*edges_per_block;
    // uintE block_end = std::min(block_start + edges_per_block, vtx_degree);
    // uintE offset = block_start;

    E* e = get_edges();

    // This is one way of decoding (check bit at a time). The other way is to
    // use a fetch_next_bit function.
    // Probably also faster to have a look-up table on the byte

    size_t k = 0;
    auto t_f = [&](const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
      if (bitsets::is_bit_set(block_bits,
                              k++)) {  // check if the k-th bit is set
        f(ngh, wgh, offset++);         // and apply f with the correct offset
      }
    };
    bytepd_amortized::template decode_block<W>(t_f, e, vtx_id, vtx_original_degree,
                                               orig_block_num);
    debug(uintE block_start = orig_block_num * edges_per_block;
    uintE block_end =
        std::min(block_start + edges_per_block, vtx_original_degree);
    assert(k == (block_end - block_start)););
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block_cond(uintE block_id,
                                                               F f) {
    metadata* block_metadata = (metadata*)blocks_start;
    uintE offset = block_metadata[block_id].offset;
    uintE orig_block_num = block_metadata[block_id].block_num;

    uint8_t* block_bits = block_data_start + bitset_bytes_per_block * block_id;

    // the following matches the perf of noop
    // uintE block_start = block_id*edges_per_block;
    // uintE block_end = std::min(block_start + edges_per_block, vtx_degree);
    // uintE offset = block_start;
    E* e = get_edges();

    size_t k = 0;
    auto t_f = [&](const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
      if (bitsets::is_bit_set(block_bits,
                              k++)) {  // check if the k-th bit is set
        return f(ngh, wgh, offset++);  // and apply f with the correct offset
      }
      return true;
    };
    bytepd_amortized::template decode_block_cond<W>(t_f, e, vtx_id, vtx_original_degree,
                                                    orig_block_num);
    debug(uintE block_start = orig_block_num * edges_per_block;
    uintE block_end =
        std::min(block_start + edges_per_block, vtx_original_degree);
    assert(k == (block_end - block_start)););
  }

  /* Only called when the discrepency between full and total blocks is large.
   * Specifically, when #full_blocks*kFullBlockPackThreshold >= vtx_num_blocks
   * Note: Defers to pack(..) to finish updating the degree/scan info */
  inline void repack_blocks_par(bool parallel) {
    uint8_t stk[bytes_per_block * kBlockAllocThreshold];  // temporary space
    uintE int_stk[kBlockAllocThreshold];
    uint8_t* tmp_space = (uint8_t*)stk;
    uintE* tmp_ints = (uintE*)int_stk;
    size_t total_bytes = vtx_num_blocks * bytes_per_block;
    if (vtx_num_blocks > kBlockAllocThreshold) {
      tmp_space = pbbs::new_array_no_init<uint8_t>(total_bytes);
      tmp_ints = pbbs::new_array_no_init<uintE>(vtx_num_blocks);
    }

    // 1. Copy all data to tmp space
    parallel_for(0, total_bytes,
                 [&](size_t i) { tmp_space[i] = blocks_start[i]; },
                 512);  // tune threshold

    // 2. Write 1 to tmp_ints (new_locs) if full, 0 if empty
    auto new_locs = pbbslib::make_sequence(tmp_ints, vtx_num_blocks);
    auto tmp_metadata = (metadata*)tmp_space;
    parallel_for(0, vtx_num_blocks, [&](size_t block_id) {
      new_locs[block_id] =
          static_cast<uintE>((tmp_metadata[block_id].offset > 0));
    });

    // 3. Scan new_locs to get new block indices for full blocks
    size_t new_num_blocks = pbbslib::scan_add_inplace(new_locs);

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
        size_t bytes_to_copy =
            (this_block_size + 8 - 1) / 8;  // ceil(this_block_size/8);

        // (b) copy bitset data
        for (size_t i = 0; i < bytes_to_copy; i++) {
          real_block_bits[i] = tmp_block_bits[i];
        }
      }
    });

    // 5. Update num_blocks info both locally and in v_infos
    uintE old_vtx_num_blocks = vtx_num_blocks;
    vtx_num_blocks = new_num_blocks;
    v_infos[vtx_id].vtx_num_blocks = new_num_blocks;

    if (old_vtx_num_blocks > kBlockAllocThreshold) {
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
    par_for(0, vtx_num_blocks, 1,
            [&](size_t block_id) {
              uintE orig_block_num = block_metadata[block_id].block_num;
              uint8_t* block_bits =
                  block_data_start + bitset_bytes_per_block * block_id;


              E* e = get_edges();

              size_t k = 0;
              size_t live_edges = 0;
              auto t_f = [&](const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
                if (bitsets::is_bit_set(block_bits, k)) {  // check if the k-th bit is set
                  if (!p(vtx_id, ngh, wgh)) {
                    // unset the k-th bit
                    bitsets::flip_bit(block_bits, k);
                  } else {
                    // otherwise increment the count of live edges
                    live_edges++;
                  }
                }
                k += 1;
              };
              bytepd_amortized::template decode_block<W>(t_f, e, vtx_id, vtx_original_degree,
                                                         orig_block_num);

              debug(uintE block_start = orig_block_num * edges_per_block;
              uintE block_end =
                  std::min(block_start + edges_per_block, vtx_original_degree);
              assert(k == (block_end - block_start)););

              // Temporarily store #live_edges in offset positions.
              block_metadata[block_id].offset = live_edges;
            },
            parallel);

      // 2. Reduce to get the #empty_blocks
      auto full_block_seq =
          pbbslib::make_sequence<size_t>(vtx_num_blocks, [&](size_t i) {
            return static_cast<size_t>(block_metadata[i].offset > 0);
          });
      size_t full_blocks = pbbslib::reduce_add(full_block_seq);

//    if (vtx_num_blocks > 500) {
//      cout << "full_blocks = " << full_blocks
//           << " vtx_blocks = " << vtx_num_blocks << endl;
//    }

    if (full_blocks * kFullBlockPackThreshold <= vtx_num_blocks) {
      repack_blocks_par(parallel);
    }

    // Update offset values.
    auto ptr_seq = pbbs::indirect_value_seq<uintE>(
        vtx_num_blocks, [&](size_t i) { return &(block_metadata[i].offset); });
    uintE sum = pbbslib::scan_add_inplace(ptr_seq);
    vtx_degree = sum;

    // Update the degree in vtx_info.
    v_infos[vtx_id].vtx_degree = sum;
    return sum;  // return the new degree
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
