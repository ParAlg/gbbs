#pragma once

#include "bitset.h"
#include "encodings/byte_pd_amortized.h"
#include "macros.h"
#include "flags.h"

/* Note that the block degree computable by differencing two starts. */
struct vtx_info {
  /* vertex's (packed) degree. */
  uintE vtx_degree;
  /* number of blocks associated with v */
  uintE vtx_num_blocks;
  /* pointer into the block structure */
  size_t vtx_block_offset;

  vtx_info(uintE degree, uintE num_blocks, size_t block_offset)
      : vtx_degree(degree),
        vtx_num_blocks(num_blocks),
        vtx_block_offset(block_offset) {}
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

  vtx_info* v_infos;

  static constexpr uintE edges_per_block = 256;
  static constexpr uintE bytes_per_block =
      edges_per_block / 8 + sizeof(metadata);
  static constexpr uintE bitset_bytes_per_block =
      bytes_per_block - sizeof(metadata);
  static constexpr uintE kFullBlockPackThreshold = 2;
  static constexpr uintE kBlockAllocThreshold = 20;

  E* e0;
#ifdef NVM
  E* e1;
#endif

  sym_bitset_manager(const uintE vtx_id, uint8_t* blocks,
                     uintE vtx_original_degree, vtx_info* v_infos,
#ifndef NVM
                     E* e0
#else
                     E* e0, E* e1
#endif
                     )
      : vtx_id(vtx_id),
        vtx_original_degree(vtx_original_degree),
        v_infos(v_infos),
#ifndef NVM
        e0(e0)
#else
        e0(e0), e1(e1)
#endif
        {
    auto& v_info = v_infos[vtx_id];
    vtx_degree = v_info.vtx_degree;
    vtx_num_blocks = v_info.vtx_num_blocks;

    blocks_start = blocks + v_info.vtx_block_offset;
    block_data_start = blocks_start + (vtx_num_blocks * sizeof(metadata));
  }

  __attribute__((always_inline)) inline E* get_edges() {
#ifndef NVM
    return e0;
#else
    if (numanode() == 0) {
      return e0;
    } else {
      return e1;
    }
#endif
  }

  __attribute__((always_inline)) inline uintE get_degree() {
    return vtx_degree;
  }

  __attribute__((always_inline)) inline size_t num_blocks() {
    return vtx_num_blocks;
  }

//  uintE vtx_num_blocks;     // number of blocks associated with v
//  size_t vtx_block_offset;  // pointer into the block structure
//
//  E* vtx_edges;  // pointer to the original edges (prevents one random read)

  inline void clear_vertex() {
    vtx_degree = 0;
    vtx_num_blocks = 0;
    v_infos[vtx_id].vtx_num_blocks = 0;
    v_infos[vtx_id].vtx_degree = 0;
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

    E* e = get_edges();

    size_t block_size = block_end - block_start;
    size_t block_size_num_longs = (block_size+64-1)/64;

    // "select" on a block
    uint64_t* long_block_bits = (uint64_t*)block_bits;
    size_t cur_offset = block_start;
    for (size_t idx = 0; idx < block_size_num_longs; idx++) {
      uint64_t cur_long = long_block_bits[idx];
      while (cur_long > 0) {
        unsigned select_idx = _tzcnt_u64(cur_long); // #trailing zeros in cur_long
        auto& ee = e[cur_offset + select_idx];
        f(std::get<0>(ee), std::get<1>(ee), offset++);
        assert((cur_long & (1UL << select_idx)) > 0);
        cur_long = _blsr_u64(cur_long); // clears lowest set bit
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

    E* e = get_edges();

    size_t block_size = block_end - block_start;
    size_t block_size_num_longs = (block_size+64-1)/64;

    // "select" on a block
    uint64_t* long_block_bits = (uint64_t*)block_bits;
    size_t cur_offset = block_start;
    for (size_t idx = 0; idx < block_size_num_longs; idx++) {
      uint64_t cur_long = long_block_bits[idx];
      while (cur_long > 0) {
        unsigned select_idx = _tzcnt_u64(cur_long); // index of first nz bit
        auto& ee = e[cur_offset + select_idx];
        if (!f(std::get<0>(ee), std::get<1>(ee), offset++)) {
          return;
        }
        assert((cur_long & (1UL << select_idx)) > 0);
        cur_long = _blsr_u64(cur_long); // reset lowest bit
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
  inline void repack_blocks_par(uint8_t* tmp, bool parallel) {
    uint8_t stk[bytes_per_block * kBlockAllocThreshold];  // temporary space
    uintE int_stk[kBlockAllocThreshold];
    uint8_t* tmp_space = (uint8_t*)stk;
    uintE* tmp_ints = (uintE*)int_stk;
    size_t total_bytes = vtx_num_blocks * bytes_per_block;
    if ((tmp == nullptr) && (vtx_num_blocks > kBlockAllocThreshold)) {
      tmp_space = pbbs::new_array_no_init<uint8_t>(total_bytes);
      tmp_ints = pbbs::new_array_no_init<uintE>(vtx_num_blocks);
    }

    // caller supplies:
    // vtx_num_blocks*sizeof(uintE) +
    // vtx_num_blocks*bytes_per_block

    if (tmp) {
      tmp_space = tmp;
    }

    metadata* block_metadata = (metadata*)blocks_start;

    size_t n_full_blocks = vtx_num_blocks - 1; // full blocks
    size_t bytes_to_copy = n_full_blocks*bytes_per_block;
    {
      // fetch original block
      size_t last_block_num = vtx_num_blocks-1;
      size_t orig_block_num = block_metadata[last_block_num].block_num;

      // get block size
      size_t block_start = orig_block_num * edges_per_block;
      size_t block_end =
          std::min(block_start + edges_per_block, static_cast<size_t>(vtx_original_degree));

      // #bytes for this block
      size_t last_block_size = block_end - block_start;
      size_t last_block_bytes = sizeof(metadata) + bitsets::get_bitset_block_size_in_bytes(last_block_size);

      bytes_to_copy += last_block_bytes;
    }
    assert(bytes_to_copy <= ((n_full_blocks+1)*bytes_per_block));

    // Is a blocked memcpy faster here?
    parallel_for(0, bytes_to_copy, [&] (size_t i) {
      tmp_space[i] = blocks_start[i];
    }, 512);

    // Tmp space for integers starts consecutively after tmp block space
    if (tmp) {
      tmp_ints = (uintE*)(tmp + bytes_to_copy);
    }

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
        size_t bytes_to_copy = bitsets::get_bitset_block_size_in_bytes(this_block_size);
        assert(bytes_to_copy <= bytes_per_block);

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

    if ((tmp == nullptr) && (old_vtx_num_blocks > kBlockAllocThreshold)) {
      pbbs::free_array(tmp_space);
      pbbs::free_array(tmp_ints);
    }
  }

  // P : (uintE, uintE, wgh) -> bool
  // only keep edges s.t. P(...) = true.
  template <class P>
  inline size_t pack_blocks(uintE vtx_id, P& p, uint8_t* tmp, bool parallel, const flags fl) {
    if (vtx_degree == 0) {
      return 0;
    }
    metadata* block_metadata = (metadata*)blocks_start;
    uintE int_stk[kBlockAllocThreshold];
    uintE* tmp_ints = (uintE*)int_stk;
    uintE old_vtx_num_blocks = vtx_num_blocks;
    if ((tmp == nullptr) && (vtx_num_blocks > kBlockAllocThreshold)) {
      tmp_ints = pbbs::new_array_no_init<uintE>(vtx_num_blocks);
    }
    if (tmp) {
      tmp_ints = (uintE*)tmp;
    }

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

              size_t block_size = block_end - block_start;
              size_t block_size_num_longs = (block_size+64-1)/64;

              // "select" on a block
              uint64_t* long_block_bits = (uint64_t*)block_bits;
              size_t cur_offset = block_start;
              size_t live_edges = 0;
              for (size_t idx = 0; idx < block_size_num_longs; idx++) {
                uint64_t cur_long = long_block_bits[idx];
                // size_t cnt = _mm_popcnt_u64(cur_long); // #bits set to one
                uint64_t long_to_write = cur_long;
                while (cur_long > 0) {
                  unsigned select_idx = _tzcnt_u64(cur_long);
                  auto& ee = e[cur_offset + select_idx];
                  if (!p(vtx_id, std::get<0>(ee), std::get<1>(ee))) {
                    long_to_write ^= (1UL << select_idx);
                  } else {
                    live_edges++;
                  }
                  assert((cur_long & (1UL << select_idx)) > 0);
                  cur_long = _blsr_u64(cur_long);
                }
                long_block_bits[idx] = long_to_write;
                cur_offset += 64; // next long
              }

//              size_t live_edges = 0;
//              for (size_t k = 0; k < (block_end - block_start); k++) {
//                if (bitsets::is_bit_set(block_bits, k)) {     // k-th edge is present in G
//                  auto& ee = e[block_start + k];  // fetch the k-th edge
//                  if (!p(vtx_id, std::get<0>(ee), std::get<1>(ee))) {
//                    // unset the k-th bit
//                    bitsets::flip_bit(block_bits, k);
//                  } else {
//                    // otherwise increment the count of live edges
//                    live_edges++;
//                  }
//                }
//              }

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

    if ((full_blocks * kFullBlockPackThreshold <= vtx_num_blocks) ||
        ((full_blocks < vtx_num_blocks) && (fl & compact_blocks))) {
      repack_blocks_par(tmp, parallel);
    }

    // Update offset values.
    auto ptr_seq = pbbs::indirect_value_seq<uintE>(
        vtx_num_blocks, [&](size_t i) { return &(block_metadata[i].offset); });
    uintE sum = pbbslib::scan_add_inplace(ptr_seq, pbbs::no_flag, tmp_ints);
    vtx_degree = sum;

    if ((tmp == nullptr) && (old_vtx_num_blocks > kBlockAllocThreshold)) {
      pbbs::free_array(tmp_ints);
    }

    // Update the degree in vtx_info.
    v_infos[vtx_id].vtx_degree = sum;
    assert(sum <= vtx_original_degree);
    return sum;  // return the new degree
  }

  std::tuple<uintE, W> ith_neighbor(size_t i) {
    metadata* block_metadata = (metadata*)blocks_start;
    auto offsets_imap = pbbslib::make_sequence<size_t>(vtx_num_blocks, [&] (size_t i) {
      return block_metadata[i].offset;
    });

    auto lte = [&](const size_t& l, const size_t& r) { return l <= r; };
    size_t block = pbbslib::binary_search(offsets_imap, i, lte);
    assert(block > 0);
    block = block-1;

    std::tuple<uintE, W> out;
    auto decode_f = [&] (const uintE& v, const W& wgh, const uintE& offset) {
      if (offset == i) {
        out = std::make_tuple(v, wgh);
        return false;
      }
      return true;
    };
    decode_block_cond(block, decode_f);
    return out;
  }

  struct iter {
    E* edges;
    uintE vtx_degree;
    uintE vtx_original_degree;
    uintE vtx_num_blocks;
    uint8_t* blocks_start;
    uint8_t* block_data_start;

    uintE cur_block; // block id
    uintE cur_block_degree; // #live edges
    uintE cur_block_start; // start offset (in edges)
    uintE cur_block_num_longs; // block size (in terms of #longs)

    uint64_t* cur_block_longs;
    uint64_t cur_long_value;
    uintE cur_long_idx;

    uintE last_ngh;
    uintE proc;
    uintE proc_cur_block;

    iter(E* edges,
        uintE vtx_degree,
        uintE vtx_original_degree,
        uintE vtx_num_blocks,
        uint8_t* blocks_start,
        uint8_t* block_data_start) :
        edges(edges),
        vtx_degree(vtx_degree),
        vtx_original_degree(vtx_original_degree),
        vtx_num_blocks(vtx_num_blocks),
        blocks_start(blocks_start),
        block_data_start(block_data_start),
        proc(0) {
      proc = 0;
      if (vtx_degree > 0) {
        next_nonempty_block(/* on_initialization = */ true);
        next(); // sets last_ngh
      }
    }

    // precondition: there is a subsequent non-empty block
    __attribute__((always_inline)) inline void next_nonempty_block(bool on_initialization=false) {
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
      assert(cur_block_degree > 0);

      // reset proc_cur_block
      proc_cur_block = 0;

      // set cur_block_start and cur_block_size
      uintE orig_block_num = block_metadata[cur_block].block_num;
      cur_block_start = orig_block_num * edges_per_block;
      uintE cur_block_end =
        std::min(cur_block_start + edges_per_block, vtx_original_degree);
      uintE cur_block_size = cur_block_end - cur_block_start;
      cur_block_num_longs = (cur_block_size+64-1)/64;

      // set cur_block_longs
      cur_block_longs = (uint64_t*)(block_data_start + bitset_bytes_per_block * cur_block);

      // initialize value to read, and the idx.
      cur_long_value = cur_block_longs[0];
      cur_long_idx = 0;
    }

    __attribute__((always_inline)) inline uintE degree() { return vtx_degree; }
    __attribute__((always_inline)) inline uintE cur() { return last_ngh; }

    // updates last_ngh
    __attribute__((always_inline)) inline uintE next() {
      while(cur_long_value == 0) {
        // done with block?
        cur_long_idx++;
        proc_cur_block += 64;
        if (cur_long_idx == cur_block_num_longs) {
          next_nonempty_block();
        }
        cur_long_value = cur_block_longs[cur_long_idx];
      }
      // cur_long_value > 0

      unsigned select_idx = _tzcnt_u64(cur_long_value); // #trailing zeros in cur_long
      cur_long_value = _blsr_u64(cur_long_value); // clears lowest set bit
      last_ngh = std::get<0>(edges[cur_block_start + proc_cur_block + select_idx]);
      proc++;
      return last_ngh;
    }

    __attribute__((always_inline)) inline bool has_next() {
      return proc < vtx_degree;
    }
  };

  auto get_iter() {
    return iter(get_edges(),
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

  vtx_info* v_infos;

  static constexpr uintE edges_per_block = PARALLEL_DEGREE;
  static constexpr uintE bytes_per_block = edges_per_block / 8 + sizeof(metadata);
  static constexpr uintE bitset_bytes_per_block =
      bytes_per_block - sizeof(metadata);

  E* e0;
#ifdef NVM
  E* e1;
#endif

  compressed_sym_bitset_manager(const uintE vtx_id, uint8_t* blocks,
                                uintE vtx_original_degree, vtx_info* v_infos,
#ifndef NVM
                     E* e0
#else
                     E* e0, E* e1
#endif
                                )
      : vtx_id(vtx_id),
      vtx_original_degree(vtx_original_degree),
      v_infos(v_infos),
#ifndef NVM
        e0(e0)
#else
        e0(e0), e1(e1)
#endif
  {
    auto& v_info = v_infos[vtx_id];
    vtx_degree = v_info.vtx_degree;
    vtx_num_blocks = v_info.vtx_num_blocks;

    blocks_start = blocks + v_info.vtx_block_offset;
    block_data_start = blocks_start + (vtx_num_blocks * sizeof(metadata));
  }

  __attribute__((always_inline)) inline E* get_edges() {
#ifndef NVM
    return e0;
#else
    if (numanode() == 0) {
      return e0;
    } else {
      return e1;
    }
#endif
  }

  __attribute__((always_inline)) inline uintE get_degree() {
    return vtx_degree;
  }

  inline void clear_vertex() {
    vtx_degree = 0;
    vtx_num_blocks = 0;
    v_infos[vtx_id].vtx_num_blocks = 0;
    v_infos[vtx_id].vtx_degree = 0;
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

    E* e = get_edges();

    // 1. decode the compressed block
    size_t k = 0;
    std::tuple<uintE, W> block_decode[edges_per_block];
    auto t_f = [&](const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
      block_decode[k++] = std::make_tuple(ngh, wgh);
    };
    bytepd_amortized::template decode_block<W>(t_f, e, vtx_id, vtx_original_degree,
                                               orig_block_num);

    uint64_t* long_block_bits = (uint64_t*)(block_data_start + bitset_bytes_per_block * block_id);
    uintE block_start = orig_block_num * edges_per_block;
    size_t block_end =
        std::min(block_start + edges_per_block, vtx_original_degree);
    size_t block_size = block_end - block_start;
    size_t block_size_num_longs = (block_size+64-1)/64;

    // "select" on a block
    size_t cur_offset = 0;
    for (size_t idx = 0; idx < block_size_num_longs; idx++) {
      uint64_t cur_long = long_block_bits[idx];
      while (cur_long > 0) {
        unsigned select_idx = _tzcnt_u64(cur_long); // #trailing zeros in cur_long
        auto& ee = block_decode[cur_offset + select_idx];
        f(std::get<0>(ee), std::get<1>(ee), offset++);
        assert((cur_long & (1UL << select_idx)) > 0);
        cur_long = _blsr_u64(cur_long); // clears lowest set bit
      }
      cur_offset += 64; // next long
    }

//    size_t k = 0;
//    uint8_t* block_bits = block_data_start + bitset_bytes_per_block * block_id;
//    auto t_f = [&](const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
//      if (bitsets::is_bit_set(block_bits,
//                              k++)) {  // check if the k-th bit is set
//        f(ngh, wgh, offset++);         // and apply f with the correct offset
//      }
//    };
//    bytepd_amortized::template decode_block<W>(t_f, e, vtx_id, vtx_original_degree,
//                                               orig_block_num);
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block_cond(uintE block_id,
                                                               F f) {
    metadata* block_metadata = (metadata*)blocks_start;
    uintE offset = block_metadata[block_id].offset;
    uintE orig_block_num = block_metadata[block_id].block_num;

    E* e = get_edges();

//    size_t k = 0;
//    std::tuple<uintE, W> block_decode[edges_per_block];
//    auto t_f = [&](const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
//      block_decode[k++] = std::make_tuple(ngh, wgh);
//    };
//    bytepd_amortized::template decode_block<W>(t_f, e, vtx_id, vtx_original_degree,
//                                               orig_block_num);
//
//    uint64_t* long_block_bits = (uint64_t*)(block_data_start + bitset_bytes_per_block * block_id);
//    uintE block_start = orig_block_num * edges_per_block;
//    size_t block_end =
//        std::min(block_start + edges_per_block, vtx_original_degree);
//    size_t block_size = block_end - block_start;
//    size_t block_size_num_longs = (block_size+64-1)/64;
//
//    // "select" on a block
//    size_t cur_offset = 0;
//    for (size_t idx = 0; idx < block_size_num_longs; idx++) {
//      uint64_t cur_long = long_block_bits[idx];
//      while (cur_long > 0) {
//        unsigned select_idx = _tzcnt_u64(cur_long); // #trailing zeros in cur_long
//        auto& ee = block_decode[cur_offset + select_idx];
//        if (!f(std::get<0>(ee), std::get<1>(ee), offset++)) {
//          return;
//        }
//        assert((cur_long & (1UL << select_idx)) > 0);
//        cur_long = _blsr_u64(cur_long); // clears lowest set bit
//      }
//      cur_offset += 64; // next long
//    }

    size_t k = 0;
    uint8_t* block_bits = block_data_start + bitset_bytes_per_block * block_id;
    auto t_f = [&](const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
     if (bitsets::is_bit_set(block_bits, // about 7% overhead
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
    assert(k == (block_end - block_start));
    );
  }

  /* Only called when the discrepency between full and total blocks is large.
   * Specifically, when #full_blocks*kFullBlockPackThreshold >= vtx_num_blocks
   * Note: Defers to pack(..) to finish updating the degree/scan info */
  inline void repack_blocks_par(uint8_t* tmp, bool parallel) {
    uint8_t stk[bytes_per_block * kBlockAllocThreshold];  // temporary space
    uintE int_stk[kBlockAllocThreshold];
    uint8_t* tmp_space = (uint8_t*)stk;
    uintE* tmp_ints = (uintE*)int_stk;
    size_t total_bytes = vtx_num_blocks * bytes_per_block;
    if ((tmp == nullptr) && (vtx_num_blocks > kBlockAllocThreshold)) {
      tmp_space = pbbs::new_array_no_init<uint8_t>(total_bytes);
      tmp_ints = pbbs::new_array_no_init<uintE>(vtx_num_blocks);
    }

    // caller supplies:
    // vtx_num_blocks*sizeof(uintE) +
    // vtx_num_blocks*bytes_per_block

    if (tmp) {
      tmp_space = tmp;
    }

    metadata* block_metadata = (metadata*)blocks_start;

    size_t n_full_blocks = vtx_num_blocks - 1; // full blocks
    size_t bytes_to_copy = n_full_blocks*bytes_per_block;
    {
      // fetch original block
      size_t last_block_num = vtx_num_blocks-1;
      size_t orig_block_num = block_metadata[last_block_num].block_num;

      // get block size
      size_t block_start = orig_block_num * edges_per_block;
      size_t block_end =
          std::min(block_start + edges_per_block, static_cast<size_t>(vtx_original_degree));

      // #bytes for this block
      size_t last_block_size = block_end - block_start;
      size_t last_block_bytes = sizeof(metadata) + bitsets::get_bitset_block_size_in_bytes(last_block_size);

      bytes_to_copy += last_block_bytes;
    }

    // Is a blocked memcpy faster here?
    parallel_for(0, bytes_to_copy, [&] (size_t i) {
      tmp_space[i] = blocks_start[i];
    }, 512);

    // Tmp space for integers starts consecutively after tmp block space
    if (tmp) {
      tmp_ints = (uintE*)(tmp + bytes_to_copy);
    }

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
        size_t block_bytes_to_copy = bitsets::get_bitset_block_size_in_bytes(this_block_size);

        // (b) copy bitset data
        for (size_t i = 0; i < block_bytes_to_copy; i++) {
          real_block_bits[i] = tmp_block_bits[i];
        }
      }
    });

    // 5. Update num_blocks info both locally and in v_infos
    uintE old_vtx_num_blocks = vtx_num_blocks;
    vtx_num_blocks = new_num_blocks;
    v_infos[vtx_id].vtx_num_blocks = new_num_blocks;

    if ((tmp == nullptr) && (old_vtx_num_blocks > kBlockAllocThreshold)) {
      pbbs::free_array(tmp_space);
      pbbs::free_array(tmp_ints);
    }
  }

  // P : (uintE, uintE, wgh) -> bool
  // only keep edges s.t. P(...) = true.
  template <class P>
  inline size_t pack_blocks(uintE vtx_id, P& p, uint8_t* tmp, bool parallel, const flags fl) {
    if (vtx_degree == 0) {
      return 0;
    }
    metadata* block_metadata = (metadata*)blocks_start;

    uintE int_stk[kBlockAllocThreshold];
    uintE* tmp_ints = (uintE*)int_stk;
    uintE old_vtx_num_blocks = vtx_num_blocks;
    if ((tmp == nullptr) && (vtx_num_blocks > kBlockAllocThreshold)) {
      tmp_ints = pbbs::new_array_no_init<uintE>(vtx_num_blocks);
    }
    if (tmp) {
      tmp_ints = (uintE*)tmp;
    }


    // 1. pack each block
    par_for(0, vtx_num_blocks, 1, [&](size_t block_id) {
        uintE orig_block_num = block_metadata[block_id].block_num;
        uint8_t* block_bits =
            block_data_start + bitset_bytes_per_block * block_id;

        E* e = get_edges();

//        // (i) decode the block

        size_t k = 0;
        std::tuple<uintE, W> block_decode[edges_per_block];
        auto t_f = [&](const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
          block_decode[k++] = std::make_tuple(ngh, wgh);
        };
        bytepd_amortized::template decode_block<W>(t_f, e, vtx_id, vtx_original_degree,
                                                   orig_block_num);

        uint64_t* long_block_bits = (uint64_t*)(block_data_start + bitset_bytes_per_block * block_id);
        uintE block_start = orig_block_num * edges_per_block;
        size_t block_end =
            std::min(block_start + edges_per_block, vtx_original_degree);
        size_t block_size = block_end - block_start;
        size_t block_size_num_longs = (block_size+64-1)/64;

        // "select" on a block
        size_t cur_offset = 0;
        size_t live_edges = 0;
        for (size_t idx = 0; idx < block_size_num_longs; idx++) {
          uint64_t cur_long = long_block_bits[idx];
          uint64_t long_to_write = cur_long;
          while (cur_long > 0) {
            unsigned select_idx = _tzcnt_u64(cur_long); // #trailing zeros in cur_long
            auto& ee = block_decode[cur_offset + select_idx];
            if (!p(vtx_id, std::get<0>(ee), std::get<1>(ee))) {
              long_to_write ^= (1UL << select_idx);
            } else {
              live_edges++;
            }
            assert((cur_long & (1UL << select_idx)) > 0);
            cur_long = _blsr_u64(cur_long); // clears lowest set bit
          }
          long_block_bits[idx] = long_to_write;
          cur_offset += 64; // next long
        }
        // Temporarily store #live_edges in offset positions.
        block_metadata[block_id].offset = live_edges;

//        size_t k = 0;
//        size_t live_edges = 0;
//        auto t_f = [&](const uintE& ngh, const W& wgh, const uintE& orig_edge_id) {
//          if (bitsets::is_bit_set(block_bits, k)) {  // check if the k-th bit is set
//            if (!p(vtx_id, ngh, wgh)) {
//              // unset the k-th bit
//              bitsets::flip_bit(block_bits, k);
//            } else {
//              // otherwise increment the count of live edges
//              live_edges++;
//            }
//          }
//          k += 1;
//        };
//        bytepd_amortized::template decode_block<W>(t_f, e, vtx_id, vtx_original_degree,
//                                                   orig_block_num);
//
//        debug(uintE block_start = orig_block_num * edges_per_block;
//        uintE block_end =
//            std::min(block_start + edges_per_block, vtx_original_degree);
//        assert(k == (block_end - block_start)););
//
//        // Temporarily store #live_edges in offset positions.
//        block_metadata[block_id].offset = live_edges;
    }, parallel);

    // 2. Reduce to get the #empty_blocks
    auto full_block_seq =
        pbbslib::make_sequence<size_t>(vtx_num_blocks, [&](size_t i) {
          return static_cast<size_t>(block_metadata[i].offset > 0);
        });
    size_t full_blocks = pbbslib::reduce_add(full_block_seq);


    if ((full_blocks * kFullBlockPackThreshold <= vtx_num_blocks) ||
        ((full_blocks < vtx_num_blocks) && (fl & compact_blocks))) {
      repack_blocks_par(tmp, parallel);
    }

    // Update offset values.
    auto ptr_seq = pbbs::indirect_value_seq<uintE>(
        vtx_num_blocks, [&](size_t i) { return &(block_metadata[i].offset); });
    uintE sum = pbbslib::scan_add_inplace(ptr_seq, pbbs::no_flag, (uintE*)tmp_ints);
    vtx_degree = sum;

    if ((tmp == nullptr) && (old_vtx_num_blocks > kBlockAllocThreshold)) {
      pbbs::free_array(tmp_ints);
    }

    // Update the degree in vtx_info.
    v_infos[vtx_id].vtx_degree = sum;
    return sum;  // return the new degree
  }

  std::tuple<uintE, W> ith_neighbor(size_t i) {
    metadata* block_metadata = (metadata*)blocks_start;
    auto offsets_imap = pbbslib::make_sequence<size_t>(vtx_num_blocks, [&] (size_t i) {
      return block_metadata[i].offset;
    });

    auto lte = [&](const size_t& l, const size_t& r) { return l <= r; };
    size_t block = pbbslib::binary_search(offsets_imap, i, lte);
    assert(block > 0);
    block = block-1;

    std::tuple<uintE, W> out;
    auto decode_f = [&] (const uintE& v, const W& wgh, const uintE& offset) {
      if (offset == i) {
        out = std::make_tuple(v, wgh);
        return false;
      }
      return true;
    };
    decode_block_cond(block, decode_f);
    return out;
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
    uintE cur_block_num_longs; // block size (in terms of #longs)

    uint64_t* cur_block_longs;
    uint64_t cur_long_value;
    uintE cur_long_idx;

    uintE block_decode[PARALLEL_DEGREE];

    uintE last_ngh;
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
        next_nonempty_block(/* on_initialization = */ true);
        next(); // sets last_ngh
      }
    }

    // precondition: there is a subsequent non-empty block
    __attribute__((always_inline)) inline void next_nonempty_block(bool on_initialization=false) {
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
      assert(cur_block_degree > 0);

      // reset proc_cur_block
      proc_cur_block = 0;

      // set cur_block_start and cur_block_size
      uintE orig_block_num = block_metadata[cur_block].block_num;
      cur_block_start = orig_block_num * edges_per_block;
      uintE cur_block_end =
        std::min(cur_block_start + edges_per_block, vtx_original_degree);
      uintE cur_block_size = cur_block_end - cur_block_start;
      cur_block_num_longs = (cur_block_size+64-1)/64;

      // set cur_block_longs
      cur_block_longs = (uint64_t*)(block_data_start + bitset_bytes_per_block * cur_block);

      // initialize value to read, and the idx.
      cur_long_value = cur_block_longs[0];
      cur_long_idx = 0;

      // decode current compressed block into block_decode
      size_t k = 0;
      auto map_f = [&] (const uintE& v, const W& wgh, size_t offset) {
        block_decode[k++] = v;
      };
      bytepd_amortized::template decode_block<W>(map_f, edges, vtx_id, vtx_original_degree, orig_block_num);
    }

    __attribute__((always_inline)) inline uintE degree() { return vtx_degree; }
    __attribute__((always_inline)) inline uintE cur() { return last_ngh; }

    // updates last_ngh
    __attribute__((always_inline)) inline uintE next() {
      while(cur_long_value == 0) {
        // done with block?
        cur_long_idx++;
        proc_cur_block += 64;
        if (cur_long_idx == cur_block_num_longs) {
          next_nonempty_block();
        }
        cur_long_value = cur_block_longs[cur_long_idx];
      }
      // cur_long_value > 0

      unsigned select_idx = _tzcnt_u64(cur_long_value); // #trailing zeros in cur_long
      cur_long_value = _blsr_u64(cur_long_value); // clears lowest set bit

      last_ngh = block_decode[proc_cur_block + select_idx];
      proc++;
      return last_ngh;
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
