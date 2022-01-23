#pragma once

#include <x86intrin.h>
#include "gbbs/macros.h"

namespace gbbs {

// see bitmagic
// precondition: rank > 0
unsigned select64_tz(uint64_t w, unsigned rank) {
  //    BM_ASSERT(w);
  //    BM_ASSERT(rank);
  //    BM_ASSERT(rank <= _mm_popcnt_u64(w));

  do {
    if ((--rank) == 0) break;
    w = _blsr_u64(w);  // w &= w - 1;
  } while (1);
  auto t = _blsi_u64(w);  // w & -w;
  unsigned count = unsigned(_tzcnt_u64(t));
  return count;
}

namespace bitsets {
// block layout:
// 4 bytes for block_num
// 4 bytes for offset
// bs/8 bytes for bits
//
// the (blck_num, offset) is the metadata. all block metadata are stored
// contiguously at the start of the array, and the block bits are stored
// contiguously.

struct metadata {
  // block in the original edge-set that this block corresponds to
  uintE block_num;
  // the starting idx of edges in this block
  uintE offset;

  metadata(uintE block_num, uintE offset)
      : block_num(block_num), offset(offset) {}
};

// returns the number of bytes for a given block_size
// note that the block_size must be a multiple of 8
// TODO: should we be rounding to the nearest multiple of 4?
static uintE bitset_bs_bytes(uintE block_size) {
  assert(block_size % 8 == 0);
  return (block_size / 8) + sizeof(metadata);
}

// #bitset bytes for a given block size. rounds up to the nearest multiple of 8
// bytes.
static inline size_t get_bitset_block_size_in_bytes(size_t block_size) {
  size_t block_bytes = (block_size + 8 - 1) / 8;  // ceil(_/8)
  return ((block_bytes + 8 - 1) & -8);            // round up to mult 8
}

// Returns the number of bytes needed for the bitset structure for a given
// degree and block size.
static uintE bytes_for_degree_and_bs(uintE degree, uintE bs,
                                     uintE bs_in_bytes) {
  if (degree == 0) {
    return 0;
  }
  uintE num_blocks = parlay::num_blocks(degree, bs);
  uintE full_blocks = num_blocks - 1;
  uintE rem = degree % bs;
  if (rem == 0) {
    full_blocks += 1;
    return full_blocks * bs_in_bytes;
  }

  // uintE rem_bytes = (rem + 8 - 1)/8; // ceil(rem/8)

  // rounds rem to nearest multiple of 8 bytes
  uintE last_block_bytes =
      sizeof(metadata) + get_bitset_block_size_in_bytes(rem);
  //    ((rem_bytes + 8 - 1) & -8);

  uintE full_block_bytes = full_blocks * bs_in_bytes;

  uintE total_bytes = full_block_bytes + last_block_bytes;
  return total_bytes;
}

__attribute__((always_inline)) static inline bool is_bit_set(uint8_t* finger,
                                                             uintE k) {
  // 8 entries/byte, so block corresponding to k is k/8 = k >> 3;
  uintE byte_id = k >> 3;

  // idx within the byte is (k & byte-mask) = (k & 0x7) (value between 0--7)
  constexpr uintE byte_mask = 0x7;
  uint8_t offset_within_byte = k & byte_mask;

  uint8_t byte_to_test = finger[byte_id];
  return byte_to_test & (static_cast<uint8_t>(1) << offset_within_byte);
}

__attribute__((always_inline)) static inline void flip_bit(uint8_t* finger,
                                                           uintE k) {
  // 8 entries/byte, so block corresponding to k is k/8 = k >> 3;
  uintE byte_id = k >> 3;

  // idx within the byte is (k & byte-mask) = (k & 0x7) (value between 0--7)
  constexpr uintE byte_mask = 0x7;
  uint8_t offset_within_byte = k & byte_mask;

  uint8_t byte_to_test = finger[byte_id];
  finger[byte_id] =
      byte_to_test ^ (static_cast<uint8_t>(1) << offset_within_byte);
}

static void bitset_init_blocks(uint8_t* finger, uintE degree, size_t num_blocks,
                               size_t bs, size_t bs_in_bytes,
                               size_t vtx_bytes) {
  if (degree != 0) {
    metadata* block_metadata = (metadata*)finger;
    parallel_for(0, num_blocks,
                 [&](size_t block_num) {
                   block_metadata[block_num] =
                       metadata(block_num, block_num * bs);
                   assert(block_metadata[block_num].block_num == block_num);
                   assert(block_metadata[block_num].offset == block_num * bs);
                 },
                 512);

    // Remaining bytes are data-bytes, all initially set to one.
    size_t data_bytes = vtx_bytes - (num_blocks * sizeof(metadata));
    uint8_t* bitset_data_start = finger + (num_blocks * sizeof(metadata));
    parallel_for(0, data_bytes,
                 [&](size_t i) {
                   bitset_data_start[i] = std::numeric_limits<uint8_t>::max();
                 },
                 512);

    uintE rem = degree % bs;
    if (rem > 0) {
      // fix bits for last block
      size_t last_block_num = num_blocks - 1;
      size_t last_block_start = last_block_num * bs;
      size_t last_block_size = degree - last_block_start;  // #set bits

      size_t last_block_bytes = get_bitset_block_size_in_bytes(last_block_size);

      size_t bs_data_bytes = bs_in_bytes - sizeof(metadata);
      //  size_t last_block_bytes = data_bytes - (num_blocks-1)*bs_data_bytes;
      size_t last_block_physical_size = last_block_bytes * 8;
      uint8_t* last_block_bits =
          bitset_data_start + bs_data_bytes * last_block_num;
      for (uintE k = last_block_size; k < last_block_physical_size; k++) {
        assert(is_bit_set(last_block_bits, k));
        flip_bit(last_block_bits, k);
      }
    }
  }
}

static inline uintE block_degree(uint8_t* finger, uintE block_num,
                                 uintE num_blocks, uintE degree) {
  metadata* block_metadata = (metadata*)finger;
  uintE offset = block_metadata[block_num].offset;
  uintE next_offset = (block_num == num_blocks - 1)
                          ? degree
                          : block_metadata[block_num + 1].offset;
  return next_offset - offset;
}
}  // namespace bitsets

}  // namespace gbbs
