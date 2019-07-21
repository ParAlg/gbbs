#pragma once

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

    metadata(uintE block_num, uintE offset) :
      block_num(block_num), offset(offset) {}
  };

  // returns the number of bytes for a given block_size
  // note that the block_size must be a multiple of 8
  static size_t bitset_bs_bytes(size_t block_size) {
    assert(block_size % 8 == 0);
    return (block_size/8) + sizeof(metadata);
  }

  static void bitset_init_all_blocks(uint8_t* finger,
      size_t num_blocks, size_t bs, size_t bs_in_bytes) {
    // uintE/block used at the front
    metadata* block_metadata = (metadata*)finger;
    size_t bitset_bytes = bs_in_bytes - sizeof(metadata); // bs_in_bytes - metadata bytes
    uint8_t* bitset_data_start = finger + (num_blocks*sizeof(metadata));
    parallel_for(0, num_blocks, [&] (size_t block_num) {
      block_metadata[block_num] = metadata(block_num, block_num*bs);
      assert(block_metadata[block_num].block_num == block_num);
      assert(block_metadata[block_num].offset == block_num*bs);
      uint8_t* start = bitset_data_start + bitset_bytes*block_num;

      for (size_t j=0; j<bitset_bytes; j++) {
        start[j] = std::numeric_limits<uint8_t>::max(); //full byte
      }
    }, 512);
    for (size_t block_num=0; block_num<num_blocks; block_num++) {
      assert(block_metadata[block_num].block_num == block_num);
      assert(block_metadata[block_num].offset == block_num*bs);
    }
  }

  static uintE block_degree(uint8_t* finger, uintE block_num, uintE num_blocks, uintE degree) {
    metadata* block_metadata = (metadata*)finger;
    uintE offset = block_metadata[block_num].offset;
    uintE next_offset = (block_num == num_blocks - 1) ? degree : block_metadata[block_num+1].offset;
    return next_offset - offset;
  }

  static bool is_bit_set(uint8_t* finger, uintE k) {
    // 8 entries/byte, so block corresponding to k is k/8 = k >> 3;
    uintE byte_id = k >> 3;

    // idx within the byte is (k & byte-mask) = (k & 0x7) (value between 0--7)
    constexpr uintE byte_mask = 0x7;
    uint8_t offset_within_byte = k & byte_mask;

    uint8_t byte_to_test = finger[byte_id];
    return byte_to_test & (static_cast<uint8_t>(1) << offset_within_byte);
  }

} // namespace bitsets
