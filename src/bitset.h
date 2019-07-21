#pragma once

// block layout:
// 4 bytes for start offset
// bs/8 bytes for bits

static size_t bitset_bs_bytes(size_t block_size) {
  assert(block_size % 8 == 0);
  return (block_size/8) + sizeof(uintE);
}

static size_t bitset_init_all_blocks(uint8_t* finger,
    size_t num_blocks, size_t bs, size_t bs_in_bytes) {
  // uintE/block used at the front
  uintE* block_offsets = (uintE*)finger;
  size_t block_bytes = bs_in_bytes - sizeof(uintE); // bs_in_bytes - offset_bytes
  uint8_t* blocks_start = finger + (num_blocks*sizeof(uintE));
  parallel_for(0, num_blocks, [&] (size_t block_num) {
    block_offsets[block_num] = block_num*bs;
    uint8_t* start = blocks_start + block_bytes*block_num;
    for (size_t j=0; j<block_bytes; j++) {
      start[j] = std::numeric_limits<uint8_t>::max(); //full byte
    }
  }, 512);
}
