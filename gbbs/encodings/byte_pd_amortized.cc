#include "byte_pd_amortized.h"

namespace gbbs {
namespace bytepd_amortized {

long compressFirstEdge(uchar* start, long offset, long source, long target) {
  long diff = target - source;
  long preCompress = diff;
  int bytesUsed = 0;
  uchar firstByte = 0;
  uintE toCompress = std::abs(preCompress);
  firstByte = toCompress & 0x3f;  // 0011|1111
  if (preCompress < 0) {
    firstByte |= 0x40;
  }
  toCompress = toCompress >> 6;
  if (toCompress > 0) {
    firstByte |= 0x80;
  }
  start[offset] = firstByte;
  offset++;

  uchar curByte = toCompress & 0x7f;
  while ((curByte > 0) || (toCompress > 0)) {
    bytesUsed++;
    uchar toWrite = curByte;
    toCompress = toCompress >> 7;
    // Check to see if there's any bits left to represent
    curByte = toCompress & 0x7f;
    if (toCompress > 0) {
      toWrite |= 0x80;
    }
    start[offset] = toWrite;
    offset++;
  }
  return offset;
}

long compressEdge(uchar* start, long curOffset, uintE e) {
  uchar curByte = e & 0x7f;
  int bytesUsed = 0;
  while ((curByte > 0) || (e > 0)) {
    bytesUsed++;
    uchar toWrite = curByte;
    e = e >> 7;
    // Check to see if there's any bits left to represent
    curByte = e & 0x7f;
    if (e > 0) {
      toWrite |= 0x80;
    }
    start[curOffset] = toWrite;
    curOffset++;
  }
  return curOffset;
}

uintE get_num_blocks(uchar* edge_start, uintE degree) {
  if (degree == 0) {
    return 0;
  }
  uintE virtual_degree = *((uintE*)edge_start);
  size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
  return num_blocks;
}

uintE get_block_degree(uchar* edge_start, uintE degree, uintE block_num) {
  if (degree == 0) {
    return 0;
  }
  uintE virtual_degree = *((uintE*)edge_start);
  size_t num_blocks = 1 + (virtual_degree - 1) / PARALLEL_DEGREE;
  uintE* block_offsets = (uintE*)(edge_start + sizeof(uintE));

  auto block_ends = [&](size_t j) {
    uintE end = (j == (num_blocks - 1))
                    ? degree
                    : (*((uintE*)(edge_start + block_offsets[j])));
    return end;
  };
  uintE block_start = (block_num == 0) ? 0 : block_ends(block_num - 1);
  uintE block_end = block_ends(block_num);
  return block_end - block_start;  // TODO: check
}

}  // namespace bytepd_amortized
}  // namespace gbbs
