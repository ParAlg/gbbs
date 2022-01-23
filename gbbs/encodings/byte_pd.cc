#include "byte_pd.h"

namespace gbbs {
namespace bytepd {

/*
  Compresses the first edge, writing target-source and a sign bit.
*/
long compressFirstEdge(uchar* start, long offset, uintE source, uintE target) {
  intE preCompress = (intE)target - source;
  int bytesUsed = 0;
  uchar firstByte = 0;
  intE toCompress = std::abs(preCompress);
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
  return parlay::num_blocks(degree, PARALLEL_DEGREE);
}

uintE get_block_degree(uchar* edge_start, uintE degree, uintE block_num) {
  if (degree == 0) {
    return 0;
  }
  size_t num_blocks = get_num_blocks(edge_start, degree);
  assert(num_blocks > 0);
  if (block_num < (num_blocks - 1)) {
    return PARALLEL_DEGREE;
  }
  return degree % PARALLEL_DEGREE;
}

}  // namespace bytepd_amortized
}  // namespace gbbs
