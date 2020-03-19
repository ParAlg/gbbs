#include "byte.h"

namespace byte {

long compressFirstEdge(uchar* start, long offset, long source,
                              long target) {
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

long compressEdge(uchar* start, long curOffset, uintE diff) {
  uchar curByte = diff & 0x7f;
  int bytesUsed = 0;
  while ((curByte > 0) || (diff > 0)) {
    bytesUsed++;
    uchar toWrite = curByte;
    diff = diff >> 7;
    // Check to see if there's any bits left to represent
    curByte = diff & 0x7f;
    if (diff > 0) {
      toWrite |= 0x80;
    }
    start[curOffset] = toWrite;
    curOffset++;
  }
  return curOffset;
}

}  // namespace byte
