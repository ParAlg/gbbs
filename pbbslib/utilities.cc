#include "utilities.h"

#include <math.h>

#ifdef USEMALLOC
void* my_alloc(size_t i) {return malloc(i);}
void my_free(void* p) {free(p);}
#else
#include "alloc.h"
void* my_alloc(size_t i) {return my_mem_pool.alloc(i);}
void my_free(void* p) {my_mem_pool.afree(p);}
#endif

namespace pbbs {

uint32_t hash32(uint32_t a) {
  a = (a+0x7ed55d16) + (a<<12);
  a = (a^0xc761c23c) ^ (a>>19);
  a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);
  a = (a+0xfd7046c5) + (a<<3);
  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}

uint32_t hash32_2(uint32_t a) {
  uint32_t z = (a + 0x6D2B79F5UL);
  z = (z ^ (z >> 15)) * (z | 1UL);
  z ^= z + (z ^ (z >> 7)) * (z | 61UL);
  return z ^ (z >> 14);
}

uint32_t hash32_3(uint32_t a) {
    uint32_t z = a + 0x9e3779b9;
    z ^= z >> 15; // 16 for murmur3
    z *= 0x85ebca6b;
    z ^= z >> 13;
    z *= 0xc2b2ae3d; // 0xc2b2ae35 for murmur3
    return z ^= z >> 16;
}

// from numerical recipes
uint64_t hash64(uint64_t u) {
  uint64_t v = u * 3935559000370003845ul + 2691343689449507681ul;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >>  4;
  v *= 4768777513237032717ul;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v <<  5;
  return v;
}

uint64_t hash64_2(uint64_t x) {
  x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
  x = x ^ (x >> 31);
  return x;
}

void free_array(void* a) {
  my_free(a);
}

size_t granularity(size_t n) {
  return (n > 100) ? ceil(pow(n,0.5)) : 100;
}

void assert_str(int cond, const std::string& s) {
  if (!cond) {std::cout << "PBBS assert error: " << s << std::endl;}
}

}  // namespace pbbs
