#pragma once

#include <ctype.h>
#include <stdlib.h>
#include <atomic>
#include <cstring>
#include <iostream>
#include <memory>
#include <type_traits>

#include "parallel.h"

#ifdef USEMALLOC
namespace pbbs {
inline void* my_alloc(size_t i) { return malloc(i); }
inline void my_free(void* p) { free(p); }
}  // namespace pbbs
#else
#include "alloc.h"
namespace pbbs {
inline void* my_alloc(size_t i) { return my_mem_pool.alloc(i); }
inline void my_free(void* p) { my_mem_pool.afree(p); }
}  // namespace pbbs
#endif

namespace pbbs {

struct empty {};

typedef uint32_t flags;
const flags no_flag = 0;
const flags fl_sequential = 1;
const flags fl_debug = 2;
const flags fl_time = 4;
const flags fl_conservative = 8;
const flags fl_inplace = 16;

template <typename T>
inline void assign_uninitialized(T& a, const T& b) {
  using TT = typename std::remove_volatile<T>::type;
  new (static_cast<void*>((TT*)std::addressof(a))) TT(b);
}

template <typename T>
inline void move_uninitialized(T& a, const T& b) {
  new (static_cast<void*>(std::addressof(a))) T(std::move(b));
}

// a 32-bit hash function
inline uint32_t hash32(uint32_t a) {
  a = (a + 0x7ed55d16) + (a << 12);
  a = (a ^ 0xc761c23c) ^ (a >> 19);
  a = (a + 0x165667b1) + (a << 5);
  a = (a + 0xd3a2646c) ^ (a << 9);
  a = (a + 0xfd7046c5) + (a << 3);
  a = (a ^ 0xb55a4f09) ^ (a >> 16);
  return a;
}

inline uint32_t hash32_2(uint32_t a) {
  uint32_t z = (a + 0x6D2B79F5UL);
  z = (z ^ (z >> 15)) * (z | 1UL);
  z ^= z + (z ^ (z >> 7)) * (z | 61UL);
  return z ^ (z >> 14);
}

inline uint32_t hash32_3(uint32_t a) {
  uint32_t z = a + 0x9e3779b9;
  z ^= z >> 15;  // 16 for murmur3
  z *= 0x85ebca6b;
  z ^= z >> 13;
  z *= 0xc2b2ae3d;  // 0xc2b2ae35 for murmur3
  return z ^= z >> 16;
}

// from numerical recipes
inline uint64_t hash64(uint64_t u) {
  uint64_t v = u * 3935559000370003845ul + 2691343689449507681ul;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >> 4;
  v *= 4768777513237032717ul;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v << 5;
  return v;
}

// a slightly cheaper, but possibly not as good version
// based on splitmix64
inline uint64_t hash64_2(uint64_t x) {
  x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
  x = x ^ (x >> 31);
  return x;
}

// Combines two hash values.
inline uint64_t hash_combine(uint64_t hash_value_1, uint64_t hash_value_2) {
  // This is the same as boost's 32-bit `hash_combine` implementation, but with
  // 2 ^ 64 / (golden ratio) chosen as an arbitrary 64-bit additive magic number
  // rather than 2 ^ 32 / (golden ratio).
  return hash_value_1 ^ (hash_value_2 + 0x9e3779b97f4a7c15 + (hash_value_1 << 6)
      + (hash_value_1 >> 2));
}

// Does not initialize the array
template <typename E>
E* new_array_no_init(size_t n, bool touch_pages = false) {  // true) {
  // pads in case user wants to allign with cache lines
  size_t line_size = 64;
  size_t bytes = ((n * sizeof(E)) / line_size + 1) * line_size;
  // E* r = (E*) aligned_alloc(line_size, bytes);
  E* r = (E*)my_alloc(bytes);
  if (r == NULL) {
    fprintf(stderr, "Cannot allocate space: %lu bytes", bytes);
    exit(1);
  }
  // parallel_for (size_t i = 0; i < bytes; i = i + (1 << 21)) ((bool*) r)[i] =
  // 0;
  return r;
}

// Initializes in parallel
template <typename E>
E* new_array(size_t n) {
  E* r = new_array_no_init<E>(n);
  if (!std::is_trivially_default_constructible<E>::value) {
    // if (!std::is_default_constructible<E>::value) {
    if (n > 2048) {
      auto f = [&](size_t i) { new ((void*)(r + i)) E; };
      parallel_for(0, n, f);
    } else
      for (size_t i = 0; i < n; i++) new ((void*)(r + i)) E;
  }
  return r;
}

inline void free_array(void* a) { my_free(a); }

// Destructs in parallel
template <typename E>
void delete_array(E* A, size_t n) {
  // C++14 -- suppored by gnu C++11
  if (!std::is_trivially_destructible<E>::value) {
    // if (!std::is_destructible<E>::value) {
    if (n > 2048) {
      auto f = [&](size_t i) { A[i].~E(); };
      parallel_for(0, n, f);
    } else
      for (size_t i = 0; i < n; i++) A[i].~E();
  }
  using NVE = typename std::remove_volatile<E>::type;
  NVE* AA = (NVE*)A;
  my_free(AA);
}

template <typename ET>
inline bool atomic_compare_and_swap(ET* a, ET oldval, ET newval) {
  if constexpr (sizeof(ET) == 1) {
    uint8_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint8_t*>(a), r_oval,
                                        r_nval);
  } else if constexpr (sizeof(ET) == 4) {
    uint32_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(a), r_oval,
                                        r_nval);
  } else if constexpr (sizeof(ET) == 8) {
    uint64_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(a), r_oval,
                                        r_nval);
  } else if constexpr (sizeof(ET) == 16) {
    __int128 r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap_16(reinterpret_cast<__int128*>(a),
                                           r_oval, r_nval);
  } else {
    std::cout << "Bad CAS Length" << sizeof(ET) << std::endl;
    exit(0);
  }
}

template <typename ET>
inline bool atomic_compare_and_swap(volatile ET* a, ET oldval, ET newval) {
  if (sizeof(ET) == 1) {
    uint8_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint8_t*>(a),
                                        r_oval, r_nval);
  } else if (sizeof(ET) == 4) {
    uint32_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint32_t*>(a),
                                        r_oval, r_nval);
  } else if (sizeof(ET) == 8) {
    uint64_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint64_t*>(a),
                                        r_oval, r_nval);
  } else if (sizeof(ET) == 16) {
    __int128 r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap_16(
        reinterpret_cast<volatile __int128*>(a), r_oval, r_nval);
  } else {
    std::cout << "Bad CAS Length" << sizeof(ET) << std::endl;
    exit(0);
  }
}

template <typename E, typename EV>
inline E fetch_and_add(E* a, EV b) {
  volatile E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!atomic_compare_and_swap(a, oldV, newV));
  return oldV;
}

template <typename E, typename EV>
inline void write_add(E* a, EV b) {
  // volatile E newV, oldV;
  E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!atomic_compare_and_swap(a, oldV, newV));
}

template <typename E, typename EV>
inline void write_add(std::atomic<E>* a, EV b) {
  // volatile E newV, oldV;
  E newV, oldV;
  do {
    oldV = a->load();
    newV = oldV + b;
  } while (!std::atomic_compare_exchange_strong(a, &oldV, newV));
}

template <typename ET, typename F>
inline bool write_min(ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = *a;
  while (less(b, c) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_min(volatile ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = *a;
  while (less(b, c) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_min(std::atomic<ET>* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = a->load();
  while (less(b, c) && !(r = std::atomic_compare_exchange_strong(a, &c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_max(ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = *a;
  while (less(c, b) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_max(volatile ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = *a;
  while (less(c, b) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_max(std::atomic<ET>* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = a->load();
  while (less(c, b) && !(r = std::atomic_compare_exchange_strong(a, &c, b)));
  return r;
}

// returns the log base 2 rounded up (works on ints or longs or unsigned
// versions)
template <class T>
size_t log2_up(T i) {
  size_t a = 0;
  T b = i - 1;
  while (b > 0) {
    b = b >> 1;
    a++;
  }
  return a;
}

size_t granularity(size_t n);

template <typename Lf, typename Rf>
static void par_do_if(bool do_parallel, Lf left, Rf right, bool cons = false) {
  if (do_parallel)
    par_do(left, right, cons);
  else {
    left();
    right();
  }
}

template <typename Lf, typename Mf, typename Rf>
inline void par_do3(Lf left, Mf mid, Rf right) {
  auto left_mid = [&]() { par_do(left, mid); };
  par_do(left_mid, right);
}

template <typename Lf, typename Mf, typename Rf>
static void par_do3_if(bool do_parallel, Lf left, Mf mid, Rf right) {
  if (do_parallel)
    par_do3(left, mid, right);
  else {
    left();
    mid();
    right();
  }
}

}  // namespace pbbs
