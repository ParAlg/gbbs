#pragma once

#include <atomic>

#include "concurrent_stack.h"

#if defined(__APPLE__)
namespace pbbs {
inline void* aligned_alloc(size_t a, size_t n) { return malloc(n); }
}  // namespace pbbs
#else
#ifdef USEMALLOC
struct __mallopt {
  __mallopt();
};

// This global variable invokes the constructor of `__mallopt` at program
// initialization. The constructor adjusts the behavior of memory-allocation
// functions like `malloc` for performance.
extern __mallopt __mallopt_var;
#endif
#endif

namespace pbbs {
struct mem_pool {
  concurrent_stack<void*>* buckets;
  static constexpr size_t header_size = 64;
  static constexpr size_t log_base = 20;
  static constexpr size_t num_buckets = 20;
  static constexpr size_t small_size_tag = 100;
  std::atomic<long> allocated{0};
  std::atomic<long> used{0};
  size_t mem_size;

  mem_pool();

  void* add_header(void* a);
  void* sub_header(void* a);
  void* alloc(size_t s);
  void afree(void* a);
  void clear();
};

static mem_pool my_mem_pool;
}  // namespace pbbs
