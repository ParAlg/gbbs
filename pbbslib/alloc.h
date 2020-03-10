#pragma once

#include <atomic>

#include "concurrent_stack.h"

#if defined(__APPLE__)
inline void* aligned_alloc(size_t a, size_t n) { return malloc(n); }
#else
#ifdef USEMALLOC
#include <malloc.h>
struct __mallopt {
  __mallopt() {
    mallopt(M_MMAP_MAX, 0);
    mallopt(M_TRIM_THRESHOLD, -1);
  }
};

__mallopt __mallopt_var;
#endif
#endif

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
