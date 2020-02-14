#pragma once

#include <atomic>

#include "concurrent_stack.h"

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

