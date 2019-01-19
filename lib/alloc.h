#pragma once
#include "concurrent_stack.h"
#include "utilities.h"

struct mem_pool {
  concurrent_stack<void*>* buckets;
  static constexpr size_t header_size = 64;
  static constexpr size_t log_base = 20;
  static constexpr size_t num_buckets = 20;
  static constexpr size_t small_size_tag = 100;
  
  mem_pool() {
    buckets = new concurrent_stack<void*>[num_buckets];
  };

  void* add_header(void* a) {
    return (void*) ((char*) a + header_size);}

  void* sub_header(void* a) {
    return (void*) ((char*) a - header_size);}

  void* alloc(size_t s) {
    size_t log_size = pbbs::log2_up((size_t) s + header_size);
    if (log_size < 20) {
      void* a = (void*) aligned_alloc(header_size, s + header_size);
      *((size_t*) a) = small_size_tag;
      return add_header(a);
    }
    size_t bucket = log_size - log_base;
    maybe<void*> r = buckets[bucket].pop();
    if (r) {
      return add_header(*r);
    }
    else {
      size_t n = ((size_t) 1) << log_size;
      void* a = (void*) aligned_alloc(header_size, n);
      if (a == NULL) std::cout << "alloc failed" << std::endl;
      // a hack to make sure pages are touched in parallel
      // not the right choice if you want processor local allocations
      size_t stride = (1 << 21); // 2 Mbytes in a huge page
      auto touch_f = [&] (size_t i) {((bool*) a)[i*stride] = 0;};
      par_for(0, n/stride, 1, touch_f);
      *((size_t*) a) = bucket;
      return add_header(a);
    }
  }

  void afree(void* a) {
    void* b = sub_header(a);
    size_t bucket = *((size_t*) b);
    if (bucket == small_size_tag) free(b);
    else if (bucket >= num_buckets) {
      std::cout << "corrupted header in free" << std::endl;
      abort();
    } else {
      buckets[bucket].push(b);
    }
  }

  void sizes() {
    for (size_t i = 0; i < num_buckets; i++) {
      std::cout << (((size_t) 1) << (i+log_base)) << " : "
		<< (buckets[i].size()) << std::endl;
    };
  }
} my_mem_pool;

void* my_alloc(size_t i) {
  return my_mem_pool.alloc(i);
}

void my_free(void* p) {
  my_mem_pool.afree(p);
}
