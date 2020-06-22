#include "alloc.h"

#include "memory_size.h"
#include "parallel.h"

#include <optional>

#if defined(__APPLE__)
#else
#ifdef USEMALLOC
#include <malloc.h>
__mallopt::__mallopt() {
  mallopt(M_MMAP_MAX, 0);
  mallopt(M_TRIM_THRESHOLD, -1);
}

__mallopt __mallopt_var = __mallopt();
#endif
#endif

namespace pbbs {
namespace {

// returns the log base 2 rounded up (works on ints or longs or unsigned
// versions)
//
// Though this is already defined in utilities.h, we copy it over here to avoid
// a circular dependency. (Consider cleaning this up.)
size_t log2_up(size_t i) {
  size_t a = 0;
  size_t b = i - 1;  // Beware of overflow...
  while (b > 0) {
    b = b >> 1;
    a++;
  }
  return a;
}

}  // namespace

mem_pool::mem_pool() : mem_size{getMemorySize()} {
  buckets = new concurrent_stack<void*>[num_buckets];
};

void* mem_pool::add_header(void* a) { return (void*)((char*)a + header_size); }

void* mem_pool::sub_header(void* a) { return (void*)((char*)a - header_size); }

void* mem_pool::alloc(size_t s) {
  size_t log_size = log2_up((size_t)s + header_size);
  if (log_size < 20) {
    void* a = (void*)aligned_alloc(header_size, s + header_size);
    *((size_t*)a) = small_size_tag;
    return add_header(a);
  }
  size_t bucket = log_size - log_base;
  std::optional<void*> r = buckets[bucket].pop();
  size_t n = ((size_t)1) << log_size;
  used += n;
  if (r.has_value()) {
    // if (n > 10000000) std::cout << "alloc: " << add_header(*r) << ", " << n <<
    // std::endl;
    return add_header(*r);
  } else {
    void* a = (void*)aligned_alloc(header_size, n);
    // if (n > 10000000) std::cout << "alloc: " << add_header(a) << ", " << n <<
    // std::endl;
    allocated += n;
    if (a == NULL) std::cout << "alloc failed" << std::endl;
    // a hack to make sure pages are touched in parallel
    // not the right choice if you want processor local allocations
    size_t stride = (1 << 21);  // 2 Mbytes in a huge page
    auto touch_f = [&](size_t i) { ((bool*)a)[i * stride] = 0; };
    parallel_for(0, n / stride, touch_f, 1);
    *((size_t*)a) = bucket;
    return add_header(a);
  }
}

void mem_pool::afree(void* a) {
  // std::cout << "free: " << a << std::endl;
  void* b = sub_header(a);
  size_t bucket = *((size_t*)b);
  if (bucket == small_size_tag)
    free(b);
  else if (bucket >= num_buckets) {
    std::cout << "corrupted header in free" << std::endl;
    abort();
  } else {
    size_t n = ((size_t)1) << (bucket + log_base);
    // if (n > 10000000) std::cout << "free: " << a << ", " << n << std::endl;
    used -= n;
    if (n > mem_size / 64) {  // fix to 64
      free(b);
      allocated -= n;
    } else {
      buckets[bucket].push(b);
    }
  }
}

void mem_pool::clear() {
  for (size_t i = 0; i < num_buckets; i++) {
    size_t n = ((size_t)1) << (i + log_base);
    std::optional<void*> r = buckets[i].pop();
    while (r.has_value()) {
      allocated -= n;
      free(*r);
      r = buckets[i].pop();
    }
  }
}
}  // namespace pbbs
