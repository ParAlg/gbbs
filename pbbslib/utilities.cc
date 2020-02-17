#include "utilities.h"

#include <math.h>

#ifdef USEMALLOC
void* my_alloc(size_t i) { return malloc(i); }
void my_free(void* p) { free(p); }
#else
#include "alloc.h"
void* my_alloc(size_t i) { return my_mem_pool.alloc(i); }
void my_free(void* p) { my_mem_pool.afree(p); }
#endif

namespace pbbs {

void free_array(void* a) { my_free(a); }

size_t granularity(size_t n) { return (n > 100) ? ceil(pow(n, 0.5)) : 100; }

void assert_str(int cond, const std::string& s) {
  if (!cond) {
    std::cout << "PBBS assert error: " << s << std::endl;
  }
}

}  // namespace pbbs
