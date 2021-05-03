#include "utilities.h"

#include <math.h>

namespace pbbs {
  long* global_cas_array = nullptr;
  size_t granularity(size_t n) { return (n > 100) ? ceil(pow(n, 0.5)) : 100; }
}  // namespace pbbs
