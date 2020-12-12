#include "utilities.h"

#include <math.h>

namespace pbbs {
  size_t granularity(size_t n) { return (n > 100) ? ceil(pow(n, 0.5)) : 100; }
}  // namespace pbbs
