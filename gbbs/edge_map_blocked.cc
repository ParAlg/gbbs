#include "gbbs/edge_map_blocked.h"

namespace gbbs {
void alloc_finish() {
  data_block_allocator::finish();
}
}  // namespace gbbs
