#pragma once

#include "gbbs/macros.h"

namespace gbbs {

/* Note that the block degree is computable by differencing two starts. */
struct vtx_info {
  /* vertex's (packed) degree. */
  uintE vtx_degree;
  /* number of blocks associated with v */
  uintE vtx_num_blocks;
  /* pointer into the block structure */
  size_t vtx_block_offset;

  vtx_info(uintE degree, uintE num_blocks, size_t block_offset)
      : vtx_degree(degree),
        vtx_num_blocks(num_blocks),
        vtx_block_offset(block_offset) {}
};

}  // namespace gbbs
