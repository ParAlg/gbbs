#pragma once

/* Represents per-vertex information necessary to extract the neighbors of a
 * vertex out of Compressed Sparse Row (CSR) */
struct vertex_data {
  /* offset into the edges in CSR */
  size_t offset;
};
