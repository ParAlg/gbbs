#pragma once

typedef uint32_t flags;
const flags no_output = 1;
const flags pack_edges = 2;
const flags sparse_blocked = 4;
const flags dense_forward = 8;
const flags dense_parallel = 16;
const flags remove_duplicates = 32;
const flags no_dense = 64;
const flags in_edges = 128;          // map over in edges instead of out edges
const flags fine_parallel = 1 << 8;  // split to a node-size of 1
const flags compact_blocks = 1 << 9;  // always call repack_blocks
inline bool should_output(const flags& fl) { return !(fl & no_output); }
