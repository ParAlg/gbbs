#pragma once

#include "benchmarks/Connectivity/UnionFind/atomic_max_counter.h"
#include "benchmarks/Connectivity/UnionFind/atomic_sum_counter.h"
//struct parent {
//  volatile uintE parent;
//};

//#define parent volatile uintE
#define parent uintE

/**************************** Framework options ****************************/
enum SamplingOption {
  sample_kout, sample_bfs, sample_ldd, no_sampling
};

/* Union-Find options */
enum FindOption {
  find_compress, find_naive, find_split, find_halve, find_atomic_split, find_atomic_halve
};
enum UniteOption {
  unite, unite_early, unite_nd, unite_rem_lock, unite_rem_cas
};

/* RemCAS-specific options */
enum SpliceOption {
  split_atomic_one, halve_atomic_one, splice_simple, splice_atomic
};

/* Jayanti-specific options */
enum JayantiFindOption {
  find_twotrysplit, find_simple
};

/* LiuTarjan-specific options */
enum LiuTarjanConnectOption {
  simple_connect, parent_connect, extended_connect
};
enum LiuTarjanUpdateOption {
  simple_update, root_update
};
enum LiuTarjanShortcutOption {
  shortcut, full_shortcut
};
enum LiuTarjanAlterOption {
  alter, no_alter
};

/* Finish Algorithm Types */
enum AlgorithmType {
  union_find_type, liu_tarjan_type, shiloach_vishkin_type, label_prop_type
};

uintE largest_comp = UINT_E_MAX;


atomic_max_counter<uintE> max_pathlen;
atomic_max_counter<uintE> max_uf_tries;

atomic_sum_counter<size_t> total_pathlen;
atomic_sum_counter<size_t> total_uf_tries;

