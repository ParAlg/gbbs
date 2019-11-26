#pragma once

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


/********************************* vertex incremental ************************/

struct ResizableParents {

  parent** arrays;
  size_t initial_size;

  ResizableParents(size_t initial_size) : initial_size(initial_size) {
     /* double each time; 64 arrays is plenty */
    arrays = pbbs::new_array_no_init<parent*>(64);

    /* base array has size = initial_size */
    arrays[0] = pbbs::new_array_no_init<parent>(initial_size);
  }

  /* idx is 0-indexed */
  uintE get_array_idx(uintE idx) {
    uintE mapped_idx = idx + 1; /* convert to 1-indexed */
    mapped_idx += 1; /* shift so that values fall in [... 2^{k}] */
    uintE array_idx = pbbs::log2_up(mapped_idx);
    array_idx -= 1; /* make array index 0-indexed */
    assert(array_idx >= 0);
    return array_idx;
  }

  uintE get_array_position(uintE idx, uintE array_idx) {
    uintE array_slots = 1 << array_idx;
    uintE one_indexed_idx = idx + 1;
    /* clear top bit */

  }

};


