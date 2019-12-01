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

enum UpdateType {
  insertion_type,
  query_type
};

using incremental_update = std::tuple<uintE, uintE, UpdateType> ;

uintE largest_comp = UINT_E_MAX;


atomic_max_counter<uintE> max_pathlen;
atomic_max_counter<uintE> max_uf_tries;

atomic_sum_counter<size_t> total_pathlen;
atomic_sum_counter<size_t> total_uf_tries;

void report_pathlen(uintE pathlen) {
#ifdef REPORT_PATH_LENGTHS
  max_pathlen.update_value(pathlen);
  total_pathlen.update_value(pathlen);
#endif
}

void report_tries(uintE tries) {
#ifdef REPORT_MAX_TRIES
  max_uf_tries.update_value(tries);
  total_uf_tries.update_value(tries);
#endif
}

template <class Seq>
std::pair<pbbs::sequence<incremental_update>, size_t> reorder_updates(Seq& updates) {
  auto bool_seq = pbbslib::make_sequence<bool>(updates.size(), [&] (size_t i) {
    return std::get<2>(updates[i]) == query_type;
  });
  return pbbs::split_two(updates, bool_seq);
}

auto annotate_updates(pbbs::sequence<std::tuple<uintE, uintE>>& updates, double insert_to_query) {
  auto result = pbbs::sequence<std::tuple<uintE, uintE, UpdateType>>(updates.size());
  auto rnd = pbbs::random();
  parallel_for(0, updates.size(), [&] (size_t i) {
    double our_rand = static_cast<double>(rnd.ith_rand(i) % UINT_E_MAX) / static_cast<double>(UINT_E_MAX);
    uintE u, v;
    std::tie(u,v) = updates[i];
    if (our_rand < insert_to_query) {
      result[i] = std::make_tuple(u, v, insertion_type);
    } else {
      result[i] = std::make_tuple(u, v, query_type);
    }
  });
  size_t num_inserts = 0;
  size_t num_queries = 0;
  for (size_t i=0; i<result.size(); i++) {
    if (std::get<2>(result[i]) == insertion_type) {
      num_inserts++;
    } else {
      num_queries++;
    }
  }
  std::cout << "Total num insertions = " << num_inserts << std::endl;
  std::cout << "Total num queries = " << num_queries << std::endl;

  return result;
}
