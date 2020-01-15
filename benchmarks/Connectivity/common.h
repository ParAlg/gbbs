#pragma once

#include "ligra/pbbslib/atomic_max_counter.h"
#include "ligra/pbbslib/atomic_sum_counter.h"

#include "connectit.h"

using parent = uintE;

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
  std::cout << "# Total num insertions = " << num_inserts << std::endl;
  std::cout << "# Total num queries = " << num_queries << std::endl;

  return result;
}
