#pragma once

#include "gbbs/pbbslib/atomic_max_counter.h"
#include "gbbs/pbbslib/atomic_sum_counter.h"

#include "connectit.h"

namespace gbbs {

using parent = uintE;
using incremental_update = std::tuple<uintE, uintE, UpdateType> ;

constexpr uintE largest_comp = UINT_E_MAX;

extern pbbslib::atomic_max_counter<uintE> max_pathlen;
extern pbbslib::atomic_sum_counter<size_t> total_pathlen;

void report_pathlen(uintE pathlen);

template <class Seq>
std::pair<pbbs::sequence<incremental_update>, size_t> reorder_updates(Seq& updates) {
  auto bool_seq = pbbslib::make_sequence<bool>(updates.size(), [&] (size_t i) {
    return std::get<2>(updates[i]) == query_type;
  });
  return pbbs::split_two(updates, bool_seq);
}

pbbs::sequence<std::tuple<uintE, uintE, UpdateType>>
annotate_updates_insert(pbbs::sequence<std::tuple<uintE, uintE>>& updates, size_t n);

pbbs::sequence<std::tuple<uintE, uintE, UpdateType>>
annotate_updates(pbbs::sequence<std::tuple<uintE, uintE>>& updates, double insert_to_query, size_t n, bool permute = false);

}  // namespace gbbs
