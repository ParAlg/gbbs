#pragma once

#include "gbbs/helpers/atomic_max_counter.h"
#include "gbbs/helpers/atomic_sum_counter.h"

#include "connectit.h"

namespace gbbs {

using parent = uintE;
using incremental_update = std::tuple<uintE, uintE, UpdateType>;

constexpr uintE largest_comp = UINT_E_MAX;

void report_pathlen(uintE pathlen);

template <class Seq>
std::pair<sequence<incremental_update>, size_t> reorder_updates(Seq& updates) {
  auto bool_seq = parlay::delayed_seq<bool>(updates.size(), [&](size_t i) {
    return std::get<2>(updates[i]) == query_type;
  });
  return parlay::split_two(updates, bool_seq);
}

template <class W>
sequence<std::tuple<uintE, uintE, UpdateType>> annotate_updates_insert(
    sequence<std::tuple<uintE, uintE, W>>& updates, size_t n) {
  auto seq = parlay::delayed_seq<std::tuple<uintE, uintE, UpdateType>>(
      n, [&](size_t i) {
        auto& ith = updates[i];
        return std::make_tuple(std::get<0>(ith), std::get<1>(ith),
                               insertion_type);
      });
  return seq;
}

template <class W>
sequence<std::tuple<uintE, uintE, UpdateType>> annotate_updates(
    sequence<std::tuple<uintE, uintE, W>>& updates, double insert_to_query,
    size_t n, bool permute = false) {
  size_t result_size = updates.size() / insert_to_query;
  if (insert_to_query > 1) {
    std::cout << "Error: 0 < insert_to_query < 1" << std::endl;
    abort();
  }
  auto result = sequence<std::tuple<uintE, uintE, UpdateType>>(result_size);
  auto rnd = parlay::random();
  parallel_for(0, updates.size(), [&](size_t i) {
    uintE u, v;
    W w;
    std::tie(u, v, w) = updates[i];
    result[i] = std::make_tuple(u, v, insertion_type);
  });
  parallel_for(updates.size(), result.size(), [&](size_t i) {
    auto our_rnd = rnd.fork(i);
    auto u = our_rnd.ith_rand(0) % n;
    auto v = our_rnd.ith_rand(1) % n;
    result[i] = std::make_tuple(u, v, query_type);
  });
  if (permute) {
    return parlay::random_shuffle(result);
  }
  return result;
}

}  // namespace gbbs
