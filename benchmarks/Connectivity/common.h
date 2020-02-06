#pragma once

#include "ligra/pbbslib/atomic_max_counter.h"
#include "ligra/pbbslib/atomic_sum_counter.h"

#include "connectit.h"

using parent = uintE;

using incremental_update = std::tuple<uintE, uintE, UpdateType> ;

uintE largest_comp = UINT_E_MAX;

atomic_max_counter<uintE> max_pathlen;
atomic_sum_counter<size_t> total_pathlen;

void report_pathlen(uintE pathlen) {
#ifdef REPORT_PATH_LENGTHS
  max_pathlen.update_value(pathlen);
  total_pathlen.update_value(pathlen);
#endif
}

template <class Seq>
std::pair<pbbs::sequence<incremental_update>, size_t> reorder_updates(Seq& updates) {
  auto bool_seq = pbbslib::make_sequence<bool>(updates.size(), [&] (size_t i) {
    return std::get<2>(updates[i]) == query_type;
  });
  return pbbs::split_two(updates, bool_seq);
}

//auto annotate_updates(pbbs::sequence<std::tuple<uintE, uintE>>& updates, double insert_to_query) {
//  auto result = pbbs::sequence<std::tuple<uintE, uintE, UpdateType>>(updates.size());
//  auto rnd = pbbs::random();
//  parallel_for(0, updates.size(), [&] (size_t i) {
//    double our_rand = static_cast<double>(rnd.ith_rand(i) % UINT_E_MAX) / static_cast<double>(UINT_E_MAX);
//    uintE u, v;
//    std::tie(u,v) = updates[i];
//    if (our_rand < insert_to_query) {
//      result[i] = std::make_tuple(u, v, insertion_type);
//    } else {
//      result[i] = std::make_tuple(u, v, query_type);
//    }
//  });
//  return result;
//}

auto annotate_updates_insert(pbbs::sequence<std::tuple<uintE, uintE>>& updates, size_t n) {
  auto seq = pbbslib::make_sequence<std::tuple<uintE, uintE, UpdateType>>(n, [&] (size_t i) {
      auto& ith = updates[i];
      return std::make_tuple(std::get<0>(ith), std::get<1>(ith), insertion_type);
  });
  return seq;
}

auto annotate_updates(pbbs::sequence<std::tuple<uintE, uintE>>& updates, double insert_to_query, size_t n, bool permute = false) {
  size_t result_size = updates.size()/insert_to_query;
  if (insert_to_query > 1) {
    std::cout << "Error: 0 < insert_to_query < 1" << std::endl;
    abort();
  }
  auto result = pbbs::sequence<std::tuple<uintE, uintE, UpdateType>>(result_size);
  auto rnd = pbbs::random();
  parallel_for(0, updates.size(), [&] (size_t i) {
    uintE u, v;
    std::tie(u,v) = updates[i];
      result[i] = std::make_tuple(u, v, insertion_type);
  });
  parallel_for(updates.size(), result.size(), [&] (size_t i) {
    auto our_rnd = rnd.fork(i);
    auto u = our_rnd.ith_rand(0) % n;
    auto v = our_rnd.ith_rand(1) % n;
    result[i] = std::make_tuple(u, v, query_type);
  });
  if (permute) {
    return pbbs::random_shuffle(result);
  }
  return result;
}
