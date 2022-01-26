#include "benchmarks/Connectivity/common.h"

namespace gbbs {

void report_pathlen(uintE pathlen) {
}

sequence<std::tuple<uintE, uintE, UpdateType>> annotate_updates_insert(
    sequence<std::tuple<uintE, uintE>>& updates, size_t n) {
  auto seq =
      parlay::sequence<std::tuple<uintE, uintE, UpdateType>>::from_function(
          n, [&](size_t i) {
            auto& ith = updates[i];
            return std::make_tuple(std::get<0>(ith), std::get<1>(ith),
                                   insertion_type);
          });
  return seq;
}

sequence<std::tuple<uintE, uintE, UpdateType>> annotate_updates(
    sequence<std::tuple<uintE, uintE>>& updates, double insert_to_query,
    size_t n, bool permute) {
  size_t result_size = updates.size() / insert_to_query;
  if (insert_to_query > 1) {
    std::cout << "Error: 0 < insert_to_query < 1" << std::endl;
    abort();
  }
  auto result = sequence<std::tuple<uintE, uintE, UpdateType>>(result_size);
  auto rnd = parlay::random();
  parallel_for(0, updates.size(), [&](size_t i) {
    uintE u, v;
    std::tie(u, v) = updates[i];
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
