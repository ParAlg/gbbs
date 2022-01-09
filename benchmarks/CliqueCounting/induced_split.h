#pragma once

#include "induced_hybrid.h"
#include "intersect.h"

/* TODO: describe what this file does */

namespace gbbs {
namespace induced_split {

template <class Graph, class F>
inline size_t CountCliques(Graph& DG, size_t k, F base_f, bool use_base = false,
                           bool label = true, long recursive_level = 0) {
  timer t;
  t.start();
  using W = typename Graph::weight_type;
  auto parallel_work = sequence<size_t>(DG.n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      return DG.get_vertex(v).out_degree();
    };
    parallel_for(0, DG.n, [&](size_t i) {
      auto monoid = parlay::addm<size_t>();
      parallel_work[i] = DG.get_vertex(i).out_neighbors().reduce(map_f, monoid);
    });
  }
  size_t total_work = parlay::scan_inplace(make_slice(parallel_work));

  size_t block_size = 50000;
  size_t n_blocks = total_work / block_size + 1;
  size_t work_per_block = total_work / n_blocks;
  n_blocks = (total_work / work_per_block) + 1;
  double tt = t.stop();
  std::cout << "##### Triangle scheduling: " << tt << std::endl;

  timer t2;
  t2.start();
  sequence<size_t> tots = sequence<size_t>::uninitialized(n_blocks);  // DG.n
  parallel_for(
      0, n_blocks,
      [&](size_t j) {
        size_t start = j * work_per_block;
        size_t end = (j + 1) * work_per_block;
        auto less_fn = std::less<size_t>();
        size_t start_ind = parlay::binary_search(parallel_work, start, less_fn);
        size_t end_ind = parlay::binary_search(parallel_work, end, less_fn);
        tots[j] = 0;
        for (size_t i = start_ind; i < end_ind; i++) {
          if (DG.get_vertex(i).out_degree() != 0) {
            HybridSpace_lw* induced = new HybridSpace_lw();
            induced->alloc(DG.get_vertex(i).out_degree(), k, DG.n, label,
                           use_base);
            induced->setup(DG, k, i);
            auto curr_counts = induced_hybrid::KCliqueDir_fast_hybrid_rec(
                DG, 1, k, induced, base_f, recursive_level);
            tots[j] += curr_counts;
            if (induced->use_base && curr_counts > 0) base_f(i, curr_counts);
            if (induced != nullptr) {
              delete induced;
            }
          }
        }
      },
      1, false);
  double tt2 = t2.stop();
  std::cout << "##### Actual counting: " << tt2 << std::endl;

  return parlay::reduce(tots);
}

}  // namespace induced_intersection
}  // namespace gbbs
