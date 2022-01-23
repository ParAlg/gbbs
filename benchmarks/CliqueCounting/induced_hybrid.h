#pragma once

/* TODO: describe what this file does */

#include "intersect.h"

namespace gbbs {
namespace induced_hybrid {

template <class Graph>
size_t get_max_deg(Graph& DG) {
  size_t max_deg = 0;
  parallel_for(0, DG.n, [&](size_t i) {
    size_t deg = DG.get_vertex(i).out_degree();
    gbbs::write_min(&max_deg, deg, std::greater<size_t>());
  });
  return max_deg;
}

template <class Graph, class F>
inline size_t KCliqueDir_fast_hybrid_rec(Graph& DG, size_t k_idx, size_t k,
                                         HybridSpace_lw* induced, F base_f,
                                         size_t recursive_level = 0) {
  size_t num_induced = induced->num_induced[k_idx - 1];
  if (num_induced == 0) return 0;
  uintE* prev_induced = induced->induced + induced->nn * (k_idx - 1);

  for (size_t i = 0; i < num_induced; i++) {
    induced->labels[prev_induced[i]] = k_idx;
  }

  if (k_idx + 1 == k) {
    size_t counts = 0;
    for (size_t i = 0; i < num_induced; i++) {
      uintE vtx = prev_induced[i];
      //  get neighbors of vtx
      uintE* intersect = induced->induced_edges + vtx * induced->nn;
      size_t tmp_counts = 0;
      for (size_t j = 0; j < induced->induced_degs[vtx]; j++) {
        if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
          tmp_counts++;
          if (induced->use_base) base_f(induced->relabel[intersect[j]], 1);
        }
      }
      if (induced->use_base && tmp_counts > 0)
        base_f(induced->relabel[vtx], tmp_counts);
      counts += tmp_counts;
    }
    for (size_t i = 0; i < num_induced; i++) {
      induced->labels[prev_induced[i]] = k_idx - 1;
    }
    return counts;
  }
  size_t total_ct = 0;
  if (recursive_level < k_idx || num_induced < 2) {
    for (size_t i = 0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      uintE* intersect = induced->induced_edges + vtx * induced->nn;
      uintE* out = induced->induced + induced->nn * k_idx;
      uintE count = 0;
      for (size_t j = 0; j < induced->induced_degs[vtx]; j++) {
        if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
          out[count] = intersect[j];
          count++;
        }
      }
      induced->num_induced[k_idx] = count;
      if (induced->num_induced[k_idx] > k - k_idx - 1) {
        auto curr_counts = KCliqueDir_fast_hybrid_rec(DG, k_idx + 1, k, induced,
                                                      base_f, recursive_level);
        total_ct += curr_counts;
        if (induced->use_base && curr_counts > 0)
          base_f(induced->relabel[vtx], curr_counts);
      }
    }
  } else {
    sequence<size_t> tots = sequence<size_t>::uninitialized(num_induced);
    auto init_induced = [&](HybridSpace_lw* induced2) {
      induced2->alloc(induced->nn, k - k_idx, DG.n, false, induced->use_base,
                      false);
    };
    auto finish_induced = [&](HybridSpace_lw* induced2) {
      if (induced2 != nullptr) {
        delete induced2;
      }
    };
    parallel_for_alloc<HybridSpace_lw>(
        init_induced, finish_induced, 0, num_induced,
        [&](size_t i, HybridSpace_lw* induced2) {
          induced2->copy(induced);
          uintE vtx = prev_induced[i];
          uintE* intersect = induced->induced_edges + vtx * induced->nn;
          uintE* out = induced2->induced;
          uintE count = 0;
          for (size_t j = 0; j < induced->induced_degs[vtx]; j++) {
            if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
              out[count] = intersect[j];
              count++;
            }
          }
          induced2->num_induced[0] = count;
          if (count > k - k_idx - 1) {
            tots[i] =
                KCliqueDir_fast_hybrid_rec(DG, 1, k - k_idx, induced2, base_f);
            if (induced->use_base && tots[i] > 0)
              base_f(induced->relabel[vtx], tots[i]);
          } else
            tots[i] = 0;
        });
    total_ct += parlay::reduce(tots);
  }
  // for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] =
  // k_idx - 1; }
  for (size_t i = 0; i < num_induced; i++) {
    induced->labels[prev_induced[i]] = k_idx - 1;
  }
  return total_ct;
}

// TODO: This function does not support recursive levels > 0
// TODO: This function is very similar to the above (can probably be
// modularized)
template <class Graph, class F>
inline size_t KCliqueDir_fast_hybrid_rec_enum(Graph& DG, size_t k_idx, size_t k,
                                              HybridSpace_lw* induced, F base_f,
                                              sequence<uintE>& base) {
  size_t num_induced = induced->num_induced[k_idx - 1];
  if (num_induced == 0) return 0;
  uintE* prev_induced = induced->induced + induced->nn * (k_idx - 1);

  for (size_t i = 0; i < num_induced; i++) {
    induced->labels[prev_induced[i]] = k_idx;
  }

  if (k_idx + 1 == k) {
    size_t counts = 0;
    for (size_t i = 0; i < num_induced; i++) {
      uintE vtx = prev_induced[i];
      //  get neighbors of vtx
      uintE* intersect = induced->induced_edges + vtx * induced->nn;
      size_t tmp_counts = 0;
      base[k_idx] = induced->relabel[vtx];
      for (size_t j = 0; j < induced->induced_degs[vtx]; j++) {
        if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
          base[k] = induced->relabel[intersect[j]];
          tmp_counts++;
          base_f(base);
        }
      }
      counts += tmp_counts;
    }
    for (size_t i = 0; i < num_induced; i++) {
      induced->labels[prev_induced[i]] = k_idx - 1;
    }
    return counts;
  }

  size_t total_ct = 0;
  for (size_t i = 0; i < num_induced; ++i) {
    uintE vtx = prev_induced[i];
    uintE* intersect = induced->induced_edges + vtx * induced->nn;
    uintE* out = induced->induced + induced->nn * k_idx;
    uintE count = 0;
    for (size_t j = 0; j < induced->induced_degs[vtx]; j++) {
      if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
        out[count] = intersect[j];
        count++;
      }
    }
    induced->num_induced[k_idx] = count;
    if (induced->num_induced[k_idx] > k - k_idx - 1) {
      base[k_idx] = induced->relabel[vtx];
      auto curr_counts = KCliqueDir_fast_hybrid_rec_enum(DG, k_idx + 1, k,
                                                         induced, base_f, base);
      total_ct += curr_counts;
    }
  }
  for (size_t i = 0; i < num_induced; i++) {
    induced->labels[prev_induced[i]] = k_idx - 1;
  }
  return total_ct;
}

template <class Graph, class F>
inline size_t CountCliques(Graph& DG, size_t k, F base_f, bool use_base = false,
                           bool label = true, long recursive_level = 0) {
  timer t2;
  t2.start();

  if (recursive_level != 1) {
    sequence<size_t> tots = sequence<size_t>::uninitialized(DG.n);
    size_t max_deg = get_max_deg(DG);
    auto init_induced = [&](HybridSpace_lw* induced) {
      induced->alloc(max_deg, k, DG.n, label, use_base);
    };
    auto finish_induced = [&](HybridSpace_lw* induced) {
      if (induced != nullptr) {
        delete induced;
      }
    };
    parallel_for_alloc<HybridSpace_lw>(
        init_induced, finish_induced, 0, DG.n,
        [&](size_t i, HybridSpace_lw* induced) {
          if (DG.get_vertex(i).out_degree() != 0) {
            induced->setup(DG, k, i);
            tots[i] = KCliqueDir_fast_hybrid_rec(DG, 1, k, induced, base_f,
                                                 recursive_level);
            if (induced->use_base && tots[i] > 0) base_f(i, tots[i]);
          } else
            tots[i] = 0;
        },
        1, false);
    double tt2 = t2.stop();
    std::cout << "##### Actual counting: " << tt2 << std::endl;
    return parlay::reduce(tots);
  }

  size_t max_deg = get_max_deg(DG);
  sequence<size_t> degs = sequence<size_t>::uninitialized(DG.n + 1);
  parallel_for(0, DG.n,
               [&](size_t i) { degs[i] = DG.get_vertex(i).out_degree(); });
  degs[DG.n] = 0;
  size_t num_edges = parlay::scan_inplace(make_slice(degs));
  num_edges = degs[DG.n];
  sequence<size_t> tots = sequence<size_t>::uninitialized(num_edges);

  auto init_induced = [&](HybridSpace_lw* induced) {
    induced->alloc(max_deg, k, DG.n, label, use_base);
  };
  auto finish_induced = [&](HybridSpace_lw* induced) {
    if (induced != nullptr) {
      delete induced;
    }
  };
  parallel_for_alloc<HybridSpace_lw>(
      init_induced, finish_induced, 0, num_edges,
      [&](size_t j, HybridSpace_lw* induced) {
        // to find i and ngh, binary search for j in degs; the index - 1 is i
        // then, DG.get_vertex(i).getOutNeighbor(j - degs[index-1]) is ngh
        auto less_fn = [&](size_t a, size_t b) { return a <= b; };
        size_t idx = parlay::binary_search(degs, j, less_fn);
        auto i = idx - 1;
        auto ngh =
            DG.get_vertex(i).out_neighbors().get_neighbor(j - degs[idx - 1]);
        induced->setup_edge(DG, k, i, ngh);
        tots[j] = KCliqueDir_fast_hybrid_rec(DG, 1, k - 1, induced, base_f, 0);
        if (use_base && tots[j] > 0) {
          base_f(ngh, tots[j]);
          base_f(i, tots[j]);
        }
      },
      1, false);
  double tt2 = t2.stop();
  std::cout << "##### Actual counting: " << tt2 << std::endl;

  return parlay::reduce(tots);
}

template <class Graph, class F>
inline size_t CountCliquesEnum(Graph& DG, size_t k, F base_f,
                               bool label = true) {
  timer t2;
  t2.start();

  sequence<size_t> tots = sequence<size_t>::uninitialized(DG.n);
  size_t max_deg = get_max_deg(DG);
  auto init_induced = [&](HybridSpace_lw* induced) {
    induced->alloc(max_deg, k, DG.n, label, true);
  };
  auto finish_induced = [&](HybridSpace_lw* induced) {
    if (induced != nullptr) {
      delete induced;
    }
  };
  parallel_for_alloc<HybridSpace_lw>(
      init_induced, finish_induced, 0, DG.n,
      [&](size_t i, HybridSpace_lw* induced) {
        if (DG.get_vertex(i).out_degree() != 0) {
          induced->setup(DG, k, i);
          auto base = sequence<uintE>(k + 1);
          base[0] = i;
          tots[i] =
              KCliqueDir_fast_hybrid_rec_enum(DG, 1, k, induced, base_f, base);
        } else
          tots[i] = 0;
      },
      1, false);
  double tt2 = t2.stop();
  std::cout << "##### Actual counting: " << tt2 << std::endl;
  return parlay::reduce(tots);
}

}  // namespace induced_neighborhood
}  // namespace gbbs
