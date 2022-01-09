#pragma once

/* TODO: describe what this file does */

#include "intersect.h"

namespace gbbs {
namespace induced_neighborhood {

template <class Graph>
size_t get_max_deg(Graph& DG) {
  size_t max_deg = 0;
  parallel_for(0, DG.n, [&](size_t i) {
    size_t deg = DG.get_vertex(i).out_degree();
    gbbs::write_min(&max_deg, deg, std::greater<size_t>());
  });
  return max_deg;
}

template <class Graph>
inline size_t KCliqueDir_fast_orig_rec(Graph& DG, size_t k_idx, size_t k,
                                       FullSpace_orig_lw* induced) {
  size_t num_induced = induced->num_induced[k_idx - 1];
  uintE* prev_induced = induced->induced + induced->nn * (k_idx - 1);
  uintE* prev_induced_degs = induced->induced_degs + induced->nn * (k_idx - 1);
  if (num_induced == 0) return 0;

  if (k_idx + 1 == k) return induced->num_edges[k_idx - 1];

  size_t total_ct = 0;
  for (size_t i = 0; i < num_induced; ++i) {
    uintE idx = prev_induced[i];
    uintE* new_induced = induced->induced + induced->nn * k_idx;
    induced->num_induced[k_idx] = prev_induced_degs[idx];
    uintE new_num_induced = induced->num_induced[k_idx];
    // parallel_for(0, new_num_induced, [&] (size_t j) {new_induced[j] =
    // induced->induced_edges[idx * induced->nn + j]; });
    for (size_t j = 0; j < new_num_induced; j++) {
      new_induced[j] = induced->induced_edges[idx * induced->nn + j];
    }
    // parallel_for(0, new_num_induced, [&] (size_t j){
    // induced->labels[new_induced[j]] = k_idx; });
    for (size_t j = 0; j < new_num_induced; j++) {
      induced->labels[new_induced[j]] = k_idx;
    }
    uintE* new_induced_degs = induced->induced_degs + induced->nn * k_idx;
    // parallel_for(0, induced->nn, [&] (size_t j) { new_induced_degs[j] = 0;
    // });
    for (size_t j = 0; j < induced->nn; j++) {
      new_induced_degs[j] = 0;
    }

    for (size_t j = 0; j < new_num_induced; j++) {
      uintE v_idx = new_induced[j];
      uintE v_deg = prev_induced_degs[v_idx];
      uintE* v_edges = induced->induced_edges + v_idx * induced->nn;
      size_t end = v_deg;
      for (size_t l = 0; l < end; l++) {
        if (induced->labels[v_edges[l]] == k_idx)
          new_induced_degs[v_idx]++;
        else {  // if (to_save)
          auto tmp = v_edges[l];
          v_edges[l--] = v_edges[--end];
          v_edges[end] = tmp;
        }
      }
    }

    /*parallel_for(0, new_num_induced, [&] (size_t j) {
      uintE v_idx = new_induced[j];
      uintE v_deg = prev_induced_degs[v_idx];
      uintE* v_edges = induced->induced_edges + v_idx * induced->nn;
      size_t end = v_deg;
      for (size_t l=0; l < end; l++) {
        if (induced->labels[v_edges[l]] == k_idx) new_induced_degs[v_idx]++;
        else { // if (to_save)
          auto tmp = v_edges[l];
          v_edges[l--] = v_edges[--end];
          v_edges[end] = tmp;
        }
      }
    });*/

    auto deg_seq = gbbs::make_slice(new_induced_degs, induced->nn);
    induced->num_edges[k_idx] = parlay::reduce(deg_seq);

    // uintE vtx = prev_induced[i];
    // induced->num_induced[k_idx] = lstintersect_set(prev_induced, num_induced,
    // (uintE*)(DG.get_vertex(vtx).getOutNeighbors()),
    // DG.get_vertex(vtx).out_degree(), true, induced->induced +
    // induced->num_induced[0] * k_idx);
    if (induced->num_induced[k_idx] > 0)
      total_ct += KCliqueDir_fast_orig_rec(DG, k_idx + 1, k, induced);
    // parallel_for(0, new_num_induced, [&] (size_t j){
    // induced->labels[new_induced[j]] = k_idx-1; });
    for (size_t j = 0; j < new_num_induced; j++) {
      induced->labels[new_induced[j]] = k_idx - 1;
    }
  }

  return total_ct;
}

template <class Graph>
inline size_t CountCliques(Graph& DG, size_t k) {
  sequence<size_t> tots = sequence<size_t>::uninitialized(DG.n);
  size_t max_deg = get_max_deg(DG);
  auto init_induced = [&](FullSpace_orig_lw* induced) {
    induced->alloc(max_deg, k, DG.n);
  };
  auto finish_induced = [&](FullSpace_orig_lw* induced) {
    if (induced != nullptr) {
      delete induced;
    }
  };
  parallel_for_alloc<FullSpace_orig_lw>(
      init_induced, finish_induced, 0, DG.n,
      [&](size_t i, FullSpace_orig_lw* induced) {
        if (DG.get_vertex(i).out_degree() != 0) {
          induced->setup(DG, k, i);
          tots[i] = KCliqueDir_fast_orig_rec(DG, 1, k, induced);
        } else
          tots[i] = 0;
      });

  return parlay::reduce(tots);
}

}  // namespace induced_neighborhood
}  // namespace gbbs
