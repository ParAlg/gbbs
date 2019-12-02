#pragma once

/* TODO: describe what this file does */

#include "intersect.h"
#include "pbbslib/seq.h"

namespace induced_hybrid {

  template<class Graph>
  size_t get_max_deg(Graph& DG) {
    size_t max_deg = 0;
    parallel_for(0, DG.n, [&] (size_t i) {
      size_t deg = DG.get_vertex(i).getOutDegree();
      pbbs::write_min(&max_deg, deg, std::greater<size_t>());
    });
    return max_deg;
  }

  template <class Graph>
  inline size_t KCliqueDir_fast_hybrid_rec(Graph& DG, size_t k_idx, size_t k, HybridSpace_lw* induced) {
    size_t num_induced = induced->num_induced[k_idx-1];
    uintE* prev_induced = induced->induced + induced->nn * (k_idx - 1);
    if (num_induced == 0) return 0;

    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx; }

    if (k_idx + 1 == k) {
      size_t counts = 0;
      for (size_t i=0; i < num_induced; i++) {
        uintE vtx = prev_induced[i];
        //  get neighbors of vtx
        uintE* intersect = induced->induced_edges + vtx * induced->nn;
        for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
          if (induced->labels[intersect[j]] == k_idx) counts++;
        }
      }
      for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
      return counts;
    }

    size_t total_ct = 0;
    for (size_t i=0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      uintE* intersect = induced->induced_edges + vtx * induced->nn;
      uintE* out = induced->induced + induced->num_induced[0] * k_idx;
      uintE count = 0;
      for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
        if (induced->labels[intersect[j]] == k_idx) {
          out[count] = intersect[j];
          count++;
        }
      }
      induced->num_induced[k_idx] = count;
      if (induced->num_induced[k_idx] > k - k_idx - 1) total_ct += KCliqueDir_fast_hybrid_rec(DG, k_idx + 1, k, induced);
    }

    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
    return total_ct;
  }

  template <class Graph>
  inline size_t CountCliques(Graph& DG, size_t k) {
    //size_t n = 0;
    sequence<size_t> tots = sequence<size_t>::no_init(DG.n);
    size_t max_deg = get_max_deg(DG);
    auto init_induced = [&]() {return new HybridSpace_lw(max_deg, k, DG.n);};
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { induced->del(); delete induced; } };
    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
      if (DG.get_vertex(i).getOutDegree() != 0) {
        induced->setup(DG, k, i);
        tots[i] = KCliqueDir_fast_hybrid_rec(DG, 1, k, induced);
      } else tots[i] = 0;
    } );
    /*
    #pragma omp parallel private(induced) reduction(+:n)
    {
    induced = new HybridSpace_lw(max_deg, k);
    #pragma omp for schedule(dynamic, 1) nowait
    for (size_t i=0; i < DG.n; ++i) {
      if (DG.get_vertex(i).getOutDegree() != 0) {
        induced->setup(DG, k, i);
        n += KCliqueDir_fast_hybrid_rec(DG, 1, k, induced);
      }
    }

    if (induced != nullptr) { induced->del(); delete induced; }

    }*/

    return pbbslib::reduce_add(tots);
  }

} // namespace induced_neighborhood
