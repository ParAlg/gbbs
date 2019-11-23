#pragma once

#include "intersect.h"
#include "pbbslib/seq.h"

/* TODO: describe what this file does */

namespace induced_intersection {
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
  inline size_t KCliqueDir_fast_rec(Graph& DG, size_t k_idx, size_t k, InducedSpace_lw* induced) {
    size_t num_induced = induced->num_induced[k_idx-1];
    uintE* prev_induced = induced->induced + induced->num_induced[0] * (k_idx - 1);
    if (num_induced == 0) return 0;

    for (size_t i=0; i < num_induced; i++) { induced->intersect[prev_induced[i]] = k_idx; }

    if (k_idx + 1 == k) {
      size_t counts = 0;
      for (size_t i=0; i < num_induced; i++) {
        uintE vtx = prev_induced[i];
        uintE* intersect = (uintE*)(DG.get_vertex(vtx).getOutNeighbors());
        for (size_t j=0; j < DG.get_vertex(vtx).getOutDegree(); j++) {
          if (induced->intersect[intersect[j]] == k_idx) counts++;
        }
        //counts += lstintersect_set(prev_induced, num_induced, (uintE*)(DG.get_vertex(vtx).getOutNeighbors()), DG.get_vertex(vtx).getOutDegree(), false, nullptr);
      }
      for (size_t i=0; i < num_induced; i++) { induced->intersect[prev_induced[i]] = k_idx - 1; }
      return counts;
    }

    size_t total_ct = 0;
    for (size_t i=0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      uintE* intersect = (uintE*)(DG.get_vertex(vtx).getOutNeighbors());
      uintE* out = induced->induced + induced->num_induced[0] * k_idx;
      uintE count = 0;
      for (size_t j=0; j < DG.get_vertex(vtx).getOutDegree(); j++) {
        if (induced->intersect[intersect[j]] == k_idx) {
          out[count] = intersect[j];
          count++;
        }
      }
      induced->num_induced[k_idx] = count;
      //induced->num_induced[k_idx] = lstintersect_set(prev_induced, num_induced, (uintE*)(DG.get_vertex(vtx).getOutNeighbors()), DG.get_vertex(vtx).getOutDegree(), true, induced->induced + induced->num_induced[0] * k_idx);
      if (induced->num_induced[k_idx] > 0) total_ct += KCliqueDir_fast_rec(DG, k_idx + 1, k, induced);
    }

    for (size_t i=0; i < num_induced; i++) { induced->intersect[prev_induced[i]] = k_idx - 1; }
    return total_ct;
  }

  template <class Graph>
  inline size_t CountCliques(Graph& DG, size_t k) {
    //sequence<size_t> tots = sequence<size_t>::no_init(DG.n);
    size_t n = 0;
    size_t max_deg = get_max_deg(DG);
    /*InducedSpace_lw** induceds = (InducedSpace_lw**) malloc(num_workers()*sizeof(InducedSpace_lw*));
    parallel_for (0, num_workers(), [&](size_t i) {
      induceds[i] = new InducedSpace_lw(k, max_deg, DG.n);
    });*/
    InducedSpace_lw* induced = nullptr;
    #pragma omp parallel private(induced) reduction(+:n)
    {
      induced = new InducedSpace_lw(max_deg, k, DG.n);
      #pragma omp for schedule(dynamic, 1) nowait
      for (size_t i=0; i < DG.n; i++) {
    //parallel_for (0, DG.n, [&](size_t i) {
      //InducedSpace_lw* induced = induceds[worker_id()];
      if (DG.get_vertex(i).getOutDegree() != 0) {
        induced->num_induced[0] = (uintE) DG.get_vertex(i).getOutDegree();
        for  (size_t j=0; j < induced->num_induced[0]; j++) {
          induced->induced[j] = ((uintE*)(DG.get_vertex(i).getOutNeighbors()))[j];
        }
        n += KCliqueDir_fast_rec(DG, 1, k, induced);
      } //else tots[i] = 0;
    }//);
    }

    if (induced != nullptr) {induced->del(); delete induced;}

    //parallel_for (0, num_workers(), [&](size_t i) {
    //  induceds[i]->del(); delete induceds[i];
    //});
    //free(induceds);

    //size_t  total = 0;
    //for (size_t i=0; i < DG.n ;i++) { total += tots[i]; }

    return n; //total; //pbbslib::reduce_add(tots);
  }
} // induced_intersection
