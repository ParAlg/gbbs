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
    using W = typename Graph::weight_type;
    size_t num_induced = induced->num_induced[k_idx-1];
    uintE* prev_induced = induced->induced + induced->num_induced[0] * (k_idx - 1);
    if (num_induced == 0) return 0;

    for (size_t i=0; i < num_induced; i++) { induced->intersect[prev_induced[i]] = k_idx; }

    if (k_idx + 1 == k) {
      size_t counts = 0;
      for (size_t i=0; i < num_induced; i++) {
        uintE vtx = prev_induced[i];
        //auto intersect = DG.get_vertex(vtx).getOutNeighbors();
        //for (size_t j=0; j < DG.get_vertex(vtx).getOutDegree(); j++) {
        auto map_intersect_f = [&] (const uintE& src, const uintE& nbhr, const W& wgh) {
          if (induced->intersect[nbhr] == k_idx) counts++;
        };
        DG.get_vertex(vtx).mapOutNgh(vtx, map_intersect_f, false);
        //counts += lstintersect_set(prev_induced, num_induced, (uintE*)(DG.get_vertex(vtx).getOutNeighbors()), DG.get_vertex(vtx).getOutDegree(), false, nullptr);
      }
      for (size_t i=0; i < num_induced; i++) { induced->intersect[prev_induced[i]] = k_idx - 1; }
      return counts;
    }

    size_t total_ct = 0;
    for (size_t i=0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      //auto intersect = DG.get_vertex(vtx).getOutNeighbors();
      uintE* out = induced->induced + induced->num_induced[0] * k_idx;
      uintE count = 0;
      //for (size_t j=0; j < DG.get_vertex(vtx).getOutDegree(); j++) {
      auto map_intersect_f = [&] (const uintE& src, const uintE& nbhr, const W& wgh) {
        if (induced->intersect[nbhr] == k_idx) {
          out[count] = nbhr;
          count++;
        }
      };
      DG.get_vertex(vtx).mapOutNgh(vtx, map_intersect_f, false);
      induced->num_induced[k_idx] = count;
      //induced->num_induced[k_idx] = lstintersect_set(prev_induced, num_induced, (uintE*)(DG.get_vertex(vtx).getOutNeighbors()), DG.get_vertex(vtx).getOutDegree(), true, induced->induced + induced->num_induced[0] * k_idx);
      if (induced->num_induced[k_idx] > 0) total_ct += KCliqueDir_fast_rec(DG, k_idx + 1, k, induced);
    }

    for (size_t i=0; i < num_induced; i++) { induced->intersect[prev_induced[i]] = k_idx - 1; }
    return total_ct;
  }

  template <class Graph>
  inline size_t CountCliques(Graph& DG, size_t k) {
    using W = typename Graph::weight_type;
    sequence<size_t> tots = sequence<size_t>::no_init(DG.n);

    size_t max_deg = get_max_deg(DG);
    auto init_induced = [&](InducedSpace_lw* induced) { induced->alloc(max_deg, k, DG.n); };
    auto finish_induced = [&](InducedSpace_lw* induced) { if (induced != nullptr) { delete induced; } };
    parallel_for_alloc<InducedSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, InducedSpace_lw* induced) {
      if (DG.get_vertex(i).getOutDegree() != 0) {
        induced->num_induced[0] = (uintE) DG.get_vertex(i).getOutDegree();
        //for (size_t j=0; j < induced->num_induced[0]; j++) {
        size_t j = 0;
        auto map_intersect_f = [&] (const uintE& src, const uintE& nbhr, const W& wgh) {
          induced->induced[j] = nbhr;
          j++;
        };
        DG.get_vertex(i).mapOutNgh(i, map_intersect_f, false);
        tots[i] = KCliqueDir_fast_rec(DG, 1, k, induced);
      } else tots[i] = 0;
    } );

    return pbbslib::reduce_add(tots);
  }
} // induced_intersection
