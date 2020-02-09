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

  template <class Graph, class F>
  inline size_t KCliqueDir_fast_rec(Graph& DG, size_t k_idx, size_t k, InducedSpace_lw* induced, F base_f, bool use_base) {
    using W = typename Graph::weight_type;
    size_t num_induced = induced->num_induced[k_idx-1];
    uintE* prev_induced = induced->induced + induced->num_induced[0] * (k_idx - 1);
    if (num_induced == 0) return 0;

    for (size_t i=0; i < num_induced; i++) { induced->intersect[prev_induced[i]] = k_idx; }

    auto tmp = sequence<std::tuple<uintE, W>>(induced->max_deg);
    auto pred = [&] (const uintE& src, const uintE& nbhr, const W& wgh) {
      return (induced->intersect[nbhr] == k_idx);
    };

    if (k_idx + 1 == k) {
      size_t counts = 0;
      for (size_t i=0; i < num_induced; i++) {
        size_t tmp_counts = 0;
        uintE vtx = prev_induced[i];

        auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
          tmp_counts++;
          if (use_base) base_f(std::get<0>(nw), 1);
        };
        DG.get_vertex(vtx).filterOutNgh(vtx, pred, out_f, tmp.begin());
        if (use_base && tmp_counts > 0) base_f(vtx, tmp_counts);
        counts += tmp_counts;
        tmp_counts = 0;

        /*auto map_intersect_f = [&] (const uintE& src, const uintE& nbhr, const W& wgh) {
          if (induced->intersect[nbhr] == k_idx) counts++;
        };
        DG.get_vertex(vtx).mapOutNgh(vtx, map_intersect_f, false);*/
        //counts += lstintersect_set(prev_induced, num_induced, (uintE*)(DG.get_vertex(vtx).getOutNeighbors()), DG.get_vertex(vtx).getOutDegree(), false, nullptr);
      }
      for (size_t i=0; i < num_induced; i++) { induced->intersect[prev_induced[i]] = k_idx - 1; }
      return counts;
    }

    size_t total_ct = 0;
    for (size_t i=0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      uintE* out = induced->induced + induced->num_induced[0] * k_idx;
      uintE count = 0;

      auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
        out[count] = std::get<0>(nw);
        count++;
      };
      DG.get_vertex(vtx).filterOutNgh(vtx, pred, out_f, tmp.begin());

      induced->num_induced[k_idx] = count;
      //induced->num_induced[k_idx] = lstintersect_set(prev_induced, num_induced, (uintE*)(DG.get_vertex(vtx).getOutNeighbors()), DG.get_vertex(vtx).getOutDegree(), true, induced->induced + induced->num_induced[0] * k_idx);
      if (induced->num_induced[k_idx] > 0) {
        auto curr_ct = KCliqueDir_fast_rec(DG, k_idx + 1, k, induced, base_f, use_base);
        if (use_base && curr_ct > 0) base_f(vtx, curr_ct);
        total_ct += curr_ct;
      }
    }

    for (size_t i=0; i < num_induced; i++) { induced->intersect[prev_induced[i]] = k_idx - 1; }
    return total_ct;
  }

  template <class Graph, class G, class F>
  inline size_t CountCliques(Graph& DG, size_t k, G use_f, F base_f, bool use_base=false) {
    using W = typename Graph::weight_type;
    sequence<size_t> tots = sequence<size_t>::no_init(DG.n);

    size_t max_deg = get_max_deg(DG);
    auto init_induced = [&](InducedSpace_lw* induced) { induced->alloc(max_deg, k, DG.n); };
    auto finish_induced = [&](InducedSpace_lw* induced) { if (induced != nullptr) { delete induced; } };
    parallel_for_alloc<InducedSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, InducedSpace_lw* induced) {
      if (DG.get_vertex(i).getOutDegree() != 0) {
        //induced->num_induced[0] = (uintE) DG.get_vertex(i).getOutDegree();
        //for (size_t j=0; j < induced->num_induced[0]; j++) {
        size_t j = 0;
        auto map_intersect_f = [&] (const uintE& src, const uintE& nbhr, const W& wgh) {
          if (use_f(i, nbhr)) {
            induced->induced[j] = nbhr;
            j++;
          }
        };
        DG.get_vertex(i).mapOutNgh(i, map_intersect_f, false);
        induced->num_induced[0] = (uintE) j;
        tots[i] = KCliqueDir_fast_rec(DG, 1, k, induced, base_f, use_base);
        if (use_base && tots[i] > 0) base_f(i, tots[i]);
      } else tots[i] = 0;
    } );

    return pbbslib::reduce_add(tots);
  }





  template <class Graph, class F>
  inline size_t KCliqueDir_simple(Graph& DG, size_t k_idx, size_t k, SimpleSpace* induced, F base_f, bool use_base) {
    using W = typename Graph::weight_type;
    size_t num_induced = induced->num_induced[k_idx-1];
    uintE* prev_induced = induced->induced + induced->num_induced[0] * (k_idx - 1);
    if (num_induced == 0) return 0;

    if (k_idx + 1 == k) {
      size_t counts = 0;
      for (size_t i=0; i < num_induced; i++) {
        size_t tmp_counts = 0;
        uintE vtx = prev_induced[i];
// TODO intersect prev_induced w/nbhrs of vtx in merge intersect
     size_t v_deg = DG.get_vertex(vtx).getOutDegree();
      auto v_iter = DG.get_vertex(vtx).getOutIter(vtx);
      size_t i_iter_idx = 0;
      size_t v_iter_idx = 0;
      while (i_iter_idx < num_induced && v_iter_idx < v_deg) {
        if (prev_induced[i_iter_idx] == std::get<0>(v_iter.cur())) {
          tmp_counts++;
          if (use_base) base_f(prev_induced[i_iter_idx], 1);
          i_iter_idx++; v_iter_idx++;
          if (v_iter.has_next()) v_iter.next();
        } else if (prev_induced[i_iter_idx] < std::get<0>(v_iter.cur())) i_iter_idx++;
        else {
          v_iter_idx++;
          if (v_iter.has_next()) v_iter.next();
        }
      }
        if (use_base && tmp_counts > 0) base_f(vtx, tmp_counts);
        counts += tmp_counts;
      }
      return counts;
    }

    size_t total_ct = 0;
    for (size_t i=0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      uintE* out = induced->induced + induced->num_induced[0] * k_idx;
      uintE count = 0;
// TODO intersect prev_induced w/nbhrs of vtx in merge intersect
     size_t v_deg = DG.get_vertex(vtx).getOutDegree();
      auto v_iter = DG.get_vertex(vtx).getOutIter(vtx);
      size_t i_iter_idx = 0;
      size_t v_iter_idx = 0;
      while (i_iter_idx < num_induced && v_iter_idx < v_deg) {
        if (prev_induced[i_iter_idx] == std::get<0>(v_iter.cur())) {
          out[count] = prev_induced[i_iter_idx];
          count++;
          i_iter_idx++; v_iter_idx++;
          if (v_iter.has_next()) v_iter.next();
        } else if (prev_induced[i_iter_idx] < std::get<0>(v_iter.cur())) i_iter_idx++;
        else {
          v_iter_idx++;
          if (v_iter.has_next()) v_iter.next();
        }
      }
      induced->num_induced[k_idx] = count;
      if (count > 0) {
        auto curr_ct = KCliqueDir_fast_rec(DG, k_idx + 1, k, induced, base_f, use_base);
        if (use_base && curr_ct > 0) base_f(vtx, curr_ct);
        total_ct += curr_ct;
      }
    }
    return total_ct;
  }
} // induced_intersection
