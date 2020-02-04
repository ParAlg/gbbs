#pragma once

#include "intersect.h"
#include "pbbslib/seq.h"

/* TODO: describe what this file does */

namespace induced_split {
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
  inline size_t KCliqueDir_fast_hybrid_rec(Graph& DG, size_t k_idx, size_t k, SplitSpace* induced) {
    //if (k == 2) return induced->num_edges;
    size_t num_induced = induced->hybrid_space->num_induced[k_idx-1];
    if (num_induced == 0) return 0;
    uintE* prev_induced = induced->hybrid_space->induced + induced->nn * (k_idx - 1);

    for (size_t i=0; i < num_induced; i++) { induced->hybrid_space->labels[prev_induced[i]] = k_idx; }

    if (k_idx + 1 == k) {
      size_t counts = 0;
      for (size_t i=0; i < num_induced; i++) {
        uintE vtx = prev_induced[i];
        //  get neighbors of vtx
        uintE* intersect = induced->hybrid_space->induced_edges + vtx * induced->nn;
        size_t tmp_counts = 0;
        for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
          if (induced->hybrid_space->labels[intersect[j]] == k_idx) {
            tmp_counts++;
            //if (induced->use_base) base_f(induced->relabel[intersect[j]], 1);
          } 
        }
        //if (induced->use_base) base_f(induced->relabel[vtx], tmp_counts);
        counts += tmp_counts;
      }
      for (size_t i=0; i < num_induced; i++) { induced->hybrid_space->labels[prev_induced[i]] = k_idx - 1; }
      return counts;
    }

    size_t total_ct = 0;
    for (size_t i=0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      //TODO problem w/storing base -- we've relabeled our vert w/relabeling: check base is correct
      uintE* intersect = induced->hybrid_space->induced_edges + vtx * induced->nn;
      uintE* out = induced->hybrid_space->induced + induced->hybrid_space->num_induced[0] * k_idx;
      uintE count = 0;
      for (size_t j=0; j < induced->hybrid_space->induced_degs[vtx]; j++) {
        if (induced->hybrid_space->labels[intersect[j]] == k_idx) {
          out[count] = intersect[j];
          count++;
        }
      }
      induced->hybrid_space->num_induced[k_idx] = count;
      if (induced->hybrid_space->num_induced[k_idx] > k - k_idx - 1) {
        auto curr_counts = KCliqueDir_fast_hybrid_rec(DG, k_idx + 1, k, induced);
        total_ct += curr_counts;
      }
    }

    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
    return total_ct;
  }

  template <class Graph>
  inline size_t KCliqueDir_fast_rec(Graph& DG, size_t k_idx, size_t k, SplitSpace* induced) {
    using W = typename Graph::weight_type;
    size_t num_induced = induced->induced_space->num_induced[k_idx-1];
    uintE* prev_induced = induced->induced_space->induced + induced->induced_space->num_induced[0] * (k_idx - 1);
    if (num_induced == 0) return 0;

    for (size_t i=0; i < num_induced; i++) { induced->induced_space->intersect[prev_induced[i]] = k_idx; }

    if (k_idx + 1 == k) {
      auto tmp = sequence<std::tuple<uintE, W>>(num_induced);
      auto pred = [&] (const uintE& src, const uintE& nbhr, const W& wgh) {
        return (induced->induced_space->intersect[nbhr] == k_idx);
      };
      size_t counts = 0;
      for (size_t i=0; i < num_induced; i++) {
        uintE vtx = prev_induced[i];

        auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
          counts++;
        };
        DG.get_vertex(vtx).filterOutNgh(vtx, pred, out_f, tmp.begin());
      }
      for (size_t i=0; i < num_induced; i++) { induced->induced_space->intersect[prev_induced[i]] = k_idx - 1; }
      return counts;
    }

    if (num_induced * (k-k_idx) <= induced->k_threshold) {
      induced->switch_alloc(DG, k-k_idx+1, DG.n);
      return KCliqueDir_fast_hybrid_rec(DG, 1, k-k_idx+1, induced);
    }

    auto tmp = sequence<std::tuple<uintE, W>>(num_induced);
    auto pred = [&] (const uintE& src, const uintE& nbhr, const W& wgh) {
      return (induced->induced_space->intersect[nbhr] == k_idx);
    };

    size_t total_ct = 0;
    for (size_t i=0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      uintE* out = induced->induced_space->induced + induced->induced_space->num_induced[0] * k_idx;
      uintE count = 0;

      auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw) {
        out[count] = std::get<0>(nw);
        count++;
      };
      DG.get_vertex(vtx).filterOutNgh(vtx, pred, out_f, tmp.begin());


      induced->induced_space->num_induced[k_idx] = count;
      if (induced->induced_space->num_induced[k_idx] > 0) total_ct += KCliqueDir_fast_rec(DG, k_idx + 1, k, induced);
    }

    for (size_t i=0; i < num_induced; i++) { induced->induced_space->intersect[prev_induced[i]] = k_idx - 1; }
    return total_ct;
  }

  template <class Graph>
  inline size_t CountCliques(Graph& DG, size_t k, F base_f, bool use_base=false, bool label=true, size_t k_threshold=0) {
    using W = typename Graph::weight_type;
    sequence<size_t> tots = sequence<size_t>::no_init(DG.n);

    size_t max_deg = get_max_deg(DG);
    auto init_induced = [&](SplitSpace* induced) { induced->alloc(max_deg, k, DG.n, label, use_base, k_threshold); };
    auto finish_induced = [&](SplitSpace* induced) { if (induced != nullptr) { delete induced; } };
    parallel_for_alloc<SplitSpace>(init_induced, finish_induced, 0, DG.n, [&](size_t i, SplitSpace* induced) {
      if (DG.get_vertex(i).getOutDegree() != 0) {
        induced->num_induced[0] = (uintE) DG.get_vertex(i).getOutDegree();
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
