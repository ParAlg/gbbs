#ifndef _KINTERSECT_
#define _KINTERSECT_

#pragma once

#include <math.h>

#include <limits>

#include "ligra/bucket.h"
#include "ligra/edge_map_reduce.h"
#include "ligra/ligra.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "ligra/pbbslib/sparse_table.h"

#include "external/simdinter/include/common.h"
#include "external/simdinter/include/intersection.h"
#include "external/graphsetinter/src/set_operation.hpp"

#define INDUCED_STACK_THR 5000

// TODO retry using lambdas for intersects
// TODO make intersects more modularized -- have some kind of wrapper that generates the vtx_seq and everything, and just have the intersects
// actually do the intersect, with the save and out_ptr options

// have a wrapper where -- if you insert a pointer, then you take responsibility for that pointer and all you get out
// is the size (and an empty sequence)
// otherwise, if you have no pointer, the sequence gets allocated for you


inline size_t lstintersect_simple(uintE* a, size_t size_a, uintE* b, size_t size_b, bool save, uintE* out) {
  // auto seq_b = pbbslib::make_sequence(b, size_b);
  size_t out_idx = 0;
  for (size_t i=0; i < size_a; i++) {
    size_t j=0;
    for (j=0; j < size_b; j++) {
      if (a[i] == b[j]) break;
    }
    //size_t j = pbbslib::binary_search(seq_b, a[i], std::less<uintE>());
    if (j < size_b) {
      out[out_idx] = b[j];
      out_idx++;
    }
  }
  return out_idx;
}
struct lstintersect_simple_struct {
  bool count_space_flag = true;
  size_t operator()(uintE* a, size_t size_a, uintE* b, size_t size_b, bool save, uintE* out) const {
    return lstintersect_simple(a, size_a, b, size_b, save, out);
  }
};

inline size_t lstintersect_par(uintE* a, size_t size_a, uintE* b, size_t size_b, bool save, uintE* out) {
  auto seq_a = pbbslib::make_sequence<uintE>(a, size_a);
  auto seq_b = pbbslib::make_sequence<uintE>(b, size_b);

  if (!save) {
    auto merge_f = [&] (uintE ngh) {};
    return intersection::seq_merge_full(seq_a, seq_b, merge_f);
  }

  size_t index = 0;
  auto merge_f = [&] (uintE ngh) {
    out[index] = ngh;
    index++;
  };
  return intersection::seq_merge_full(seq_a, seq_b, merge_f);
}

struct lstintersect_par_struct {
  bool count_space_flag = false;
  size_t operator()(uintE* a, size_t size_a, uintE* b, size_t size_b, bool save, uintE* out) const {
    return lstintersect_par(a, size_a, b, size_b, save, out);
  }
};

// TODO radix sort in place
inline size_t lstintersect_vec(uintE* a, size_t size_a, uintE* b, size_t size_b, bool save, uintE* out) {
  SIMDCompressionLib::intersectionfunction inter = SIMDCompressionLib::IntersectionFactory::getFromName("simd");
  size_t out_size = inter(a, size_a, b, size_b, out);
  return out_size;
}

struct lstintersect_vec_struct {
  bool count_space_flag = true;
  size_t operator()(uintE* a, size_t size_a, uintE* b, size_t size_b, bool save, uintE* out) const {
    return lstintersect_vec(a, size_a, b, size_b, save, out);
  }
};

inline size_t lstintersect_set(uintE* a, size_t size_a, uintE* b, size_t size_b, bool save, uintE* out) {
  if (!save) {
#if SIMD_STATE == 2
    return (size_t) intersect_scalar2x_count((int*) a, (int) size_a, (int*) b, (int) size_b);
#elif SIMD_STATE == 4
    return (size_t) intersect_simd4x_count((int*) a, (int) size_a, (int*) b, (int) size_b);
#else
    return (size_t) intersect_count((int*) a, (int) size_a, (int*) b, (int) size_b);
#endif
  }
#if SIMD_STATE == 2
  int out_size = intersect_scalar2x((int*) a, (int) size_a, (int*) b, (int) size_b, (int*) out);
#elif SIMD_STATE == 4
  int out_size = intersect_simd4x((int*) a, (int) size_a, (int*) b, (int) size_b, (int*) out);
#else
  int out_size = intersect((int*) a, (int) size_a, (int*) b, (int) size_b, (int*) out);
#endif
  return (size_t) out_size;
}

struct lstintersect_set_struct {
  bool count_space_flag = false;
  size_t operator()(uintE* a, size_t size_a, uintE* b, size_t size_b, bool save, uintE* out) const {
    return lstintersect_set(a, size_a, b, size_b, save, out);
  }
};



struct InducedSpace_lw {
  // Array storing size of induced neighbor list, for each level of recursion.
  uintE* num_induced = nullptr;
  // Array storing induced neighbor list, for all levels of recursion. Levels go
  // from i=0 ... k-1. For each i, the induced neighbor list at level i is
  // induced + num_induced[0]*i (num_induced[0] is an upper-bound on the maximum
  // number of induced neighbors at any level).
  uintE* induced = nullptr;

  int* intersect = nullptr;
  InducedSpace_lw() {}

  void alloc(size_t max_deg, size_t k, size_t n) {
    if (!induced) induced = (uintE*) malloc(k*max_deg*sizeof(uintE));
    if (!num_induced) num_induced = (uintE*) malloc(k*sizeof(uintE));
    if (!intersect) {
      intersect = (int*) malloc(n*sizeof(int));
      for (size_t  i=0; i < n; i++) {
        intersect[i] = 0;
      }
    }
  }
  void del() {
    if (induced) { free(induced); induced = nullptr; }
    if (num_induced) { free(num_induced); num_induced = nullptr; }
    if (intersect) { free(intersect); intersect = nullptr; }
  }

  ~InducedSpace_lw() { del(); }
};

struct FullSpace_orig_lw {
  uintE* num_induced = nullptr;
  uintE* induced = nullptr;
  uintE* num_edges = 0;
  uintE* induced_edges = nullptr;
  uintE* induced_degs = nullptr;
  uintE* labels = nullptr;
  uintE* old_labels = nullptr;
  size_t nn = 0;

  FullSpace_orig_lw() {}

  void alloc(size_t max_induced, size_t k, size_t n) {
    induced = (uintE*) malloc(sizeof(uintE)*k*max_induced);
    induced_degs = (uintE*) malloc(sizeof(uintE)*k*max_induced);
    labels = (uintE*) malloc(sizeof(uintE)*max_induced);
    induced_edges = (uintE*) malloc(sizeof(uintE)*max_induced*max_induced);
    num_induced = (uintE*) malloc(sizeof(uintE)*k);
    num_edges = (uintE*) malloc(sizeof(uintE)*k);
    if (old_labels == nullptr) old_labels = (uintE*) calloc(n, sizeof(uintE));
  }

  template <class Graph>
  void setup(Graph& DG, size_t k, size_t i) {
    using W = typename Graph::weight_type;
    num_induced[0] = DG.get_vertex(i).getOutDegree();
    nn = num_induced[0];
    
    for (size_t  j=0; j < nn; j++) { induced[j] = j; }
    for (size_t j=0; j < nn; j++) { induced_degs[j] = 0; }
    for (size_t j=0; j < nn; j++)  { labels[j] = 0; }
  
    //auto induced_g = DG.get_vertex(i).getOutNeighbors();
    //for (size_t o=0; o < num_induced[0]; o++) { old_labels[std::get<0>(induced_g[o])] = o + 1; }
    size_t o = 0;
    auto map_label_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
      old_labels[ngh] = o + 1;
      o++;
    };
    DG.get_vertex(i).mapOutNgh(i, map_label_f, false);

    size_t j = 0;
    auto map_f = [&] (const uintE& src, const uintE& v, const W& wgh) {
    //for (size_t j=0; j < nn; j++) {
      //uintE v = std::get<0>(induced_g[j]);
      //auto v_nbhrs = DG.get_vertex(v).getOutNeighbors();
      //size_t v_deg = DG.get_vertex(v).getOutDegree();
      // intersect v_nbhrs from 0 to v_deg with induced_g from 0 to num_induced[0]
      // store result in induced_edges[j*nn]
      // store size in induced_degs[j]
      //for (size_t l=0; l < v_deg; l++) {
      auto map_nbhrs_f = [&] (const uintE& src_v, const uintE& v_nbhr, const W& wgh_v) {
        if (old_labels[v_nbhr] > 0) { //std::get<0>(v_nbhrs[l])
          induced_edges[j*nn + induced_degs[j]] = old_labels[v_nbhr] - 1; //std::get<0>(v_nbhrs[l])
          induced_degs[j]++;
        }
      };
      DG.get_vertex(v).mapOutNgh(v, map_nbhrs_f, false);
      //}
    //}
      j++;
    };
    DG.get_vertex(i).mapOutNgh(i, map_f, false);

    //for (size_t o=0; o < num_induced[0]; o++) { old_labels[std::get<0>(induced_g[o])] = 0; }
    auto map_relabel_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
      old_labels[ngh] = 0;
    };
    DG.get_vertex(i).mapOutNgh(i, map_relabel_f, false);

    auto deg_seq = pbbslib::make_sequence(induced_degs, nn);
    num_edges[0] = pbbslib::reduce_add(deg_seq);
  }

  static void init(){}
  static void finish(){}

  void del() {
    if (labels) { free(labels); labels = nullptr; }
    if (induced) { free(induced); induced = nullptr; }
    if (induced_edges) { free(induced_edges); induced_edges = nullptr; }
    if (induced_degs) { free(induced_degs); induced_degs = nullptr; }
    if (num_edges) { free(num_edges); num_edges = nullptr; }
    if (num_induced) { free(num_induced); num_induced = nullptr; }
    if (old_labels) { free(old_labels); old_labels=nullptr; }
  }

  ~FullSpace_orig_lw() { del(); }

};



struct HybridSpace_lw {
  uintE* num_induced = nullptr;
  uintE* induced = nullptr;
  uintE num_edges = 0;
  uintE* induced_edges = nullptr;
  uintE* induced_degs = nullptr;
  char* labels = nullptr;
  uintE* old_labels = nullptr;
  size_t nn = 0;
  bool use_old_labels = true;
  bool free_relabel = true;
  uintE* relabel = nullptr;
  bool use_base = false;
  HybridSpace_lw () {}

  void alloc(size_t max_induced, size_t k, size_t n, bool _use_old_labels, bool _use_base, bool _free_relabel=true) {
    use_old_labels = _use_old_labels;
    use_base = _use_base;
    free_relabel = _free_relabel;
    if (induced == nullptr) induced = (uintE*) malloc(sizeof(uintE)*k*max_induced);
    if (free_relabel && induced_degs == nullptr) induced_degs = (uintE*) malloc(sizeof(uintE)*max_induced);
    if (labels == nullptr) labels = (char*) calloc(max_induced, sizeof(char));
    if (free_relabel && induced_edges == nullptr) induced_edges = (uintE*) malloc(sizeof(uintE)*max_induced*max_induced);
    if (num_induced == nullptr) num_induced = (uintE*) malloc(sizeof(uintE)*k);
    if (free_relabel && use_old_labels && old_labels == nullptr) old_labels = (uintE*) calloc(n, sizeof(uintE));
    if (free_relabel && use_base && relabel == nullptr) { relabel = (uintE*) malloc(sizeof(uintE)*max_induced);}
  }
  
  void copy(HybridSpace_lw* space) {
    if (use_base) relabel = space->relabel;
    induced_edges = space->induced_edges;
    induced_degs = space->induced_degs;
    nn = space->nn;
  }


  // f should denote if a vert is active or not
  template <class Graph, class F>
  void setup(Graph& DG, size_t k, size_t i, F f) {
    //if (use_base) base[0] = i;
    if (use_old_labels) setup_labels(DG, k, i, f);
    else setup_intersect(DG, k, i, f);
  }

  template <class Graph>
  void setup(Graph& DG, size_t k, size_t i) {
    //if (use_base) base[0] = i;
    auto f = [&](const uintE& src, const uintE& u) { return true; };
    if (use_old_labels) setup_labels(DG, k, i, f);
    else setup_intersect(DG, k, i, f);
  }

  template <class Graph, class F>
  void setup_intersect(Graph& DG, size_t k, size_t i, F f) {
    using W = typename Graph::weight_type;
    if (use_base) {
      size_t j = 0;
      auto map_base_f = [&] (const uintE& src, const uintE& v, const W& wgh) { relabel[j] = v; j++;};
      DG.get_vertex(i).mapOutNgh(i, map_base_f, false);
    }

    nn = DG.get_vertex(i).getOutDegree();

    //for (size_t j=0; j < nn; j++) { induced_degs[j] = 0; }
    parallel_for(0, nn, [&] (size_t j) { induced_degs[j] = 0; });
  
    num_induced[0] = nn;
    //for (size_t  j=0; j < nn; j++) { induced[j] = j; }
    parallel_for(0, nn, [&] (size_t j) { induced[j] = j; });

    size_t j = 0;
    auto map_f = [&] (const uintE& src, const uintE& v, const W& wgh) {
      if (!f(src, v)) { j++; return; }
      size_t v_deg = DG.get_vertex(v).getOutDegree();
      // intersect v_nbhrs from 0 to v_deg with induced_g from 0 to num_induced[0]
      // store result in induced_edges[j*nn]
      // store size in induced_degs[j]
      auto i_iter = DG.get_vertex(i).getOutIter(i);
      auto v_iter = DG.get_vertex(v).getOutIter(v);
      size_t i_iter_idx = 0;
      size_t v_iter_idx = 0;

      while (i_iter_idx < nn && v_iter_idx < v_deg) {
        if (std::get<0>(i_iter.cur()) == std::get<0>(v_iter.cur())) {
          if (f(i, std::get<0>(i_iter.cur())) && f(v, std::get<0>(i_iter.cur()))) {
            induced_edges[j*nn + induced_degs[j]] = i_iter_idx;
            induced_degs[j]++;
          }
          i_iter_idx++; v_iter_idx++;
          if (i_iter.has_next()) i_iter.next();
          if (v_iter.has_next()) v_iter.next();
        } else if (std::get<0>(i_iter.cur()) < std::get<0>(v_iter.cur())) {
          i_iter_idx++;
          if (i_iter.has_next()) i_iter.next();
        }
        else {
          v_iter_idx++;
          if (v_iter.has_next()) v_iter.next();
        }
      }
      j++;
    };
    DG.get_vertex(i).mapOutNgh(i, map_f, false);

    auto deg_seq = pbbslib::make_sequence(induced_degs, nn);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  template <class Graph, class F>
  void setup_labels(Graph& DG, size_t k, size_t i, F f) {
    using W = typename Graph::weight_type;
    nn = DG.get_vertex(i).getOutDegree();
    //for (size_t j=0; j < nn; j++) { induced_degs[j] = 0; }
    parallel_for(0, nn, [&] (size_t j) { induced_degs[j] = 0; });
  
    num_induced[0] = nn;
    //for (size_t  j=0; j < nn; j++) { induced[j] = j; }
    parallel_for(0, nn, [&] (size_t j) { induced[j] = j; });

    size_t o = 0;
    auto map_label_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
      if (!f(src, ngh)) {o++; return;}
      old_labels[ngh] = o + 1;
      if (use_base) { relabel[o] = ngh; }
      o++;
    };
    DG.get_vertex(i).mapOutNgh(i, map_label_f, false);
  

    size_t j = 0;
    auto map_f = [&] (const uintE& src, const uintE& v, const W& wgh) {
      if (!f(src, v)) { j++; return; }
      // intersect v_nbhrs from 0 to degree(v) with induced_g from 0 to num_induced[0]
      // store result in induced_edges[j*nn]
      // store size in induced_degs[j]
      auto map_nbhrs_f = [&] (const uintE& src_v, const uintE& v_nbhr, const W& wgh_v) {
        if (!f(src_v, v_nbhr)) return;
        if (old_labels[v_nbhr] > 0) {
          induced_edges[j*nn + induced_degs[j]] = old_labels[v_nbhr] - 1;
          induced_degs[j]++;
        }
      };
      DG.get_vertex(v).mapOutNgh(v, map_nbhrs_f, false);
      j++;
    };
    DG.get_vertex(i).mapOutNgh(i, map_f, false);

    auto map_relabel_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
      old_labels[ngh] = 0;
    };
    DG.get_vertex(i).mapOutNgh(i, map_relabel_f, false);

    auto deg_seq = pbbslib::make_sequence(induced_degs, nn);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  template <class Graph, class Graph2, class F>
  void setup(Graph& G, Graph2& DG, size_t k, size_t i, F f, char* still_active) {
    setup_labels(G, DG, k, i, f, still_active);
  }

  template <class Graph, class Graph2, class F>
  void setup_labels(Graph& G, Graph2& DG, size_t k, size_t i, F f, char* still_active) {
    using W = typename Graph::weight_type;
    nn = G.get_vertex(i).getOutDegree();
    for (size_t j=0; j < nn; j++) { induced_degs[j] = 0; }
    num_induced[0] = nn;
    for (size_t  j=0; j < nn; j++) { induced[j] = j; }
    size_t o = 0;
    auto map_label_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
      if (!f(src, ngh)) {o++; return;}
      old_labels[ngh] = o + 1;
      if (use_base) { relabel[o] = ngh; }
      o++;
    };
    G.get_vertex(i).mapOutNgh(i, map_label_f, false);
  
    size_t j = 0;
    auto map_f = [&] (const uintE& src, const uintE& v, const W& wgh) {
      if (!f(src, v)) { j++; return; }
      auto map_nbhrs_f = [&] (const uintE& src_v, const uintE& v_nbhr, const W& wgh_v) {
        if (old_labels[v_nbhr] > 0) { //still_active[v_nbhr] != 2 && 
          induced_edges[j*nn + induced_degs[j]] = old_labels[v_nbhr] - 1;
          induced_degs[j]++;
        }
      };
      DG.get_vertex(v).mapOutNgh(v, map_nbhrs_f, false);
      j++;
    };
    G.get_vertex(i).mapOutNgh(i, map_f, false);

    auto map_relabel_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
      old_labels[ngh] = 0;
    };
    G.get_vertex(i).mapOutNgh(i, map_relabel_f, false);

    auto deg_seq = pbbslib::make_sequence(induced_degs, nn);
    num_edges = pbbslib::reduce_add(deg_seq);
  }





  template <class Graph, class F>
  void setup_edge(Graph& DG, size_t k, size_t i, size_t l, F f) {
    //if (use_base) base[0] = i;
    if (use_old_labels) setup_labels_edge(DG, k, i, l, f);
    else setup_intersect_edge(DG, k, i, l, f);
  }

  template <class Graph>
  void setup_edge(Graph& DG, size_t k, size_t i, size_t l) {
    auto f = [&](const uintE& src, const uintE& u) { return true; };
    if (use_old_labels) setup_labels_edge(DG, k, i, l, f);
    else setup_intersect_edge(DG, k, i, l, f);
  }

 template <class Graph, class F>
  void setup_labels_edge(Graph& DG, size_t k, size_t i, size_t l, F f) {
  using W = typename Graph::weight_type;
    auto map_label_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
      if (!f(src, ngh)) return;
      old_labels[ngh] = UINT_E_MAX;
    };
    DG.get_vertex(i).mapOutNgh(i, map_label_f, false);
    size_t o = 0;
    auto lmap_label_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
      if (!f(src, ngh)) {return;}
      if (old_labels[ngh] == UINT_E_MAX) {
        induced[o] = ngh;
        old_labels[ngh] = o + 1;
        if (use_base) { relabel[o] = ngh; }
        o++;
      }
    };
    DG.get_vertex(l).mapOutNgh(l, lmap_label_f, false);
    auto remap_label_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) { if (old_labels[ngh] == UINT_E_MAX) old_labels[ngh] = 0; };
    DG.get_vertex(i).mapOutNgh(i, remap_label_f, true);
    nn = o; num_induced[0] = o;
    parallel_for (0, nn, [&] (size_t p) { induced_degs[p] = 0; });
  
    for (size_t j=0; j < nn; j++) {
      // intersect each vertex in induced against indced
      auto v = induced[j];
      auto map_nbhrs_f = [&] (const uintE& src_v, const uintE& v_nbhr, const W& wgh_v) {
        if (!f(src_v, v_nbhr)) return;
        if (old_labels[v_nbhr] > 0) {
          induced_edges[j*nn + induced_degs[j]] = old_labels[v_nbhr] - 1;
          induced_degs[j]++;
        }
      };
      DG.get_vertex(v).mapOutNgh(v, map_nbhrs_f, false);
    }

    for (size_t p=0; p < nn; p++) { induced[p] = p; }
    auto lremap_label_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) { old_labels[ngh] = 0; };
    DG.get_vertex(l).mapOutNgh(l, lremap_label_f, true);
    auto deg_seq = pbbslib::make_sequence(induced_degs, nn);
    num_edges = pbbslib::reduce_add(deg_seq);

  }

template <class Graph, class F>
  void setup_intersect_edge(Graph& DG, size_t k, size_t i, size_t l, F f) {
    size_t j = 0;
    size_t v_deg = DG.get_vertex(l).getOutDegree();
    size_t i_deg = DG.get_vertex(i).getOutDegree();
    auto i_iter = DG.get_vertex(i).getOutIter(i);
    auto v_iter = DG.get_vertex(l).getOutIter(l);
    size_t i_iter_idx = 0;
    size_t v_iter_idx = 0;
    while (i_iter_idx < i_deg && v_iter_idx < v_deg) {
      if (std::get<0>(i_iter.cur()) == std::get<0>(v_iter.cur())) {
        if (f(i, std::get<0>(i_iter.cur())) && f(l, std::get<0>(i_iter.cur()))) {
          induced[j] = i_iter_idx;
          j++;
        }
        i_iter_idx++; v_iter_idx++;
        if (i_iter.has_next()) i_iter.next();
        if (v_iter.has_next()) v_iter.next();
      } else if (std::get<0>(i_iter.cur()) < std::get<0>(v_iter.cur())) {
        i_iter_idx++;
        if (i_iter.has_next()) i_iter.next();
      }
      else {
        v_iter_idx++;
        if (v_iter.has_next()) v_iter.next();
      }
    }
    nn = j; num_induced[0] = j;
    if (use_base) {
      parallel_for(0, nn, [&] (size_t p) { relabel[p] = induced[p]; });
    }
    parallel_for (0, nn, [&] (size_t p) { induced_degs[p] = 0; });
  
    for (size_t p=0; p < nn; p++) {
      // intersect each vertex in induced against indced
      size_t u_deg = DG.get_vertex(induced[p]).getOutDegree();
      auto u_iter = DG.get_vertex(induced[p]).getOutIter(induced[p]);
      size_t u_iter_idx = 0;
      i_iter_idx = 0;
      while (i_iter_idx < nn && u_iter_idx < u_deg) {
        if (induced[i_iter_idx] == std::get<0>(u_iter.cur())) {
          if (f(induced[p], induced[i_iter_idx])) {
            induced_edges[p*nn + induced_degs[p]] = i_iter_idx;
            induced_degs[p]++;
          }
          i_iter_idx++; u_iter_idx++;
          if (u_iter.has_next()) u_iter.next();
        } else if (induced[i_iter_idx] < std::get<0>(u_iter.cur())) {
          i_iter_idx++;
        }
        else {
          u_iter_idx++;
          if (u_iter.has_next()) u_iter.next();
        }
      }
    }

    for (size_t p=0; p < nn; j++) { induced[p] = p; }

    auto deg_seq = pbbslib::make_sequence(induced_degs, nn);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  static void init(){}
  static void finish(){}

  void del() {
    if (labels) {free(labels); labels=nullptr;}
    if (induced) {free(induced); induced=nullptr;}
    if (free_relabel && induced_edges) {free(induced_edges); induced_edges=nullptr;}
    if (free_relabel && induced_degs) {free(induced_degs); induced_degs=nullptr;}
    if (num_induced) {free(num_induced); num_induced=nullptr;}
    if (use_old_labels && old_labels) {free(old_labels); old_labels=nullptr;}
    if (use_base && relabel && free_relabel) {free(relabel); relabel=nullptr;}
    if (!free_relabel) {relabel = nullptr; induced_edges=nullptr; induced_degs=nullptr;}
  }

  ~HybridSpace_lw() { del(); }

};







typedef struct {
	uintE key;
	long value;
} keyvalueLLU;

typedef struct {
	size_t n_max;// max number of nodes.
	size_t n;// number of nodes.
	size_t *pt;// pointers to nodes.
	keyvalueLLU *kv;// (node,nck)
} bheapLLU;

bheapLLU *constructLLU(size_t n_max){
	bheapLLU *heap=(bheapLLU*)malloc(sizeof(bheapLLU));

	heap->n_max=n_max;
	heap->n=0;
	heap->pt=(size_t*)malloc(n_max*sizeof(size_t));
	for (size_t i=0;i<n_max;i++) heap->pt[i]=std::numeric_limits<size_t>::max();
	heap->kv=(keyvalueLLU*)malloc(n_max*sizeof(keyvalueLLU));
	return heap;
}

inline void swapLLU(bheapLLU *heap,unsigned i, unsigned j) {
	keyvalueLLU kv_tmp=heap->kv[i];
	auto pt_tmp=heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key]=heap->pt[heap->kv[j].key];
	heap->kv[i]=heap->kv[j];
	heap->pt[heap->kv[j].key]=pt_tmp;
	heap->kv[j]=kv_tmp;
}

inline void bubble_upLLU(bheapLLU *heap,unsigned i) {
	unsigned j=(i-1)/2;
	while (i>0) {
		if (heap->kv[j].value>heap->kv[i].value) {
			swapLLU(heap,i,j);
			i=j;
			j=(i-1)/2;
		}
		else break;
	}
}

inline void bubble_downLLU(bheapLLU *heap) {
	unsigned i=0,j1=1,j2=2,j;
	while (j1<heap->n) {
		j=( (j2<heap->n) && (heap->kv[j2].value<heap->kv[j1].value) ) ? j2 : j1 ;
		if (heap->kv[j].value < heap->kv[i].value) {
			swapLLU(heap,i,j);
			i=j;
			j1=2*i+1;
			j2=j1+1;
			continue;
		}
		break;
	}
}

inline void insertLLU(bheapLLU *heap,keyvalueLLU kv){
	heap->pt[kv.key]=(heap->n)++;
	heap->kv[heap->n-1]=kv;
	bubble_upLLU(heap,heap->n-1);
}

inline void updateLLU(bheapLLU *heap,unsigned key,unsigned long long delta){
	unsigned i=heap->pt[key];
	if (i!=std::numeric_limits<size_t>::max()){
		((heap->kv[i]).value)-=delta;
		bubble_upLLU(heap,i);
	}
}

inline keyvalueLLU popminLLU(bheapLLU *heap){
	keyvalueLLU min=heap->kv[0];
	heap->pt[min.key]=std::numeric_limits<size_t>::max();
	heap->kv[0]=heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key]=0;
	bubble_downLLU(heap);
	return min;
}

//Building the heap structure with (key,value)=(node,k-clique degree) for each node
bheapLLU* mkheapLLU(size_t* nck, char* still_active, size_t n){
	keyvalueLLU kv;
	bheapLLU* heap=constructLLU(n);
	for (size_t i=0;i<n;i++){
    if (still_active[i] != 2) {
		  kv.key=i;
		  kv.value=nck[i];
		  insertLLU(heap,kv);
    }
	}
	return heap;
}

void freeheapLLU(bheapLLU *heap){
	free(heap->pt);
	free(heap->kv);
	free(heap);
}



#endif
