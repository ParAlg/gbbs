#ifndef _KINTERSECT_
#define _KINTERSECT_

#pragma once

#include <math.h>

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
      }
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
  HybridSpace_lw () {}

  void alloc(size_t max_induced, size_t k, size_t n) {
    if (induced == nullptr && k > 2) induced = (uintE*) malloc(sizeof(uintE)*k*max_induced);
    if (induced_degs == nullptr) induced_degs = (uintE*) malloc(sizeof(uintE)*max_induced);
    if (labels == nullptr && k > 2) labels = (char*) calloc(max_induced, sizeof(char));
    if (induced_edges == nullptr && k > 2) induced_edges = (uintE*) malloc(sizeof(uintE)*max_induced*max_induced);
    if (num_induced == nullptr && k > 2) num_induced = (uintE*) malloc(sizeof(uintE)*k);
    if (old_labels == nullptr) old_labels = (uintE*) calloc(n, sizeof(uintE));
  }

  template <class Graph>
  void setup(Graph& DG, size_t k, size_t i) {
    nn = DG.get_vertex(i).getOutDegree();
    auto induced_g = DG.get_vertex(i).getOutNeighbors(); //((uintE*)(DG.get_vertex(i).getOutNeighbors()));
    for (size_t j=0; j < nn; j++) { induced_degs[j] = 0; }
  
    if (k > 2) {
      num_induced[0] = nn;
      for (size_t  j=0; j < nn; j++) { induced[j] = j; }
    }

    for (size_t o=0; o < nn; o++) { old_labels[std::get<0>(induced_g[o])] = o + 1; }


    for (size_t j=0; j < nn; j++) {
      uintE v = std::get<0>(induced_g[j]);
      auto v_nbhrs = DG.get_vertex(v).getOutNeighbors();
      size_t v_deg = DG.get_vertex(v).getOutDegree();
      // intersect v_nbhrs from 0 to v_deg with induced_g from 0 to num_induced[0]
      // store result in induced_edges[j*nn]
      // store size in induced_degs[j]
      for (size_t l=0; l < v_deg; l++) {
        if (old_labels[std::get<0>(v_nbhrs[l])] > 0) {
          if (k > 2) induced_edges[j*nn + induced_degs[j]] = old_labels[std::get<0>(v_nbhrs[l])] - 1;
          induced_degs[j]++;
        }
      }
    }

    for (size_t o=0; o < nn; o++) { old_labels[std::get<0>(induced_g[o])] = 0; }

    auto deg_seq = pbbslib::make_sequence(induced_degs, nn);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  static void init(){}
  static void finish(){}

  void del() {
    if (labels) {free(labels); labels=nullptr;}
    if (induced) {free(induced); induced=nullptr;}
    if (induced_edges) {free(induced_edges); induced_edges=nullptr;}
    if (induced_degs) {free(induced_degs); induced_degs=nullptr;}
    if (num_induced) {free(num_induced); num_induced=nullptr;}
    if (old_labels) {free(old_labels); old_labels=nullptr;}
  }

  ~HybridSpace_lw() { del(); }

};







#endif