#ifndef _KINTERSECT_
#define _KINTERSECT_

#pragma once

#include <math.h>

#include "ligra/bucket.h"
#include "ligra/edge_map_reduce.h"
#include "ligra/ligra.h"
#include "ligra/pbbslib/dyn_arr.h"

#include "external/simdinter/include/common.h"
#include "external/simdinter/include/intersection.h"
#include "external/graphsetinter/src/set_operation.hpp"

#define INDUCED_STACK_THR 1000

// TODO retry using lambdas for intersects
// TODO make intersects more modularized -- have some kind of wrapper that generates the vtx_seq and everything, and just have the intersects
// actually do the intersect, with the save and out_ptr options

// have a wrapper where -- if you insert a pointer, then you take responsibility for that pointer and all you get out
// is the size (and an empty sequence)
// otherwise, if you have no pointer, the sequence gets allocated for you

inline size_t lstintersect_par(uintE* a, size_t size_a, uintE* b, size_t size_b, bool save, uintE* out) {
  auto seq_a = pbbslib::make_sequence<uintE>(a, size_a);
  auto seq_b = pbbslib::make_sequence<uintE>(b, size_b);

  if (!save) {
    auto merge_f = [&] (uintE ngh) {};
    return intersection::merge(seq_a, seq_b, merge_f);
  }

  size_t index = 0;
  auto merge_f = [&] (uintE ngh) {
    out[index] = ngh;
    index++;
  };
  intersection::merge(seq_a, seq_b, merge_f);
  return index;
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

template <class F, class Graph>
std::tuple<uintE*, size_t> lstintersect(F f, Graph& DG, uintE vtx, uintE* induced, size_t induced_size, bool save = true, uintE* out_ptr = nullptr) {
  assert (vtx < DG.n);
  auto vtx_ptr = (uintE*)(DG.get_vertex(vtx).getOutNeighbors());
  auto vtx_size = DG.get_vertex(vtx).getOutDegree();
  size_t min_size = std::min((size_t) induced_size, (size_t) vtx_size);
  if (min_size == 0) return std::make_tuple(nullptr, 0);
  bool out_ptr_flag = false;
  if (!out_ptr && (save || f.count_space_flag)) {
    out_ptr_flag = true;
    out_ptr = pbbs::new_array_no_init<uintE>(min_size);
  }

  size_t out_size = f(vtx_ptr, vtx_size, induced, induced_size, save, out_ptr);

  //if (out_ptr && out_size > 0) assert(out[0] < DG.n);
  //if (out_ptr || !save) out.to_array();
  if (out_ptr_flag && (!save || out_size == 0)) {
    pbbs::delete_array<uintE>(out_ptr, min_size);
    out_ptr = nullptr;
  } //else if (out_ptr_flag && out_size < min_size) {
    //pbbs::delete_array<uintE>(out_ptr + out_size, min_size - out_size);
  //}
  return std::make_tuple(out_ptr, out_size);
}

template <class Graph>
inline sequence<uintE> kintersect(Graph& DG, sequence<uintE> base, size_t num) {
  if (num == 1) {
    uintT deg = DG.get_vertex(base[0]).getOutDegree();
    uintE* ngh = (uintE*)(DG.get_vertex(base[0]).getOutNeighbors());
    return pbbslib::make_sequence<uintE>(ngh, deg);
  }

  // find base index with min outdeg
  auto base_idxs = sequence<size_t>::no_init(num);
  parallel_for (0,num,[&] (size_t i) { base_idxs[i] = i; });
  auto base_deg_f = [&](size_t i, size_t j) -> size_t {
    return DG.get_vertex(base[i]).getOutDegree() < DG.get_vertex(base[j]).getOutDegree() ? i : j;
  };
  size_t min_base = pbbslib::reduce(base_idxs, pbbslib::make_monoid(base_deg_f, 0));
  size_t min_deg = DG.get_vertex(base[min_base]).getOutDegree();
  auto min_seq = pbbslib::make_sequence<uintE>((uintE*)(DG.get_vertex(base[min_base]).getOutNeighbors()), min_deg);

  // set up array to mark where min seq intersects with other outneighbors in base
  auto marks = sequence<bool>(min_deg * (num - 1), false);
  // mark intersections
  parallel_for (0,num,[&] (size_t i) {
    if (i != min_base) {
    size_t idx = i > min_base ? i - 1 : i;
    auto merge_f = [&] (uintE ngh) {
      size_t j = pbbslib::binary_search(min_seq, ngh, std::less<size_t>());
      marks[min_deg * idx + j] = true;
    };
    auto a_seq = pbbslib::make_sequence<uintE>((uintE*)(DG.get_vertex(base[i]).getOutNeighbors()), DG.get_vertex(base[i]).getOutDegree());
    intersection::merge(min_seq, a_seq, merge_f);
    }
  });

  // merge marks to determine which intersects all (can make more efficient here)
  auto merge_seq = sequence<bool>(min_deg, true);
  parallel_for (0,min_deg,[&] (size_t i) {
    for (size_t j=0; j < num - 1; ++j) {
      merge_seq[i] = merge_seq[i] && marks[min_deg * j + i];
    }
  });
  // pack out intersected vertices
  auto out = pbbslib::pack(min_seq, merge_seq);
  return out;
  // want to intersect V[base[i]] outneighbors
  // figure out which one represents the smallest
  // intersect each against the smallest
  // any element in the smallest that has num-1 marks is in our merged list

}




struct InducedSpace_dyn {
  size_t num_induced;
  uintE* induced;
  bool protected_flag;

  InducedSpace_dyn() : num_induced(0), induced(nullptr), protected_flag(false) {}
  InducedSpace_dyn(uintE* _induced, size_t _num_induced) : induced(_induced), num_induced(_num_induced), protected_flag(true) {}

  void alloc_induced(size_t size) {
    if (!induced) induced = pbbs::new_array_no_init<uintE>(size);
  }

  void set_num_induced(size_t size) {num_induced = size;}

  void del() {
    if (!protected_flag && induced) {pbbs::delete_array<uintE>(induced, num_induced); induced = nullptr;}
  }

  static void init() {}
  static void finish() {}
};

struct InducedSpace_alloc {
  using induced_alloc = list_allocator<uintE[INDUCED_STACK_THR]>;
  size_t num_induced;
  uintE* induced;
  uintE(* ptr_induced)[INDUCED_STACK_THR];

  InducedSpace_alloc() : num_induced(0), induced(nullptr) {}

  void alloc_induced(size_t size) {
    if (!induced) {
      ptr_induced = induced_alloc::alloc();
      induced = *ptr_induced;
    }
  }

  void set_num_induced(size_t size) {num_induced = size;}

  void del() {
    if (induced) { induced_alloc::free(ptr_induced); induced = nullptr; }
  }

  static void init() {
    induced_alloc::init();
  }

  static void finish() {
    induced_alloc::finish();
  }
};

struct InducedSpace_stack {
  size_t num_induced;
  uintE induced[INDUCED_STACK_THR];

  InducedSpace_stack() : num_induced(0) {}

  void alloc_induced(size_t size) {}

  void set_num_induced(size_t size) {num_induced = size;}

  void del() {}

  static void init() {}
  static void finish() {}
};

// induced_space must have: num_induced, .clear(), induced (array)
template <class Graph, class I, class F, class IN>
size_t lstintersect_induced(Graph& DG, size_t k_idx, size_t k, size_t i, I& induced_space, F intersect_op_type, 
  sequence<uintE>& base, bool count_only, bool to_save, IN& new_induced_space) {
  if (!count_only) base[k_idx] = induced_space.induced[i];
  uintE vtx = induced_space.induced[i];
  //assert (vtx < DG.n);

  auto vtx_ptr = (uintE*)(DG.get_vertex(vtx).getOutNeighbors());
  auto vtx_size = DG.get_vertex(vtx).getOutDegree();
  size_t min_size = std::min((size_t) induced_space.num_induced, (size_t) vtx_size);
  if (min_size < k - k_idx) return 0;

  bool out_ptr_flag = false;
  if (!new_induced_space.induced && (to_save || intersect_op_type.count_space_flag)) {
    out_ptr_flag = true;
    new_induced_space.alloc_induced(min_size);
  }
  auto out_ptr = new_induced_space.induced;

  size_t out_size = intersect_op_type(vtx_ptr, vtx_size, induced_space.induced, induced_space.num_induced, to_save, out_ptr);

  if (out_ptr_flag && (!to_save || out_size == 0)) {
    new_induced_space.del();
  }
  new_induced_space.num_induced = out_size;
  return out_size;
}

/*template <class F>
struct lstintersect_induced_struct {
  size_t k;
  F intersect_op_type;
  bool count_only;

  lstintersect_induced_struct (size_t  _k, F _intersect_op_type, bool _count_only) : k(_k), intersect_op_type(_intersect_op_type), count_only(_count_only) {}
  
  template <class Graph, class I, class IN>
  size_t operator()(Graph& DG, size_t k_idx, size_t i, I& induced_space, sequence<uintE>& base, bool to_save, IN& new_induced_space) const {
    return lstintersect_induced(DG, k_idx, k, i, induced_space, intersect_op_type, base, count_only, to_save, new_induced_space);
  }
};*/

#endif