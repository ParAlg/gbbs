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

template <class A, class S, class O>
inline size_t lstintersect_par(A& vtx_seq, S& induced, bool save, O& out) {
  size_t index = 0;
  if (!save) {
    auto merge_f = [&] (uintE ngh) {};
    return intersection::merge(vtx_seq, induced, merge_f);
  }

  auto merge_f = [&] (uintE ngh) {
    out[index] = ngh;
    index++;
  };
  intersection::merge(vtx_seq, induced, merge_f);
  out.shrink(index);
  return out.size();
}

struct lstintersect_par_struct {
  template <class A, class S, class O>
  size_t operator()(A& vtx_seq, S& induced, bool save, O& out) const {
    return lstintersect_par(vtx_seq, induced, save, out);
  }
};

// TODO radix sort in place
template <class A, class S, class O>
inline size_t lstintersect_vec(A& vtx_seq, S& induced, bool save, O& out) {
  SIMDCompressionLib::intersectionfunction inter = SIMDCompressionLib::IntersectionFactory::getFromName("simd");
  size_t out_size = inter(induced.begin(), induced.size(), vtx_seq.begin(), vtx_seq.size(), out.begin());
  out.shrink(out_size);
  return out_size;
}

struct lstintersect_vec_struct {
  template <class A, class S, class O>
  size_t operator()(A& vtx_seq, S& induced, bool save, O& out) const {
    return lstintersect_vec(vtx_seq, induced, save, out);
  }
};


template <class A, class S, class O>
inline size_t lstintersect_set(A& vtx_seq, S& induced, bool save, O& out) {
  if (!save) {
#if SIMD_STATE == 2
    return (size_t) intersect_scalar2x_count((int*) induced.begin(), (int) induced.size(), (int*) vtx_seq.begin(), (int) vtx_seq.size());
#elif SIMD_STATE == 4
    return (size_t) intersect_simd4x_count((int*) induced.begin(), (int) induced.size(), (int*) vtx_seq.begin(), (int) vtx_seq.size());
#else
    return (size_t) intersect_count((int*) induced.begin(), (int) induced.size(), (int*) vtx_seq.begin(), (int) vtx_seq.size());
#endif
  }
#if SIMD_STATE == 2
  int out_size = intersect_scalar2x((int*) induced.begin(), (int) induced.size(), (int*) vtx_seq.begin(), (int) vtx_seq.size(), (int*) out.begin());
#elif SIMD_STATE == 4
  int out_size = intersect_simd4x((int*) induced.begin(), (int) induced.size(), (int*) vtx_seq.begin(), (int) vtx_seq.size(), (int*) out.begin());
#else
  int out_size = intersect((int*) induced.begin(), (int) induced.size(), (int*) vtx_seq.begin(), (int) vtx_seq.size(), (int*) out.begin());
#endif
  out.shrink(out_size);
  return (size_t) out_size;
}

struct lstintersect_set_struct {
  template <class A, class S, class O>
  size_t operator()(A& vtx_seq, S& induced, bool save, O& out) const {
    return lstintersect_set(vtx_seq, induced, save, out);
  }
};

template <class F, class Graph, class S>
std::tuple<sequence<uintE>, size_t> lstintersect(F f, Graph& DG, uintE vtx, S& induced, bool save = true, uintE* out_ptr = nullptr) {
  assert (vtx < DG.n);
  auto vtx_seq = pbbslib::make_sequence<uintE>((uintE*)(DG.get_vertex(vtx).getOutNeighbors()), DG.get_vertex(vtx).getOutDegree());
  size_t min_size = std::min(induced.size(), vtx_seq.size());

  sequence<uintE> out;
  if (out_ptr) out = pbbslib::make_sequence<uintE>(out_ptr, min_size);
  else out = sequence<uintE>::no_init(min_size);

  size_t out_size = f(vtx_seq, induced, save, out);
  /*SIMDCompressionLib::intersectionfunction inter = SIMDCompressionLib::IntersectionFactory::getFromName("simd");
  size_t out_size = inter(induced.begin(), induced.size(), vtx_seq.begin(), vtx_seq.size(), out.begin());
  out.shrink(out_size);*/

  //if (out_ptr || !save) out.to_array();
  return std::make_tuple(out, out_size);
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

#endif