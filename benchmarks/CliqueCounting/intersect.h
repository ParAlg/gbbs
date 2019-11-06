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

#define INDUCED_STACK_THR 500

// TODO retry using lambdas for intersects
// TODO make intersects more modularized -- have some kind of wrapper that generates the vtx_seq and everything, and just have the intersects
// actually do the intersect, with the save and out_ptr options

// have a wrapper where -- if you insert a pointer, then you take responsibility for that pointer and all you get out
// is the size (and an empty sequence)
// otherwise, if you have no pointer, the sequence gets allocated for you


inline size_t lstintersect_simple(uintE* a, size_t size_a, uintE* b, size_t size_b, bool save, uintE* out) {
  auto seq_b = pbbslib::make_sequence(b, size_b);
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
  bool full_flag = false;
  size_t num_edges = 0;
  bool pruned = false;
  size_t running_sum = 0;

  InducedSpace_dyn() : num_induced(0), induced(nullptr), protected_flag(false) {}

  template <class Graph>
  InducedSpace_dyn(Graph& DG, size_t k, size_t i) : protected_flag(true) {
    induced = (uintE*)(DG.get_vertex(i).getOutNeighbors());
    num_induced = DG.get_vertex(i).getOutDegree();
    running_sum = num_induced;
  }

  template <class I>
  void alloc_induced(size_t size, I prev) {
    if (!induced) induced = pbbs::new_array_no_init<uintE>(size);
  }

  void del() {
    timer t; t.start();
    if (!protected_flag && induced) {pbbs::delete_array<uintE>(induced, num_induced);}
    induced = nullptr;
    double time = t.stop();
  }

  //~InducedSpace_dyn() {del();}

  static void init() {}
  static void finish() {}
};

struct InducedSpace_alloc {
  using induced_alloc = list_allocator<uintE[INDUCED_STACK_THR]>;
  size_t num_induced;
  uintE* induced;
  uintE(* ptr_induced)[INDUCED_STACK_THR];
  bool full_flag = false;
  size_t num_edges = 0;
  bool pruned = false;
  size_t running_sum = 0;

  InducedSpace_alloc() : num_induced(0), induced(nullptr) {}

  template <class I>
  void alloc_induced(size_t size, I prev) {
    if (!induced) {
      ptr_induced = induced_alloc::alloc();
      induced = *ptr_induced;
    }
  }

  void del() {
    timer t; t.start();
    if (induced) { induced_alloc::free(ptr_induced); }
    induced = nullptr; 
    double time = t.stop();
  }

  //~InducedSpace_alloc() {del();}

  static void init() {
    timer t; t.start();
    induced_alloc::init();
    double time = t.stop();
  }

  static void finish() {
    timer t; t.start();
    induced_alloc::finish();
    double time = t.stop();
  }
};

struct InducedSpace_stack {
  size_t num_induced;
  uintE* induced = nullptr;
  uintE induced_stack[INDUCED_STACK_THR];
  bool full_flag = false;
  size_t num_edges = 0;
  size_t running_sum = 0;

  InducedSpace_stack() : num_induced(0) {}

  template <class I>
  void alloc_induced(size_t size, I prv) {induced =  (uintE*) induced_stack;}

  void del() {induced = nullptr;}

  //~InducedSpace_stack() {}

  static void init() {}
  static void finish() {}
};

struct InducedSpace_stack_setup {
  size_t num_induced;
  uintE* induced = nullptr;
  uintE induced_stack[INDUCED_STACK_THR];
  bool full_flag = false;
  size_t num_edges = 0;
  size_t running_sum = 0;

  InducedSpace_stack_setup() : num_induced(0) {}

  template <class Graph>
  InducedSpace_stack_setup(Graph& DG, size_t k, size_t i) {
    num_induced = DG.get_vertex(i).getOutDegree();
    auto tmp_induced = (uintE*)(DG.get_vertex(i).getOutNeighbors());
    parallel_for (0, num_induced, [&] (size_t j) {
      induced_stack[j] = tmp_induced[j];
    });
    induced =  (uintE*) induced_stack;
    running_sum = num_induced;
  }

  template <class I>
  void alloc_induced(size_t size, I prev) {}

  void del() { induced = nullptr; }

  //~InducedSpace_stack_setup() {del();}

  static void init() {}
  static void finish() {}
};

struct InducedSpace_dyn_setup {
  size_t num_induced;
  uintE* induced = nullptr;
  bool full_flag = false;
  size_t num_edges = 0;
  size_t running_sum = 0;

  InducedSpace_dyn_setup() : num_induced(0) {}

  template <class Graph>
  InducedSpace_dyn_setup(Graph& DG, size_t k, size_t i) {
    num_induced = DG.get_vertex(i).getOutDegree();
    induced = pbbs::new_array_no_init<uintE>(k*num_induced);
    auto tmp_induced = (uintE*)(DG.get_vertex(i).getOutNeighbors());
    parallel_for (0, num_induced, [&] (size_t j) {
      induced[j] = tmp_induced[j];
    });
    running_sum = num_induced;
  }

  template <class I>
  void alloc_induced(size_t size, I prev) {}

  void del() { if (induced) pbbs::free_array(induced); induced = nullptr; }

  //~InducedSpace_dyn_setup() {del();}

  static void init() {}
  static void finish() {}
};

struct InducedSpace_rec {
  size_t num_induced;
  uintE* induced = nullptr;
  bool full_flag = false;
  size_t num_edges = 0;
  size_t running_sum = 0;

  InducedSpace_rec() : num_induced(0) {}

  template <class I>
  void alloc_induced(size_t size, I prev) { induced = prev.induced + prev.num_induced; }

  void del() { induced = nullptr; }

  //~InducedSpace_rec() {del();}

  static void init() {}
  static void finish() {}
};


// induced_space must have: num_induced, .clear(), induced (array)
template <class Graph, class I, class F, class IN>
inline size_t lstintersect_induced(Graph& DG, size_t k_idx, size_t k, size_t i, I& induced_space, F intersect_op_type, 
  sequence<uintE>& base, bool count_only, bool to_save, IN& new_induced_space) {
  if (!count_only) base[k_idx] = induced_space.induced[i];
  uintE vtx = induced_space.induced[i];
  //assert (vtx < DG.n);

  auto vtx_ptr = (uintE*)(DG.get_vertex(vtx).getOutNeighbors());
  auto vtx_size = DG.get_vertex(vtx).getOutDegree();
  size_t min_size = std::min((size_t) induced_space.num_induced, (size_t) vtx_size);
  if (min_size < k - k_idx) {new_induced_space.running_sum = induced_space.running_sum; return 0;}

  bool out_ptr_flag = false;
  if (!new_induced_space.induced && (to_save || intersect_op_type.count_space_flag)) {
    out_ptr_flag = true;
    new_induced_space.alloc_induced(min_size, induced_space);
  }
  auto out_ptr = new_induced_space.induced;
  new_induced_space.num_induced = min_size;

  size_t out_size = intersect_op_type(vtx_ptr, vtx_size, induced_space.induced, induced_space.num_induced, to_save, out_ptr);

  if (out_ptr_flag && (!to_save || out_size == 0)) {
    new_induced_space.del();
  }
  new_induced_space.num_induced = out_size;
  new_induced_space.running_sum = induced_space.running_sum + new_induced_space.num_induced;
  return out_size;
}

template <class F>
struct lstintersect_induced_struct {
  size_t k;
  F intersect_op_type;
  bool count_only;

  lstintersect_induced_struct (size_t  _k, F _intersect_op_type, bool _count_only) : k(_k-1), intersect_op_type(_intersect_op_type), count_only(_count_only) {}
  
  template <class Graph, class I, class IN>
  size_t operator()(Graph& DG, size_t k_idx, size_t i, I& induced_space, sequence<uintE>& base, bool to_save, IN& new_induced_space) const {
    return lstintersect_induced(DG, k_idx, k, i, induced_space, intersect_op_type, base, count_only, to_save, new_induced_space);
  }
};

struct lstintersect_induced_struct2 {
  template <class Graph, class I, class F, class IN>
  size_t operator()(Graph& DG, size_t k_idx, size_t k, size_t i, I& induced_space, F intersect_op_type, sequence<uintE>& base, bool count_only, bool to_save, IN& new_induced_space)const {
    return lstintersect_induced(DG, k_idx, k, i, induced_space, intersect_op_type, base, count_only, to_save, new_induced_space);
  }
};

// for each node u in induced, and each node v adj to u, if label of v is l, then set  to l-1 and add v to our new level
// then, for each v that we've added to our new subgraph, for each w that was adjacent to it in G, check if w was labeled with l-1; if so, add  to adj list of v
// rearrange and keep track of an endpoint

// space must have: getDegree(i), getNeighbors(i), induced, num_induced, prune(induced_space, min_deg), alloc_induced(size)
template <class Graph, class I, class F, class IN>
inline size_t lstintersect_orig(Graph& DG, size_t k_idx, size_t k, size_t i, I& induced_space, F intersect_op_type, 
  sequence<uintE>& base, bool count_only, bool to_save, IN& new_induced_space) {
  if (!count_only) base[k_idx] = induced_space.induced[i];
  if (induced_space.getDegree(i) < k - k_idx) {new_induced_space.running_sum = induced_space.running_sum; return 0;}

  new_induced_space.prune(DG, i, induced_space, (k-k_idx > 2), k_idx);
  return new_induced_space.num_induced;
}

struct lstintersect_orig_struct {
  template <class Graph, class I, class F, class IN>
  size_t operator()(Graph& DG, size_t k_idx, size_t k, size_t i, I& induced_space, F intersect_op_type, sequence<uintE>& base, bool count_only, bool to_save, IN& new_induced_space)const {
    return lstintersect_orig(DG, k_idx, k, i, induced_space, intersect_op_type, base, count_only, to_save, new_induced_space);
  }
};

struct FullSpace_orig {
  size_t num_induced = 0;
  uintE* induced = nullptr;
  bool protected_flag = false;
  bool full_flag = true;
  bool orig_flag = false;
  size_t num_edges = 0;
  uintE* induced_edges = nullptr;
  size_t* induced_degs = nullptr;
  uintE* labels = nullptr;
  size_t nn = 0;
  size_t step = 0;
  size_t running_sum = 0;

  FullSpace_orig() {}

  size_t getDegree(size_t i) { return induced_degs[induced[i]]; }

  template <class Graph>
  FullSpace_orig(Graph& DG, size_t k, size_t i) : protected_flag(false), orig_flag(true) {
    num_induced = DG.get_vertex(i).getOutDegree();
    running_sum = num_induced;
    if (num_induced == 0) return;
    nn = num_induced;
    induced = pbbs::new_array_no_init<uintE>(num_induced);
    uintE* induced_g = ((uintE*)(DG.get_vertex(i).getOutNeighbors()));
    parallel_for(0, num_induced, [&] (size_t j) { induced[j] = j; });
    induced_degs = pbbs::new_array_no_init<size_t>(nn);
    parallel_for(0, nn, [&] (size_t j) { induced_degs[j] = 0; });
    labels = pbbs::new_array_no_init<uintE>(nn);
    parallel_for(0, nn, [&] (size_t j) { labels[j] = 0; });

    /*auto idxs = sequence<size_t>::no_init(num_induced);
    parallel_for (0, num_induced, [&] (size_t j) { idxs[j] = DG.get_vertex(induced_g[j]).getOutDegree(); });
    auto base_deg_f = [&](size_t l, size_t j) -> size_t {
      return idxs[l] > idxs[j] ? idxs[l] : idxs[j];
    };
    step = pbbslib::reduce(idxs, pbbslib::make_monoid(base_deg_f, 0));*/
    /*step = 0;
    for (size_t j=0; j < num_induced; j++) {
      if (DG.get_vertex(induced_g[j]).getOutDegree() > step) step = DG.get_vertex(induced_g[j]).getOutDegree();
      if (step > num_induced) {
        step = num_induced;
        break;
      }
    }*/
    step = num_induced;
    //auto intersect_op_type = lstintersect_vec_struct{};

    induced_edges = pbbs::new_array_no_init<uintE>(nn*step);
    parallel_for(0, num_induced, [&] (size_t j) {
      uintE v = induced_g[j];
      uintE* v_nbhrs = (uintE*)(DG.get_vertex(v).getOutNeighbors());
      size_t v_deg = DG.get_vertex(v).getOutDegree();
      for (size_t l=0; l < v_deg; l++) {
        for (size_t o=0; o < num_induced; o++) {
          if (v_nbhrs[l] == induced_g[o]) {
            induced_edges[j*step + induced_degs[j]] = o;
            induced_degs[j]++;
            break;
          }
        }
      }
      /*induced_degs[j] = intersect_op_type(v_nbhrs, v_deg, induced_g, num_induced, true, induced_edges+j*step);
      // map vert we put into induced_edges to induced
      parallel_for(0, induced_degs[j], [&] (size_t l) {
        for (size_t o=0; o<num_induced; o++) {
          if (induced_edges[j*step + l] == induced_g[o]) {
            induced_edges[j*step + l] = o;
            break;
          }
        }
      });*/
    });
    auto deg_seq = pbbslib::make_sequence(induced_degs, nn);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  template<class Graph, class O>
  void prune(Graph& DG, size_t i, O& orig, bool to_save, size_t k_idx) {
    protected_flag = true;
    orig_flag = false;
    uintE idx = orig.induced[i];
    num_induced = orig.induced_degs[idx];
    if (num_induced == 0) {running_sum = orig.running_sum; return;}
    induced = orig.induced_edges + (idx*orig.step);
    nn = orig.nn;
    step = orig.step;
    induced_edges = orig.induced_edges;
    labels = orig.labels;
    parallel_for(0, num_induced, [&] (size_t j){
      labels[induced[j]] = k_idx;
    });
    
    
    induced_degs = pbbs::new_array_no_init<size_t>(nn);
    parallel_for(0, nn, [&] (size_t j) { induced_degs[j] = 0; });
    
    parallel_for(0, num_induced, [&] (size_t j) {
      uintE v_idx = induced[j];
      uintE v_deg = orig.induced_degs[v_idx];
      uintE* v_edges = induced_edges + v_idx*step;
      size_t end = v_deg;
      for (size_t l=0; l < end; l++) {
        if (labels[v_edges[l]] == k_idx) induced_degs[v_idx]++;
        else if (to_save){
          auto tmp = v_edges[l];
          v_edges[l--] = v_edges[--end];
          v_edges[end] = tmp;
        }
      }
    });
    
    auto deg_seq = pbbslib::make_sequence(induced_degs, nn);
    num_edges = pbbslib::reduce_add(deg_seq);
    running_sum = orig.running_sum + num_induced;
  }

  static void init(){}
  static void finish(){}

  void del() {
    if (orig_flag && labels) {pbbs::delete_array<uintE>(labels, nn);}
    else if (labels && induced) {
      parallel_for(0, num_induced, [&] (size_t j){
        labels[induced[j]] = 0;
      });
    }
    if (!protected_flag && induced) {pbbs::delete_array<uintE>(induced, num_induced);}
    induced = nullptr;
    if (orig_flag && induced_edges) {pbbs::delete_array<uintE>(induced_edges, nn*step);}
    induced_edges = nullptr;
    labels = nullptr;
    if (induced_degs) {pbbs::delete_array<size_t>(induced_degs, nn);}
    induced_degs = nullptr;
  }

};








  struct hash_struct {
    inline size_t operator () (const uintE& t) {
      return pbbslib::hash64_2((size_t) t);
    }
  };

struct FullSpace_bool_dyn {
  size_t num_induced;
  uintE* induced;
  sparse_table<uintE, uintE, hash_struct> induced_table = sparse_table<uintE, uintE, hash_struct>(0, std::make_tuple((uintE) UINT_E_MAX, (uintE) UINT_E_MAX), hash_struct());
  bool protected_flag;
  bool full_flag = true;
  bool del_flag = false;
  size_t num_edges = 0;
  uintE** induced_edges;
  uintE* induced_degs;
  size_t nn = 0;

  FullSpace_bool_dyn() : num_induced(0), induced(nullptr), protected_flag(false), induced_edges(nullptr), induced_degs(nullptr) {}

  size_t getDegree(size_t i) { return induced_degs[i]; }
  uintE* getNeighbors(size_t i) { return induced_edges[i]; }
  void hash_induced() {
    induced_table.resize_no_copy((size_t)1 << pbbslib::log2_up((size_t)(1.5 * num_induced)));
    parallel_for (0, num_induced, [&] (size_t j) {
      induced_table.insert(std::make_tuple(induced[j], (uintE) j));
    });
  }

  template <class Graph>
  FullSpace_bool_dyn(Graph& DG, size_t k, size_t i) : protected_flag(true), del_flag(true) {
    num_induced = DG.get_vertex(i).getOutDegree();
    if (num_induced == 0) return;
    induced = (uintE*)(DG.get_vertex(i).getOutNeighbors());
    hash_induced();
  
    induced_degs = pbbs::new_array_no_init<uintE>(num_induced);
    induced_edges = pbbs::new_array_no_init<uintE*>(num_induced);

    auto intersect_op_type = lstintersect_vec_struct{};

    parallel_for(0, num_induced, [&] (size_t j) {
      // intersect neighbors of each vert in induced with induced
      uintE v = induced[j];
      uintE* nbhrs = (uintE*)(DG.get_vertex(v).getOutNeighbors());
      size_t deg = DG.get_vertex(v).getOutDegree();
      induced_edges[j] = pbbs::new_array_no_init<uintE>(std::min((size_t) deg, (size_t) num_induced));
      induced_degs[j] = intersect_op_type(nbhrs, deg, induced, num_induced, true, induced_edges[j]);
      if (induced_degs[j] == 0) {
        pbbs::delete_array<uintE>(induced_edges[j], std::min((size_t) deg, (size_t) num_induced));
        induced_edges[j] = nullptr;
      }
    });
    auto deg_seq = pbbslib::make_sequence(induced_degs, num_induced);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  template<class Graph, class O, class F>
  void prune(Graph& DG, size_t idx, O& orig, bool to_save, F intersect_op_type) {
    protected_flag = true;
    num_induced = orig.getDegree(idx);
    if (num_induced == 0) return;
    induced = orig.getNeighbors(idx);
    hash_induced();

    induced_degs = pbbs::new_array_no_init<uintE>(num_induced);
    if (to_save) induced_edges = pbbs::new_array_no_init<uintE*>(num_induced);

    parallel_for(0, num_induced, [&] (size_t j) {
    //for (size_t j=0; j < num_induced; j++) {
      uintE v = induced[j];
      //size_t find_idx = pbbs::binary_search(orig_induced_seq, v, lte);
      size_t find_idx = orig.induced_table.find(v, UINT_E_MAX);
      /*for (find_idx=0; find_idx < orig.num_induced; find_idx++) {
          if (orig.induced[find_idx] == v) break;
        }*/
      uintE v_deg = orig.induced_degs[find_idx];
      uintE* v_edges = orig.induced_edges[find_idx];
      induced_degs[j] = 0;

      if (to_save) induced_edges[j] = v_edges;
      size_t end = v_deg;
      // set induced_degs here, and move to front if needed
      for (size_t l=0; l < end; l++) {
        /*size_t xx  = 0;
        for (xx=0; xx < num_induced; xx++) {
          if (induced[xx] == v_edges[l]) break;
        }*/
        bool is_in = induced_table.contains(v_edges[l]);
        if (is_in) induced_degs[j]++;
        else if (to_save) {
          auto tmp = v_edges[l];
          end--;
          v_edges[l] = v_edges[end];
          l--;
          v_edges[end] = tmp;
        }
        // if v_edges[l] is in induced, add 1 to induced degs
        // if not and if to_save, swap with v_edges[end]
      }
      if (to_save && induced_degs[j] == 0) induced_edges[j] = nullptr;
    });
    auto deg_seq = pbbslib::make_sequence(induced_degs, num_induced);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  void alloc_induced(size_t size) {
    if (!induced) induced = pbbs::new_array_no_init<uintE>(size);
  }

  void del() {
    if (!protected_flag && induced) {pbbs::delete_array<uintE>(induced, num_induced);}
    induced = nullptr;
    if (del_flag && induced_edges && induced_degs) {
      for (size_t i=0; i < num_induced; ++i) {
        if (induced_edges[i]) pbbs::delete_array<uintE>(induced_edges[i], induced_degs[i]);
        induced_edges[i] = nullptr;
      }
    }
    if (induced_edges) {pbbs::delete_array<uintE*>(induced_edges, num_induced);}
    induced_edges = nullptr;
    if (induced_degs) {pbbs::delete_array<uintE>(induced_degs, num_induced);}
    induced_degs = nullptr;
    induced_table.del();
  }

  ~FullSpace_bool_dyn() {del();}

  static void init() {}
  static void finish() {}
};










struct FullSpace_csv_dyn {
  size_t num_induced;
  uintE* induced;
  bool protected_flag;
  bool full_flag = true;
  size_t num_edges = 0;
  uintE** induced_edges;
  uintE* induced_degs;
  size_t nn = 0;

  FullSpace_csv_dyn() : num_induced(0), induced(nullptr), protected_flag(false), induced_edges(nullptr), induced_degs(nullptr) {}

  size_t getDegree(size_t i) { return induced_degs[i]; }
  uintE* getNeighbors(size_t i) { return induced_edges[i]; }

  template <class Graph>
  FullSpace_csv_dyn(Graph& DG, size_t k, size_t i) : protected_flag(true) {
    induced = (uintE*)(DG.get_vertex(i).getOutNeighbors());
    num_induced = DG.get_vertex(i).getOutDegree();
    if (num_induced == 0) return;
  
    induced_degs = pbbs::new_array_no_init<uintE>(num_induced);
    induced_edges = pbbs::new_array_no_init<uintE*>(num_induced);

    auto intersect_op_type = lstintersect_vec_struct{};

    parallel_for(0, num_induced, [&] (size_t j) {
    //for (size_t j = 0; j < num_induced; j++) {
      // intersect neighbors of each vert in induced with induced
      uintE v = induced[j];
      uintE* nbhrs = (uintE*)(DG.get_vertex(v).getOutNeighbors());
      size_t deg = DG.get_vertex(v).getOutDegree();
      induced_edges[j] = pbbs::new_array_no_init<uintE>(std::min((size_t) deg, (size_t) num_induced));
      induced_degs[j] = intersect_op_type(nbhrs, deg, induced, num_induced, true, induced_edges[j]);
      if (induced_degs[j] == 0) {
        pbbs::delete_array<uintE>(induced_edges[j], std::min((size_t) deg, (size_t) num_induced));
        induced_edges[j] = nullptr;
      }
    });
    auto deg_seq = pbbslib::make_sequence(induced_degs, num_induced);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  template<class Graph, class O, class F>
  void prune(Graph& DG, size_t idx, O& orig, bool to_save, F intersect_op_type) {
    protected_flag = true;
    induced = orig.getNeighbors(idx);
    num_induced = orig.getDegree(idx);
    if (num_induced == 0 || orig.num_edges == 0 || !induced) return;
    induced_degs = pbbs::new_array_no_init<uintE>(num_induced);
    if (to_save || intersect_op_type.count_space_flag) {
      induced_edges = pbbs::new_array_no_init<uintE*>(num_induced);
    }
    auto lte = [] (const uintE& l, const uintE& r) { return l < r; };
    auto orig_induced_seq = pbbslib::make_sequence<uintE>(orig.induced, orig.num_induced);

    parallel_for(0, num_induced, [&] (size_t j) {
    //for (size_t j=0; j < num_induced; j++) {
      uintE v = induced[j];
      size_t find_idx = pbbs::binary_search(orig_induced_seq, v, lte);
      uintE v_deg = orig.induced_degs[find_idx];
      uintE* v_edges = orig.induced_edges[find_idx];

      if (to_save || intersect_op_type.count_space_flag) {
        induced_edges[j] = pbbs::new_array_no_init<uintE>(std::min((size_t) v_deg, (size_t) num_induced));
        induced_degs[j] = intersect_op_type(v_edges, v_deg, induced, num_induced, true, induced_edges[j]);
        if (induced_degs[j] == 0) {
          pbbs::delete_array<uintE>(induced_edges[j], std::min((size_t) v_deg, (size_t) num_induced));
          induced_edges[j] = nullptr;
        }
      } else induced_degs[j] = intersect_op_type(v_edges, v_deg, induced, num_induced, false, nullptr);
    });
    auto deg_seq = pbbslib::make_sequence(induced_degs, num_induced);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  void alloc_induced(size_t size) {
    if (!induced) induced = pbbs::new_array_no_init<uintE>(size);
  }

  void del() {
    if (!protected_flag && induced) {pbbs::delete_array<uintE>(induced, num_induced);}
    induced = nullptr;
    if (induced_edges && induced_degs) {
      for (size_t i=0; i < num_induced; ++i) {
        if (induced_edges[i]) pbbs::delete_array<uintE>(induced_edges[i], induced_degs[i]);
        induced_edges[i] = nullptr;
      }
      pbbs::delete_array<uintE*>(induced_edges, num_induced);
    }
    induced_edges = nullptr;
    if (induced_degs) {pbbs::delete_array<uintE>(induced_degs, num_induced);}
    induced_degs = nullptr;
  }

  ~FullSpace_csv_dyn() {del();}

  static void init() {}
  static void finish() {}
};



// space must have: getDegree(i), getNeighbors(i), induced, num_induced, prune(induced_space, min_deg), alloc_induced(size)
template <class Graph, class I, class F, class IN>
inline size_t lstintersect_full(Graph& DG, size_t k_idx, size_t k, size_t i, I& induced_space, F intersect_op_type, 
  sequence<uintE>& base, bool count_only, bool to_save, IN& new_induced_space) {
  if (!count_only) base[k_idx] = induced_space.induced[i];
  if (induced_space.getDegree(i) < k - k_idx) return 0;

  new_induced_space.prune(DG, i, induced_space, (k-k_idx > 2), intersect_op_type);
  return new_induced_space.num_induced;
}





/*

  struct hash_struct {
    inline size_t operator () (const uintE& t) {
      return pbbslib::hash64_2((size_t) t);
    }
  };

struct FullSpace_csv_hashhash_dyn {
  size_t num_induced = 0;
  uintE* induced = nullptr;
  sparse_table<uintE, uintE, hash_struct> induced_table = sparse_table(0, std::make_tuple((uintE) UINT_E_MAX, (uintE) UINT_E_MAX), hash_struct());
  bool protected_flag = false;
  bool full_flag = true;
  size_t num_edges = 0;
  uintE** induced_edges = nullptr;
  uintE* induced_degs = nullptr;
  size_t nn = 0;

  FullSpace_csv_hashhash_dyn() {}
  size_t getDegree(size_t i) { return induced_degs[i]; }
  uintE* getNeighbors(size_t i) { return induced_edges[i]; }
  void hash_induced() {
    induced_table.resize_no_copy(num_induced);
    parallel_for (0, num_induced, [&] (size_t j) {
      induced_table.insert(std::make_tuple(induced[j], (uintE) j));
    });
  }

  template <class Graph>
  FullSpace_csv_hashhash_dyn(Graph& DG, size_t k, size_t i) : protected_flag(true) {
    induced = (uintE*)(DG.get_vertex(i).getOutNeighbors());
    num_induced = DG.get_vertex(i).getOutDegree();
    if (num_induced == 0) return;
    nn = DG.n;
    hash_induced();
  
    induced_degs = pbbs::new_array_no_init<uintE>(num_induced);
    induced_edges = pbbs::new_array_no_init<uintE*>(num_induced);

    parallel_for (0, num_induced, [&] (size_t j) {
      // intersect neighbors of each vert in induced with induced
      uintE v = induced[j];
      uintE* nbhrs = (uintE*)(DG.get_vertex(v).getOutNeighbors());
      size_t deg = DG.get_vertex(v).getOutDegree();
      induced_degs[j] = 0;
      induced_edges[j] = nullptr;
      if (deg > 0) {
        induced_edges[j] = pbbs::new_array_no_init<uintE>(std::min((size_t) deg, (size_t) num_induced));
        for (size_t l = 0; l < deg; l++) {
          if (induced_table.contains(nbhrs[l])) {
            induced_edges[j][induced_degs[j]] = nbhrs[l];
            induced_degs[j]++;
          }
        }
      //induced_degs[j] = intersect_op_type(nbhrs, deg, induced, num_induced, true, induced_edges[j]);
      if (induced_degs[j] == 0) {
        pbbs::delete_array<uintE>(induced_edges[j], std::min((size_t) deg, (size_t) num_induced));
        induced_edges[j] = nullptr;
      }
      }
    });

    auto deg_seq = pbbslib::make_sequence(induced_degs, num_induced);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  template<class O, class F>
  void prune(O& orig, bool to_save, F intersect_op_type, size_t k_idx, size_t k) {
    if (num_induced == 0 || orig.num_edges == 0 || !induced) return;
    nn = orig.nn;
    hash_induced();

    induced_degs = pbbs::new_array_no_init<uintE>(num_induced);
    if (to_save) induced_edges = pbbs::new_array_no_init<uintE*>(num_induced);

    parallel_for (0, num_induced, [&] (size_t j) {
      uintE v = induced[j];
      uintE orig_v_idx = orig.induced_table.find(v, nn);
      uintE v_deg = orig.induced_degs[orig_v_idx];
      uintE* nbhrs = orig.induced_edges[orig_v_idx];
      induced_degs[j] = 0;
      if (to_save) induced_edges[j] = nullptr;
      if (v_deg > 0) {
      if (to_save) {
        induced_edges[j] = pbbs::new_array_no_init<uintE>(std::min((size_t) v_deg, (size_t) num_induced));
        for (size_t l = 0; l < v_deg; l++) {
          if (induced_table.contains(nbhrs[l])) { induced_edges[j][induced_degs[j]] = nbhrs[l]; induced_degs[j]++; }
        }
        if (induced_degs[j] == 0) {
          pbbs::delete_array<uintE>(induced_edges[j], std::min((size_t) v_deg, (size_t) num_induced));
          induced_edges[j] = nullptr;
        }
      } else {
        for (size_t l = 0; l < v_deg; l++) {
          if (induced_table.contains(nbhrs[l])) induced_degs[j]++;
        }
      }
      }
    });
    auto deg_seq = pbbslib::make_sequence(induced_degs, num_induced);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  void alloc_induced(size_t size) {
    if (!induced) induced = pbbs::new_array_no_init<uintE>(size);
  }

  void del() {
    if (!protected_flag && induced) {pbbs::delete_array<uintE>(induced, num_induced);}
    induced = nullptr;
    if (induced_edges && induced_degs) {
      for (size_t i=0; i < num_induced; ++i) {
        if (induced_edges[i]) pbbs::delete_array<uintE>(induced_edges[i], induced_degs[i]);
        induced_edges[i] = nullptr;
      }
      pbbs::delete_array<uintE*>(induced_edges, num_induced);
    }
    induced_edges = nullptr;
    if (induced_degs) {pbbs::delete_array<uintE>(induced_degs, num_induced);}
    induced_degs = nullptr;
    induced_table.del();
  }

  ~FullSpace_csv_hashhash_dyn() {del();}

  static void init() {}
  static void finish() {}
};
*/
// graph: *cdj = offsets, *adj = neighbor list, core = int (core val?)
// subgraph: *n = n[i] is # nodes in G_i, **d = d[i] is degrees of G_i, *adj = adj list, *lab = lab[i] is label of node i, **nodes = nodes[l] is nodes of G_l
// core = max deg
// for each level, (k-1 to 2)
//  the use CSV
// for each node u in induced, and each node v adj to u, if label of v is l, then set  to l-1 and add v to our new level
// then, for each v that we've added to our new subgraph, for each w that was adjacent to it in G, check if w was labeled with l-1; if so, add  to adj list of v

/*

struct FullSpace_csv_hash_dyn {
  size_t num_induced;
  uintE* induced;
  bool protected_flag;
  bool full_flag = true;
  size_t num_edges = 0;
  uintE** induced_edges;
  uintE* induced_degs;
  uintE* induced_labels;
  size_t nn;
  bool orig_flag = false;

  FullSpace_csv_hash_dyn() : num_induced(0), induced(nullptr), protected_flag(false), induced_edges(nullptr), induced_degs(nullptr), induced_labels(nullptr), nn(0) {}

  size_t getDegree(size_t i) { assert(induced_degs); return induced_degs[i]; }
  uintE* getNeighbors(size_t i) { assert(induced_edges); return induced_edges[i]; }

  template <class Graph>
  FullSpace_csv_hash_dyn(Graph& DG, size_t k, size_t i) : protected_flag(true), orig_flag(true) {
    assert (k*DG.n < UINT_E_MAX);
    induced = (uintE*)(DG.get_vertex(i).getOutNeighbors());
    num_induced = DG.get_vertex(i).getOutDegree();
    if (num_induced == 0) return;
    nn = DG.n;

    induced_labels = pbbs::new_array_no_init<uintE>(nn);
    auto induced_labels_f = [=] (size_t j) {new ((void*) (induced_labels+j)) uintE(k*nn);};
    parallel_for(0, nn, induced_labels_f);
    parallel_for(0, num_induced, [&] (size_t j) {induced_labels[induced[j]] = j;});
  
    induced_degs = pbbs::new_array_no_init<uintE>(num_induced);
    induced_edges = pbbs::new_array_no_init<uintE*>(num_induced);

    //for (size_t j=0; j < num_induced; j++) {
    parallel_for (0, num_induced, [&] (size_t j) {
      // intersect neighbors of each vert in induced with induced
      uintE v = induced[j];
      uintE* nbhrs = (uintE*)(DG.get_vertex(v).getOutNeighbors());
      size_t deg = DG.get_vertex(v).getOutDegree();
      induced_degs[j] = 0;
      induced_edges[j] = nullptr;
      if (deg > 0) {
        induced_edges[j] = pbbs::new_array_no_init<uintE>(std::min((size_t) deg, (size_t) num_induced));
        for (size_t l = 0; l < deg; l++) {
          if (induced_labels[nbhrs[l]] < k*nn) {
            induced_edges[j][induced_degs[j]] = nbhrs[l];
            induced_degs[j]++;
          }
        }
      //induced_degs[j] = intersect_op_type(nbhrs, deg, induced, num_induced, true, induced_edges[j]);
      if (induced_degs[j] == 0) {
        pbbs::delete_array<uintE>(induced_edges[j], std::min((size_t) deg, (size_t) num_induced));
        induced_edges[j] = nullptr;
      }
      }
    });
    auto deg_seq = pbbslib::make_sequence(induced_degs, num_induced);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  template<class O, class F>
  void prune(O& orig, bool to_save, F intersect_op_type, size_t k_idx, size_t k) {
    if (num_induced == 0 || orig.num_edges == 0 || !induced) return;
    nn = orig.nn;
    induced_labels = orig.induced_labels;
    parallel_for(0, num_induced, [&] (size_t j) {induced_labels[induced[j]] = k_idx*nn + j;});

    induced_degs = pbbs::new_array_no_init<uintE>(num_induced);
    if (to_save) induced_edges = pbbs::new_array_no_init<uintE*>(num_induced);
    //for (size_t j=0; j < num_induced; j++) {
    parallel_for (0, num_induced, [&] (size_t j) {
      uintE v = induced[j];
      uintE v_deg = orig.induced_degs[orig.induced_labels[v] - (k_idx-1)*nn];
      uintE* nbhrs = orig.induced_edges[orig.induced_labels[v] - (k_idx-1)*nn];
      induced_degs[j] = 0;
      if (to_save) induced_edges[j] = nullptr;
      if (v_deg > 0) {
      if (to_save) {
        induced_edges[j] = pbbs::new_array_no_init<uintE>(std::min((size_t) v_deg, (size_t) num_induced));
        for (size_t l = 0; l < v_deg; l++) {
          if (induced_labels[nbhrs[l]] < (k_idx+1)*nn && induced_labels[nbhrs[l]] >= k_idx*nn) { induced_edges[j][induced_degs[j]] = nbhrs[l]; induced_degs[j]++; }
        }
        if (induced_degs[j] == 0) {
          pbbs::delete_array<uintE>(induced_edges[j], std::min((size_t) v_deg, (size_t) num_induced));
          induced_edges[j] = nullptr;
        }
      } else {
        for (size_t l = 0; l < v_deg; l++) {
          if (induced_labels[nbhrs[l]] < (k_idx+1)*nn && induced_labels[nbhrs[l]] >= k_idx*nn)
            induced_degs[j]++;
        }
      }
      }
    });
    auto deg_seq = pbbslib::make_sequence(induced_degs, num_induced);
    num_edges = pbbslib::reduce_add(deg_seq);
  }

  void alloc_induced(size_t size) {
    if (!induced) induced = pbbs::new_array_no_init<uintE>(size);
  }

  void del() {
    // TODO must reset labels somehow
    if (induced) {
      parallel_for(0, num_induced, [&] (size_t j) {induced_labels[induced[j]] = k_idx*nn + j;});
    }
    if (!protected_flag && induced) {pbbs::delete_array<uintE>(induced, num_induced);}
    induced = nullptr;
    if (induced_edges && induced_degs) {
      for (size_t i=0; i < num_induced; ++i) {
        if (induced_edges[i]) pbbs::delete_array<uintE>(induced_edges[i], induced_degs[i]);
        induced_edges[i] = nullptr;
      }
      pbbs::delete_array<uintE*>(induced_edges, num_induced);
    }
    induced_edges = nullptr;
    if (induced_degs) {pbbs::delete_array<uintE>(induced_degs, num_induced);}
    induced_degs = nullptr;
    if (orig_flag && induced_labels) {pbbs::delete_array<uintE>(induced_labels, nn);}
    induced_labels = nullptr;
  }

  ~FullSpace_csv_hash_dyn() {del();}

  static void init() {}
  static void finish() {}
};
*/
/*
struct FullSpace_csv_dyn {
  bool full_flag = true;
  size_t num_edges = 0;
  bool protected_flag = false;
  size_t num_induced = 0;
  uintE* induced = nullptr;
  uintE* induced_offsets = nullptr;
  uintE* induced_edges = nullptr;
  static constexpr auto intersect_op_type = lstintersect_par_struct{};
  bool pruned = false;
  size_t nn = 0;

  FullSpace_csv_dyn() {}

  template <class Graph>
  FullSpace_csv_dyn(Graph& DG, size_t min_deg, size_t i) {
    protected_flag = true;
    induced = (uintE*)(DG.get_vertex(i).getOutNeighbors());
    num_induced = DG.get_vertex(i).getOutDegree();
    prune_init(DG, min_deg);
    nn = DG.n;
  }

  void alloc_induced(size_t size) {
    if (!induced) induced = pbbs::new_array_no_init<uintE>(size);
  }

  void del() {
    if (!protected_flag && induced) {pbbs::delete_array<uintE>(induced, num_induced); induced = nullptr;}
    if (induced_edges) {pbbs::delete_array<uintE>(induced_edges, num_edges); induced_edges = nullptr;}
    if (induced_offsets) {pbbs::delete_array<uintE>(induced_offsets, num_induced + 1); induced_offsets = nullptr;}
  }

  ~FullSpace_csv_dyn() {del();}

  uintE* getNeighbors(size_t i) {
    assert (i < num_induced);
    assert (induced_offsets);
    assert (induced_edges);
    return induced_edges + induced_offsets[i];
  }

  size_t getDegree(size_t i) {
    return induced_offsets[i+1] - induced_offsets[i];
  }

  template <class O>
  void prune(O& orig, size_t min_deg) {
    if (num_induced == 0) return;
  
    nn = orig.nn;
    if (nn > 0) {
    for(size_t i=0; i < num_induced; i++) {
        uintE orig_v_offset = orig.induced_offsets[i];
        uintE orig_v_deg = orig.induced_offsets[i + 1] - orig_v_offset;
        uintE* orig_v_edges = orig.induced_edges + orig_v_offset;
        for(size_t j=0; j < orig_v_deg; j++) {
          assert (orig_v_edges[j] < nn);
        }
      }
    }

    if (orig.num_edges == 0) {del(); num_induced = 0; return;}
    prune_seq_orig(orig, orig.num_edges, min_deg);
    //prune_par(get_edges, min_deg);
    if (num_induced > 0) num_edges = induced_offsets[num_induced];
    else del();
  }

  template <class Graph>
  void prune_init(Graph& DG, size_t min_deg = 0) {
    if (num_induced == 0) return;

    auto deg_f = [&](size_t i) { return DG.get_vertex(induced[i]).getOutDegree(); };
    size_t orig_num_edges = reduce(pbbslib::make_sequence<uintE>(num_induced, deg_f), pbbslib::addm<uintE>());

    if (orig_num_edges == 0) {del(); num_induced = 0; return;}
    prune_seq(DG, orig_num_edges, min_deg);

    if (num_induced > 0) {
      num_edges = induced_offsets[num_induced];
  
      for(size_t i=0; i < num_induced; i++) {
        uintE orig_v_offset = induced_offsets[i];
        uintE orig_v_deg = induced_offsets[i + 1] - orig_v_offset;
        uintE* orig_v_edges = induced_edges + orig_v_offset;
        for(size_t j=0; j < orig_v_deg; j++) {
          assert (orig_v_edges[j] < DG.n);
        }
      }
    } else del();
  }

  template <class O>
  void prune_seq_orig(O& orig, size_t orig_num_edges, size_t min_deg = 0) {
    if (!induced || num_induced == 0) return;
    if (orig_num_edges == 0) {num_induced = 0; return;}

    // first, construct induced_offsets
    induced_offsets = pbbs::new_array_no_init<uintE>(num_induced + 1);
    induced_edges = pbbs::new_array_no_init<uintE>(orig_num_edges);
    induced_offsets[0] = 0;
    uintE offsets = 0;

    for (size_t i = 0; i < num_induced; ++i) {
      uintE v = induced[i];
      size_t orig_v_deg;
      uintE* orig_v_edges;
        auto lte = [] (const uintE& l, const uintE& r) { return l < r; };
      auto orig_induced_seq = pbbslib::make_sequence<uintE>(orig.induced, orig.num_induced);
      size_t find_idx = pbbs::binary_search(orig_induced_seq, v, lte);
      assert (find_idx < orig.num_induced);
      uintE orig_v_offset = orig.induced_offsets[find_idx];
      size_t out_size = intersect_op_type(orig.induced_edges + orig_v_offset, orig.induced_offsets[find_idx + 1] - orig_v_offset, induced, num_induced, true, induced_edges + offsets);
      offsets = out_size >= min_deg ? out_size + offsets : offsets;
      induced_offsets[i + 1] = offsets;
    }

    if (induced_offsets[num_induced] == 0) {num_induced = 0; return;}
  }

  template <class Graph>
  void prune_seq(Graph& DG, size_t orig_num_edges, size_t min_deg = 0) {
    if (!induced || num_induced == 0) return;
    if (orig_num_edges == 0) {num_induced = 0; return;}

    // first, construct induced_offsets
    induced_offsets = pbbs::new_array_no_init<uintE>(num_induced + 1);
    induced_edges = pbbs::new_array_no_init<uintE>(orig_num_edges+1);
    induced_offsets[0] = 0;
    uintE offsets = 0;

    for (size_t i = 0; i < num_induced; ++i) {
      uintE v = induced[i];
      size_t out_size = intersect_op_type((uintE*)(DG.get_vertex(v).getOutNeighbors()), DG.get_vertex(v).getOutDegree(), induced, num_induced, true, induced_edges + offsets);
      offsets = (out_size >= min_deg) ? out_size + offsets : offsets;
      induced_offsets[i + 1] = offsets;
    }

    if (induced_offsets[num_induced] == 0) {num_induced = 0; return;}
  }

  template <class E>
  void prune_par(E& get_edges, size_t min_deg = 0) {
    if (!induced || num_induced == 0) return;

    // first, construct induced_offsets
    // TODO: alternatively, if we do everything sequentially, we can skip this and just take the hit
    induced_offsets = pbbs::new_array_no_init<uintE>(num_induced + 1);
    induced_offsets[num_induced] = 0;

    parallel_for (0, num_induced, [&] (size_t i) { 
    //for (size_t i = 0; i < num_induced; ++i) {
      uintE v = induced[i];
      auto edges_tup = get_edges(v);
      auto orig_v_deg = std::get<1>(edges_tup);
      if (orig_v_deg > 0) {
        uintE* orig_v_edges = std::get<0>(edges_tup);

        uintE* out_ptr = nullptr;
        if (intersect_op_type.count_space_flag) out_ptr = pbbs::new_array_no_init<uintE>(std::min((size_t) orig_v_deg, (size_t) num_induced));
        size_t out_size = intersect_op_type(orig_v_edges, orig_v_deg, induced, num_induced, false, out_ptr);
        if (intersect_op_type.count_space_flag) pbbs::delete_array<uintE>(out_ptr, std::min((size_t) orig_v_deg, (size_t) num_induced));
        induced_offsets[i] = out_size >= min_deg ? out_size : 0;
      }
      else induced_offsets[i] = 0;
    });
    scan_inplace((pbbslib::make_sequence<uintE>(induced_offsets, num_induced + 1)).slice(), pbbslib::addm<uintE>());
  
    if (induced_offsets[num_induced] == 0) {del(); num_induced = 0; return;}
  
    induced_edges = pbbs::new_array_no_init<uintE>(induced_offsets[num_induced]);
    parallel_for (0, num_induced, [&] (size_t i) { 
    //for (size_t i=0; i < num_induced; ++i) {
      if (induced_offsets[i+1] - induced_offsets[i] > 0) {
        uintE v = induced[i];
        auto edges_tup = get_edges(v);
        auto orig_v_deg = std::get<1>(edges_tup);
        auto orig_v_edges = std::get<0>(edges_tup);

        uintE v_offset = induced_offsets[i];
        uintE* v_edges = induced_edges + v_offset;
        intersect_op_type(orig_v_edges, orig_v_deg, induced, num_induced, true, v_edges);
      }
    });
    // for each element in induced
    // find out it's edges in orig_space
    // delete all edges in orig_space that do not appear in induced
    // if < min_deg, then mark element to be deleted (just set deg to 0, keep no edges)
    // otherwise, we'll record its edges
  }

  

  static void init() {}
  static void finish() {}
};*/




#endif