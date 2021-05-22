#pragma once

#include <math.h>
#include <limits>

// Library dependencies
#include "gbbs/bucket.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"
#include "gbbs/pbbslib/dyn_arr.h"
#include "gbbs/pbbslib/sparse_table.h"
#include "gbbs/pbbslib/sparse_additive_map.h"
#include "pbbslib/assert.h"
#include "pbbslib/list_allocator.h"
#include "pbbslib/integer_sort.h"

// Clique files
#include "benchmarks/CliqueCounting/intersect.h"
#include "benchmarks/CliqueCounting/induced_intersection.h"
#include "benchmarks/CliqueCounting/induced_neighborhood.h"
#include "benchmarks/CliqueCounting/induced_hybrid.h"
#include "benchmarks/CliqueCounting/induced_split.h"
#include "benchmarks/CliqueCounting/relabel.h"
  
namespace gbbs {

  extern int nd_global_shift_factor;
  
  struct hash128 {
    inline size_t operator () (unsigned __int128 t) const {
      size_t l = t >> 64;
      unsigned __int128 mask = l << 64;
      size_t r = t & (~mask);
      return pbbs::hash_combine(pbbslib::hash64_2(l), pbbslib::hash64_2(r));
    }
  };

  struct nhash64 {
    inline size_t operator () (unsigned long long t) const {
      return pbbslib::hash64_2(t);
    }
  };

  struct nhash32 {
    inline size_t operator () (unsigned int t) const {
      return pbbslib::hash32(t);
    }
  };

  inline uintE middleOfThree(uintE a, uintE b, uintE c){
    // Checking for b
    if ((a < b && b < c) || (c < b && b < a)) return b;
    // Checking for a
    else if ((b < a && a < c) || (c < a && a < b)) return a;
    else return c;
  }
  
  template <class F>
  inline uintE middleOfThree(uintE a, uintE b, uintE c, F& func){
    // Checking for b
    if (( func(a, b) && func(b, c)) || (func(c, b) && func(b, a))) return b;
    // Checking for a
    else if ((func(b, a) && func(a, c)) || (func(c, a) && func(a, b))) return a;
    else return c;
  }

  template <class F>
  inline uintE minOfThree(uintE a, uintE b, uintE c, F& func){
    if ( func(a,b) && func(a, c)) return a;
    if (func(b,a) && func(b,c)) return b;
    return c;
  }

  template <class F>
  inline uintE maxOfThree(uintE a, uintE b, uintE c, F& func){
    if ( func(b,a) && func(c,a)) return a;
    if (func(a,b) && func(c,b)) return b;
    return c;
  }
 
 


  template<class Graph>
  size_t get_max_deg3(Graph& DG) {
    size_t max_deg = 0;
    parallel_for(0, DG.n, [&] (size_t i) {
      size_t deg = DG.get_vertex(i).getOutDegree();
      pbbs::write_min(&max_deg, deg, std::greater<size_t>());
    });
    return max_deg;
  }

template <class Graph>
bool is_edge3(Graph& DG, uintE v, uintE u) {
  using W = typename Graph::weight_type;
  bool is = false;
  auto map_f = [&] (const uintE& src, const uintE& vv, const W& wgh) {
    if (vv == u) is = true;
  };
  DG.get_vertex(v).mapOutNgh(v, map_f, false);
  auto map_f2 = [&] (const uintE& src, const uintE& vv, const W& wgh) {
    if (vv == v) is = true;
  };
  DG.get_vertex(u).mapOutNgh(u, map_f2, false);
  return is;
}

  template <class Graph, class F, class B>
  inline size_t NKCliqueDir_fast_hybrid_rec(Graph& DG, size_t k_idx, size_t k, HybridSpace_lw* induced, F& base_f,
  B& base) {
    size_t num_induced = induced->num_induced[k_idx-1];
    if (num_induced == 0) return 0;
    uintE* prev_induced = induced->induced + induced->nn * (k_idx - 1);

    if (k_idx == k) {
      size_t tmp_counts = 0;
      //assert(num_induced == induced->nn);
      //std::cout << "num_indced: "<< num_induced << std::endl; fflush(stdout);
      for (size_t j=0; j < num_induced; j++) {
        uintE vtx = prev_induced[j];
        //assert(vtx == j);
        base[k] = induced->relabel[vtx];
        
        if (base[k] != UINT_E_MAX) {
          /*for (int i = 0; i < k + 1; i++) {
            int i1 = (i + 1) % (k + 1);
            if (!is_edge3(DG, base[i], base[i1])) {
              std::cout << "Flip: " << is_edge3(DG, base[i1], base[i]) << std::endl;
              std::cout << "i: " << i << ", i1: " << i1 << ", base i: "<< base[i] << ", base i1: " << base[i1] << std::endl;
              fflush(stdout);
            }
            assert(is_edge3(DG, base[i], base[i1]));
          }*/
          //std::cout << "Exists base" << std::endl; fflush(stdout);
          base_f(base);
          tmp_counts++;
        }
      }
      return tmp_counts;
    }

    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx; }

    if (k_idx + 1 == k) {
      size_t counts = 0;
      for (size_t i=0; i < num_induced; i++) {
        uintE vtx = prev_induced[i];
        //  get neighbors of vtx
        uintE* intersect = induced->induced_edges + vtx * induced->nn;
        size_t tmp_counts = 0;
        base[k_idx] = induced->relabel[vtx];
        if (base[k_idx] != UINT_E_MAX) {
        for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
          if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
            base[k] = induced->relabel[intersect[j]];
            if (base[k] != UINT_E_MAX) {
              tmp_counts++;
              base_f(base);
            }
          }
        }
        counts += tmp_counts;
        }
      }
      for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
      return counts;
    }

size_t total_ct = 0;
    for (size_t i=0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      uintE* intersect = induced->induced_edges + vtx * induced->nn;
      uintE* out = induced->induced + induced->nn * k_idx;
      uintE count = 0;
      for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
        if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
          out[count] = intersect[j];
          count++;
        }
      }
      induced->num_induced[k_idx] = count;
      if (induced->num_induced[k_idx] > k - k_idx - 1) {
        base[k_idx] = induced->relabel[vtx];
        if (base[k_idx] != UINT_E_MAX) {
          auto curr_counts = NKCliqueDir_fast_hybrid_rec(DG, k_idx + 1, k, induced, base_f, base);
          total_ct += curr_counts;
        }
      }
    }
    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
    return total_ct;
  }

} // end namespace gbbs