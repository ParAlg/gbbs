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


  template<class Graph>
  size_t get_max_deg3(Graph& DG) {
    size_t max_deg = 0;
    parallel_for(0, DG.n, [&] (size_t i) {
      size_t deg = DG.get_vertex(i).getOutDegree();
      pbbs::write_min(&max_deg, deg, std::greater<size_t>());
    });
    return max_deg;
  }

  template <class Graph, class F>
  inline size_t NKCliqueDir_fast_hybrid_rec(Graph& DG, size_t k_idx, size_t k, HybridSpace_lw* induced, F base_f,
  sequence<uintE>& base) {
    size_t num_induced = induced->num_induced[k_idx-1];
    if (num_induced == 0) return 0;
    uintE* prev_induced = induced->induced + induced->nn * (k_idx - 1);

    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx; }

    if (k_idx + 1 == k) {
      size_t counts = 0;
      for (size_t i=0; i < num_induced; i++) {
        uintE vtx = prev_induced[i];
        //  get neighbors of vtx
        uintE* intersect = induced->induced_edges + vtx * induced->nn;
        size_t tmp_counts = 0;
        base[k_idx] = induced->relabel[vtx];
        for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
          if (static_cast<size_t>(induced->labels[intersect[j]]) == k_idx) {
            base[k] = induced->relabel[intersect[j]];
            tmp_counts++;
            base_f(base);
          }
        }
        counts += tmp_counts;
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
        auto curr_counts = NKCliqueDir_fast_hybrid_rec(DG, k_idx + 1, k, induced, base_f, base);
        total_ct += curr_counts;
      }
    }
    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
    return total_ct;
  }

} // end namespace gbbs