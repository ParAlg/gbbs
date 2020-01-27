#pragma once

/* TODO: describe what this file does */

#include "intersect.h"
#include "pbbslib/seq.h"

namespace induced_hybrid {

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
  inline size_t KCliqueDir_fast_hybrid_rec(Graph& DG, size_t k_idx, size_t k, HybridSpace_lw* induced, F base_f) {
    //if (k == 2) return induced->num_edges;
    size_t num_induced = induced->num_induced[k_idx-1];
    if (num_induced == 0) return 0;
    uintE* prev_induced = induced->induced + induced->nn * (k_idx - 1);

    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx; }

    if (k_idx + 1 == k) {
      size_t counts = 0;
      for (size_t i=0; i < num_induced; i++) {
        uintE vtx = prev_induced[i];
        if (induced->use_base) induced->base[k_idx] = induced->relabel[vtx];
        //  get neighbors of vtx
        uintE* intersect = induced->induced_edges + vtx * induced->nn;
        for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
          if (induced->labels[intersect[j]] == k_idx) {
            counts++;
            if (induced->use_base) {
              induced->base[k] = induced->relabel[intersect[j]];
              base_f(induced->base);
            }
          } 
        }
      }
      for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
      return counts;
    }

    size_t total_ct = 0;
    for (size_t i=0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      if (induced->use_base) induced->base[k_idx] = induced->relabel[vtx]; // TODO problem w/storing base -- we've relabeled our vert w/relabeling: check base is correct
      uintE* intersect = induced->induced_edges + vtx * induced->nn;
      uintE* out = induced->induced + induced->num_induced[0] * k_idx;
      uintE count = 0;
      for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
        if (induced->labels[intersect[j]] == k_idx) {
          out[count] = intersect[j];
          count++;
        }
      }
      induced->num_induced[k_idx] = count;
      if (induced->num_induced[k_idx] > k - k_idx - 1) total_ct += KCliqueDir_fast_hybrid_rec(DG, k_idx + 1, k, induced, base_f);
    }

    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
    return total_ct;
  }

  /*template <class Graph>
  inline size_t CountCliques_unbalanced(Graph& DG, size_t k) {
    sequence<size_t> tots = sequence<size_t>::no_init(DG.n);
    size_t max_deg = get_max_deg(DG);
    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, DG.n); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del(); 
    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
      if (DG.get_vertex(i).getOutDegree() != 0) {
        induced->setup(DG, k, i);
        tots[i] = KCliqueDir_fast_hybrid_rec(DG, 1, k, induced);
      } else tots[i] = 0;
    } );

    return pbbslib::reduce_add(tots);
  }*/


  template <class Graph, class F>
  inline size_t CountCliques(Graph& DG, size_t k, F base_f, bool use_base=false, bool label=true) {
    timer t; t.start();
    using W = typename Graph::weight_type;
    auto parallel_work = sequence<size_t>(DG.n);
    {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      return DG.get_vertex(v).getOutDegree();
    };
    par_for(0, DG.n, [&] (size_t i) {
      auto monoid = pbbslib::addm<size_t>();
      parallel_work[i] = DG.get_vertex(i).template reduceOutNgh<size_t>(i, map_f, monoid);
    });
    }
    size_t total_work = pbbslib::scan_add_inplace(parallel_work.slice());

    size_t block_size = 50000;
    size_t n_blocks = total_work/block_size + 1;
    size_t work_per_block = total_work / n_blocks;
    n_blocks = (total_work/work_per_block) + 1;
    double tt = t.stop();
    std::cout << "##### Triangle scheduling: " << tt << std::endl;

    timer t2; t2.start();
    sequence<size_t> tots = sequence<size_t>::no_init(n_blocks); //DG.n
    size_t max_deg = get_max_deg(DG);
    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, DG.n, label, use_base); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del(); 
    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, n_blocks, [&](size_t j, HybridSpace_lw* induced) {
      size_t start = j * work_per_block;
      size_t end = (j + 1) * work_per_block;
      auto less_fn = std::less<size_t>();
      size_t start_ind = pbbslib::binary_search(parallel_work, start, less_fn);
      size_t end_ind = pbbslib::binary_search(parallel_work, end, less_fn);
      tots[j] = 0;
      for (size_t i=start_ind; i < end_ind; i++) {
        if (DG.get_vertex(i).getOutDegree() != 0) {
          induced->setup(DG, k, i);
          tots[j] += KCliqueDir_fast_hybrid_rec(DG, 1, k, induced, base_f);
        }
      }
    }, 1, false);
    double tt2 = t2.stop();
    std::cout << "##### Actual counting: " << tt2 << std::endl;


    /*sequence<size_t> tots = sequence<size_t>::no_init(DG.n);
    size_t max_deg = get_max_deg(DG);
    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, DG.n); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del(); 
    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
      if (DG.get_vertex(i).getOutDegree() != 0) {
        induced->setup(DG, k, i);
        tots[i] = KCliqueDir_fast_hybrid_rec(DG, 1, k, induced);
      } else tots[i] = 0;
    }, 1, false);*/
    /*
    #pragma omp parallel private(induced) reduction(+:n)
    {
    induced = new HybridSpace_lw(max_deg, k);
    #pragma omp for schedule(dynamic, 1) nowait
    for (size_t i=0; i < DG.n; ++i) {
      if (DG.get_vertex(i).getOutDegree() != 0) {
        induced->setup(DG, k, i);
        n += KCliqueDir_fast_hybrid_rec(DG, 1, k, induced);
      }
    }

    if (induced != nullptr) { induced->del(); delete induced; }

    }*/

    return pbbslib::reduce_add(tots);
  }

} // namespace induced_neighborhood
