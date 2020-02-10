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
  inline size_t KCliqueDir_fast_hybrid_rec(Graph& DG, size_t k_idx, size_t k, HybridSpace_lw* induced, F base_f,
    size_t recursive_level=0) {
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
        for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
          if (induced->labels[intersect[j]] == k_idx) {
            tmp_counts++;
            if (induced->use_base) base_f(induced->relabel[intersect[j]], 1);
          } 
        }
        if (induced->use_base && tmp_counts > 0) base_f(induced->relabel[vtx], tmp_counts);
        counts += tmp_counts;
      }
      for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
      return counts;
    }
size_t total_ct = 0;
if (recursive_level < k_idx || num_induced < 2) {
    for (size_t i=0; i < num_induced; ++i) {
      uintE vtx = prev_induced[i];
      uintE* intersect = induced->induced_edges + vtx * induced->nn;
      uintE* out = induced->induced + induced->nn * k_idx;
      uintE count = 0;
      for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
        if (induced->labels[intersect[j]] == k_idx) {
          out[count] = intersect[j];
          count++;
        }
      }
      induced->num_induced[k_idx] = count;
      if (induced->num_induced[k_idx] > k - k_idx - 1) {
        auto curr_counts = KCliqueDir_fast_hybrid_rec(DG, k_idx + 1, k, induced, base_f, recursive_level);
        total_ct += curr_counts;
        if (induced->use_base && curr_counts > 0) base_f(induced->relabel[vtx], curr_counts);
      }
    }
} else {
    sequence<size_t> tots = sequence<size_t>::no_init(num_induced);
    auto init_induced = [&](HybridSpace_lw* induced2) { induced2->alloc(induced->nn, k-k_idx, DG.n, false, induced->use_base, false); };
    auto finish_induced = [&](HybridSpace_lw* induced2) { if (induced2 != nullptr) { delete induced2; } };
    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, num_induced, [&](size_t i, HybridSpace_lw* induced2) {
      induced2->copy(induced);
      uintE vtx = prev_induced[i];
      uintE* intersect = induced->induced_edges + vtx * induced->nn;
      uintE* out = induced2->induced;
      uintE count = 0;
      for (size_t j=0; j < induced->induced_degs[vtx]; j++) {
        if (induced->labels[intersect[j]] == k_idx) {
          out[count] = intersect[j];
          count++;
        }
      }
      induced2->num_induced[0] = count;
      if (count > k - k_idx - 1) {
        tots[i] = KCliqueDir_fast_hybrid_rec(DG, 1, k-k_idx, induced2, base_f);
        if (induced->use_base && tots[i] > 0) base_f(induced->relabel[vtx], tots[i]);
      } else tots[i] = 0;
    });
    total_ct += pbbslib::reduce_add(tots);
}
    //for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
    for (size_t i=0; i < num_induced; i++) { induced->labels[prev_induced[i]] = k_idx - 1; }
    return total_ct;
  }


  template <class Graph, class F>
  inline size_t CountCliques(Graph& DG, size_t k, F base_f, bool use_base=false, bool label=true, long recursive_level=0) {
    timer t2; t2.start();
    using W = typename Graph::weight_type;

if (recursive_level != 1) {
    sequence<size_t> tots = sequence<size_t>::no_init(DG.n);
    size_t max_deg = get_max_deg(DG);
    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, DG.n, label, use_base); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del(); 
    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
      if (DG.get_vertex(i).getOutDegree() != 0) {
        induced->setup(DG, k, i);
        tots[i] = KCliqueDir_fast_hybrid_rec(DG, 1, k, induced, base_f, recursive_level);
        if (induced->use_base && tots[i] > 0) base_f(i, tots[i]);
      } else tots[i] = 0;
    }, 1, false);
    double tt2 = t2.stop();
    std::cout << "##### Actual counting: " << tt2 << std::endl;

    return pbbslib::reduce_add(tots);
}

   
   size_t max_deg = get_max_deg(DG);
   sequence<size_t> degs = sequence<size_t>::no_init(DG.n+1);
    parallel_for(0, DG.n, [&] (size_t i) { degs[i] = DG.get_vertex(i).getOutDegree();});
    degs[DG.n] = 0;
    size_t num_edges = pbbslib::scan_add_inplace(degs.slice());
    num_edges = degs[DG.n];
    sequence<size_t> tots = sequence<size_t>::no_init(num_edges);

    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, DG.n, label, use_base); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };
  parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, num_edges, [&](size_t j, HybridSpace_lw* induced) {
    // to find i and ngh, binary search for j in degs; the index - 1 is i
    // then, DG.get_vertex(i).getOutNeighbor(j - degs[index-1]) is ngh
    auto less_fn = [&](size_t a, size_t b){ return a <= b; };
    size_t idx = pbbslib::binary_search(degs, j, less_fn);
    auto i = idx - 1;
    auto ngh = DG.get_vertex(i).getOutNeighbor(j - degs[idx - 1]);
    induced->setup_edge(DG,k,i,ngh);
    tots[j] = KCliqueDir_fast_hybrid_rec(DG, 1, k-1, induced, base_f, 0);
    if (use_base && tots[j] > 0) { base_f(ngh, tots[j]); base_f(i, tots[j]); } 
   }, 1, false);
    double tt2 = t2.stop();
    std::cout << "##### Actual counting: " << tt2 << std::endl;

    return pbbslib::reduce_add(tots);
  }




  template <class Graph, class F>
  inline size_t CountCliques_gen(Graph& DG, size_t k, F base_f, bool use_base=false, bool label=true, long recursive_level=0) {
    timer t2; t2.start();
    using W = typename Graph::weight_type;

    sequence<size_t> tots = sequence<size_t>::no_init(DG.n);
    size_t max_deg = get_max_deg(DG);
    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, DG.n, label, use_base); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; 

 switch (k) {
 case 3:  {
 parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
 if (DG.get_vertex(i).getOutDegree() != 0) {
 induced->setup(DG, k, i);
 size_t k_idx = 1;
 size_t counts; uintE count; uintE* out; uintE* prev_induced; size_t num_induced; uintE vtx; uintE* intersect=nullptr;
 k_idx = 1;
 num_induced = induced->num_induced[k_idx-1];
 size_t acurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t b=0; b < num_induced; b++) { induced->labels[prev_induced[b]] = k_idx; }
 for (size_t b=0; b < num_induced; ++b) {
 auto bvtx = prev_induced[b];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t c=0; c < induced->induced_degs[vtx]; c++) {
 if (induced->labels[intersect[c]] == k_idx) {
 out[count] = intersect[c];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 2;
 num_induced = induced->num_induced[k_idx-1];
 size_t bcurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t c=0; c < num_induced; c++) { induced->labels[prev_induced[c]] = k_idx; }
 for (size_t c=0; c < num_induced; ++c) {
 auto cvtx = prev_induced[c];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t d=0; d < induced->induced_degs[vtx]; d++) {
 if (induced->labels[intersect[d]] == k_idx) {
 out[count] = intersect[d];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 3;
 num_induced = induced->num_induced[k_idx-1];
 size_t ccurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t d=0; d < num_induced; d++) { induced->labels[prev_induced[d]] = k_idx; }
 for (size_t d=0; d < num_induced; ++d) {
 auto dvtx = prev_induced[d];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t e=0; e < induced->induced_degs[vtx]; e++) {
 if (induced->labels[intersect[e]] == k_idx) {
 out[count] = intersect[e];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 4;
 num_induced = induced->num_induced[k_idx-1];
 size_t dcurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t e=0; e < num_induced; e++) { induced->labels[prev_induced[e]] = k_idx; }
 for (size_t e=0; e < num_induced; ++e) {
 auto evtx = prev_induced[e];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t f=0; f < induced->induced_degs[vtx]; f++) {
 if (induced->labels[intersect[f]] == k_idx) {
 out[count] = intersect[f];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 5;
 num_induced = induced->num_induced[k_idx-1];
 if (num_induced == 0) {continue;}
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 size_t ecurr_counts=0;
 counts = 0;
 for (size_t x=0; x < num_induced; i++) {
 vtx = prev_induced[x];
 intersect = induced->induced_edges + vtx * induced->nn;
 size_t tmp_counts = 0;
 for (size_t y=0; y < induced->induced_degs[vtx]; y++) {
 if (induced->labels[intersect[y]] == k_idx) {
 tmp_counts++;
 } 
 }
 etots[i] += tmp_counts;
 }
 for (size_t z=0; z < num_induced; z++) { induced->labels[prev_induced[z]] = k_idx - 1; }
 } }
 k_idx = 4;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t f=0; f<num_induced; f++) { induced->labels[prev_induced[f]] = k_idx - 1; }
  }
 } }
 k_idx = 3;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t e=0; e<num_induced; e++) { induced->labels[prev_induced[e]] = k_idx - 1; }
  }
 } }
 k_idx = 2;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t d=0; d<num_induced; d++) { induced->labels[prev_induced[d]] = k_idx - 1; }
  }
 } }
 k_idx = 1;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t c=0; c<num_induced; c++) { induced->labels[prev_induced[c]] = k_idx - 1; }
  }if (induced->use_base && tots[i] > 0) base_f(i, tots[i]);} else tots[i] = 0;}, 1, false);
 break; }
 case 4:  {
 parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
 if (DG.get_vertex(i).getOutDegree() != 0) {
 induced->setup(DG, k, i);
 size_t k_idx = 1;
 size_t counts; uintE count; uintE* out; uintE* prev_induced; size_t num_induced; uintE vtx; uintE* intersect=nullptr;
 k_idx = 1;
 num_induced = induced->num_induced[k_idx-1];
 size_t acurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t b=0; b < num_induced; b++) { induced->labels[prev_induced[b]] = k_idx; }
 for (size_t b=0; b < num_induced; ++b) {
 auto bvtx = prev_induced[b];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t c=0; c < induced->induced_degs[vtx]; c++) {
 if (induced->labels[intersect[c]] == k_idx) {
 out[count] = intersect[c];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 2;
 num_induced = induced->num_induced[k_idx-1];
 size_t bcurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t c=0; c < num_induced; c++) { induced->labels[prev_induced[c]] = k_idx; }
 for (size_t c=0; c < num_induced; ++c) {
 auto cvtx = prev_induced[c];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t d=0; d < induced->induced_degs[vtx]; d++) {
 if (induced->labels[intersect[d]] == k_idx) {
 out[count] = intersect[d];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 3;
 num_induced = induced->num_induced[k_idx-1];
 size_t ccurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t d=0; d < num_induced; d++) { induced->labels[prev_induced[d]] = k_idx; }
 for (size_t d=0; d < num_induced; ++d) {
 auto dvtx = prev_induced[d];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t e=0; e < induced->induced_degs[vtx]; e++) {
 if (induced->labels[intersect[e]] == k_idx) {
 out[count] = intersect[e];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 4;
 num_induced = induced->num_induced[k_idx-1];
 size_t dcurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t e=0; e < num_induced; e++) { induced->labels[prev_induced[e]] = k_idx; }
 for (size_t e=0; e < num_induced; ++e) {
 auto evtx = prev_induced[e];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t f=0; f < induced->induced_degs[vtx]; f++) {
 if (induced->labels[intersect[f]] == k_idx) {
 out[count] = intersect[f];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 5;
 num_induced = induced->num_induced[k_idx-1];
 if (num_induced == 0) {continue;}
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 size_t ecurr_counts=0;
 counts = 0;
 for (size_t x=0; x < num_induced; i++) {
 vtx = prev_induced[x];
 intersect = induced->induced_edges + vtx * induced->nn;
 size_t tmp_counts = 0;
 for (size_t y=0; y < induced->induced_degs[vtx]; y++) {
 if (induced->labels[intersect[y]] == k_idx) {
 tmp_counts++;
 } 
 }
 etots[i] += tmp_counts;
 }
 for (size_t z=0; z < num_induced; z++) { induced->labels[prev_induced[z]] = k_idx - 1; }
 } }
 k_idx = 4;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t f=0; f<num_induced; f++) { induced->labels[prev_induced[f]] = k_idx - 1; }
  }
 } }
 k_idx = 3;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t e=0; e<num_induced; e++) { induced->labels[prev_induced[e]] = k_idx - 1; }
  }
 } }
 k_idx = 2;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t d=0; d<num_induced; d++) { induced->labels[prev_induced[d]] = k_idx - 1; }
  }
 } }
 k_idx = 1;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t c=0; c<num_induced; c++) { induced->labels[prev_induced[c]] = k_idx - 1; }
  }if (induced->use_base && tots[i] > 0) base_f(i, tots[i]);} else tots[i] = 0;}, 1, false);
 break; }
 case 5:  {
 parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
 if (DG.get_vertex(i).getOutDegree() != 0) {
 induced->setup(DG, k, i);
 size_t k_idx = 1;
 size_t counts; uintE count; uintE* out; uintE* prev_induced; size_t num_induced; uintE vtx; uintE* intersect=nullptr;
 k_idx = 1;
 num_induced = induced->num_induced[k_idx-1];
 size_t acurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t b=0; b < num_induced; b++) { induced->labels[prev_induced[b]] = k_idx; }
 for (size_t b=0; b < num_induced; ++b) {
 auto bvtx = prev_induced[b];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t c=0; c < induced->induced_degs[vtx]; c++) {
 if (induced->labels[intersect[c]] == k_idx) {
 out[count] = intersect[c];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 2;
 num_induced = induced->num_induced[k_idx-1];
 size_t bcurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t c=0; c < num_induced; c++) { induced->labels[prev_induced[c]] = k_idx; }
 for (size_t c=0; c < num_induced; ++c) {
 auto cvtx = prev_induced[c];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t d=0; d < induced->induced_degs[vtx]; d++) {
 if (induced->labels[intersect[d]] == k_idx) {
 out[count] = intersect[d];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 3;
 num_induced = induced->num_induced[k_idx-1];
 size_t ccurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t d=0; d < num_induced; d++) { induced->labels[prev_induced[d]] = k_idx; }
 for (size_t d=0; d < num_induced; ++d) {
 auto dvtx = prev_induced[d];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t e=0; e < induced->induced_degs[vtx]; e++) {
 if (induced->labels[intersect[e]] == k_idx) {
 out[count] = intersect[e];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 4;
 num_induced = induced->num_induced[k_idx-1];
 size_t dcurr_counts=0;
 if (num_induced != 0) {
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t e=0; e < num_induced; e++) { induced->labels[prev_induced[e]] = k_idx; }
 for (size_t e=0; e < num_induced; ++e) {
 auto evtx = prev_induced[e];
 intersect = induced->induced_edges + vtx * induced->nn;
 out = induced->induced + induced->nn * k_idx;
 count = 0;
 for (size_t f=0; f < induced->induced_degs[vtx]; f++) {
 if (induced->labels[intersect[f]] == k_idx) {
 out[count] = intersect[f];
 count++;
 }
 }
 induced->num_induced[k_idx] = count;
 if (induced->num_induced[k_idx] > k - k_idx - 1) {
 k_idx = 5;
 num_induced = induced->num_induced[k_idx-1];
 if (num_induced == 0) {continue;}
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 size_t ecurr_counts=0;
 counts = 0;
 for (size_t x=0; x < num_induced; i++) {
 vtx = prev_induced[x];
 intersect = induced->induced_edges + vtx * induced->nn;
 size_t tmp_counts = 0;
 for (size_t y=0; y < induced->induced_degs[vtx]; y++) {
 if (induced->labels[intersect[y]] == k_idx) {
 tmp_counts++;
 } 
 }
 etots[i] += tmp_counts;
 }
 for (size_t z=0; z < num_induced; z++) { induced->labels[prev_induced[z]] = k_idx - 1; }
 } }
 k_idx = 4;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t f=0; f<num_induced; f++) { induced->labels[prev_induced[f]] = k_idx - 1; }
  }
 } }
 k_idx = 3;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t e=0; e<num_induced; e++) { induced->labels[prev_induced[e]] = k_idx - 1; }
  }
 } }
 k_idx = 2;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t d=0; d<num_induced; d++) { induced->labels[prev_induced[d]] = k_idx - 1; }
  }
 } }
 k_idx = 1;
 num_induced = induced->num_induced[k_idx-1];
 prev_induced = induced->induced + induced->nn * (k_idx - 1);
 for (size_t c=0; c<num_induced; c++) { induced->labels[prev_induced[c]] = k_idx - 1; }
  }if (induced->use_base && tots[i] > 0) base_f(i, tots[i]);} else tots[i] = 0;}, 1, false);
 break; } 
 }

    double tt2 = t2.stop();
    std::cout << "##### Actual counting: " << tt2 << std::endl;

    return pbbslib::reduce_add(tots);
  }
} // namespace induced_neighborhood
