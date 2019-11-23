// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <math.h>

#include "ligra/bucket.h"
#include "ligra/edge_map_reduce.h"
#include "ligra/ligra.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "pbbslib/list_allocator.h"
#include "pbbslib/integer_sort.h"

#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

#include "intersect.h"
#include "induced_intersection.h"
#include "relabel.h"

#define SIMD_STATE 4

template <class Graph>
inline uintE* rankNodes(Graph& G, size_t n) {
  uintE* r = pbbslib::new_array_no_init<uintE>(n);
  sequence<uintE> o(n);

  timer t;
  t.start();
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) { o[i] = i; });
  pbbslib::sample_sort_inplace(o.slice(), [&](const uintE u, const uintE v) {
    return G.get_vertex(u).getOutDegree() < G.get_vertex(v).getOutDegree();
  });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { r[o[i]] = i; });
  t.stop();
  debug(t.reportTotal("Rank time"););
  return r;
}


template <class Graph>
inline size_t KCliqueDir_fast_orig_rec(Graph& DG, size_t k_idx, size_t k, FullSpace_orig_lw* induced) {
  size_t num_induced = induced->num_induced[k_idx-1];
  uintE* prev_induced = induced->induced + induced->nn * (k_idx - 1);
  uintE* prev_induced_degs = induced->induced_degs + induced->nn * (k_idx - 1);
  if (num_induced == 0) return 0;

  if (k_idx + 1 == k) return induced->num_edges[k_idx - 1];

  size_t total_ct = 0;
  for (size_t i=0; i < num_induced; ++i) {
    uintE idx = prev_induced[i];
    uintE* new_induced = induced->induced + induced->nn * k_idx;
    induced->num_induced[k_idx] = prev_induced_degs[idx];
    uintE new_num_induced = induced->num_induced[k_idx];
    //parallel_for(0, new_num_induced, [&] (size_t j) {new_induced[j] = induced->induced_edges[idx * induced->nn + j]; });
    for (size_t j=0; j < new_num_induced; j++) { new_induced[j] = induced->induced_edges[idx * induced->nn + j]; }
    //parallel_for(0, new_num_induced, [&] (size_t j){ induced->labels[new_induced[j]] = k_idx; });
    for (size_t j=0; j < new_num_induced; j++) { induced->labels[new_induced[j]] = k_idx; }
    uintE* new_induced_degs = induced->induced_degs + induced->nn * k_idx;
    //parallel_for(0, induced->nn, [&] (size_t j) { new_induced_degs[j] = 0; });
    for (size_t j=0; j < induced->nn; j++) { new_induced_degs[j] = 0; }

    for (size_t j=0; j < new_num_induced; j++) {
      uintE v_idx = new_induced[j];
      uintE v_deg = prev_induced_degs[v_idx];
      uintE* v_edges = induced->induced_edges + v_idx * induced->nn;
      size_t end = v_deg;
      for (size_t l=0; l < end; l++) {
        if (induced->labels[v_edges[l]] == k_idx) new_induced_degs[v_idx]++;
        else { // if (to_save)
          auto tmp = v_edges[l];
          v_edges[l--] = v_edges[--end];
          v_edges[end] = tmp;
        }
      }
    }

    /*parallel_for(0, new_num_induced, [&] (size_t j) {
      uintE v_idx = new_induced[j];
      uintE v_deg = prev_induced_degs[v_idx];
      uintE* v_edges = induced->induced_edges + v_idx * induced->nn;
      size_t end = v_deg;
      for (size_t l=0; l < end; l++) {
        if (induced->labels[v_edges[l]] == k_idx) new_induced_degs[v_idx]++;
        else { // if (to_save)
          auto tmp = v_edges[l];
          v_edges[l--] = v_edges[--end];
          v_edges[end] = tmp;
        }
      }
    });*/

    auto deg_seq = pbbslib::make_sequence(new_induced_degs, induced->nn);
    induced->num_edges[k_idx] = pbbslib::reduce_add(deg_seq);

    //uintE vtx = prev_induced[i];
    //induced->num_induced[k_idx] = lstintersect_set(prev_induced, num_induced, (uintE*)(DG.get_vertex(vtx).getOutNeighbors()), DG.get_vertex(vtx).getOutDegree(), true, induced->induced + induced->num_induced[0] * k_idx);
    if (induced->num_induced[k_idx] > 0) total_ct += KCliqueDir_fast_orig_rec(DG, k_idx + 1, k, induced);
    //parallel_for(0, new_num_induced, [&] (size_t j){ induced->labels[new_induced[j]] = k_idx-1; });
    for (size_t j=0; j < new_num_induced; j++) { induced->labels[new_induced[j]] = k_idx-1; }
  }

  return total_ct;
}

template <class Graph>
inline size_t KCliqueDir_fast_orig(Graph& DG, size_t k) {
  size_t n = 0;
  size_t max_deg = get_max_deg(DG);
  FullSpace_orig_lw* induced = nullptr;
  #pragma omp parallel private(induced) reduction(+:n)
  {
  induced = new FullSpace_orig_lw(max_deg, k);
  #pragma omp for schedule(dynamic, 1) nowait
  for (size_t i=0; i < DG.n; ++i) {
    if (DG.get_vertex(i).getOutDegree() != 0) {
      induced->setup(DG, k, i);
      n += KCliqueDir_fast_orig_rec(DG, 1, k, induced);
    }
  }

  if (induced != nullptr) { induced->del(); delete induced; }

  }

  return n;
}

template <class Graph>
pbbs::sequence<uintE> get_ordering(Graph& GA, long order_type, double epsilon = 0.1, uintE* rankfile = nullptr) {
  if (order_type == 0) return goodrichpszona_degen::DegeneracyOrder_intsort(GA, epsilon);
  else if (order_type == 1) return barenboimelkin_degen::DegeneracyOrder(GA, epsilon);
  else if (order_type == 2) {
    auto rank = sequence<uintE>(GA.n, [&](size_t i) { return i; });
    auto kcore = KCore(GA);
    auto get_core = [&](uintE& p) -> uintE { return kcore[p]; };
    pbbs::integer_sort_inplace(rank.slice(), get_core);
    return rank;
  }
  else if (order_type == 3) return pbbslib::make_sequence(rankNodes(GA, GA.n), GA.n);
  else if (order_type == 4) return pbbslib::make_sequence(rankfile, GA.n);
}


// induced
// generated
// -i 0 (simple gbbs intersect), -i 2 (simd intersect), -i 1 (graph set inter)
// -o 0 (goodrich), 1 (barnboimelkin approx), 2 (barenboimelkin exact)

// todo approx work and do some kind of break in gen if too much
// TODO get rid of duplicates in edge lists????
template <class Graph>
inline size_t KClique(Graph& GA, size_t k, long order_type = 0, double epsilon = 0.1,
long space_type = 2, uintE* rankfile = nullptr) {
  using W = typename Graph::weight_type;
  assert (k >= 1);
  if (k == 1) return GA.n;
  else if (k == 2) return GA.m;

  timer t_rank; t_rank.start();
  sequence<uintE> rank = get_ordering(GA, order_type, epsilon, rankfile);
  double tt_rank = t_rank.stop();
  std::cout << "### Rank Running Time: " << tt_rank << std::endl;

  timer t_filter; t_filter.start();
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return (rank[u] < rank[v]) && GA.get_vertex(u).getOutDegree() >= k-1 && GA.get_vertex(v).getOutDegree() >= k-1;
  };
  auto DG = relabel_graph(GA, rank.begin(), pack_predicate); //filter_graph(GA, pack_predicate);
  double tt_filter = t_filter.stop();
  std::cout << "### Filter Graph Running Time: " << tt_filter << std::endl;

  // Done preprocessing
  timer t; t.start();
  size_t count = 0;
  if (space_type == 2) {
    count = KCliqueDir_fast(DG, k-1);
  }
  else if (space_type == 3) {
    count = KCliqueDir_fast_orig(DG, k-1);
  }

  double tt = t.stop();
  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << k << " cliques = " << count << "\n";
  return count;
}
