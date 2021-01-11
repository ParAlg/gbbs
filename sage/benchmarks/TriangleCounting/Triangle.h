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

#include <algorithm>

#include "pbbslib/sample_sort.h"
#include "pbbslib/monoid.h"
#include "gbbs/gbbs.h"
#include "sage/sage.h"

namespace gbbs {

template <class Graph>
inline uintE* rankNodes(Graph& G, size_t n) {
  uintE* r = pbbslib::new_array_no_init<uintE>(n);
  sequence<uintE> o(n);

  timer t;
  t.start();
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) { o[i] = i; });
  pbbslib::sample_sort_inplace(o.slice(), [&](const uintE u, const uintE v) {
    return G.get_vertex(u).out_degree() < G.get_vertex(v).out_degree();
  });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { r[o[i]] = i; });
  t.stop();
  t.reportTotal("Rank time");
  return r;
}


// Returns the number of directed triangles in the input graph of the following
// orientation:
//        w
//       ^ ^
//      /   \.
//     u --> v
//
// Arguments:
//   DG
//     Graph on which we'll count triangles.
//   f: (uintE, uintE, uintE) -> void
//     Function that's run each triangle. On a directed triangle like the one
//     pictured above, we run `f(u, v, w)`.
template <class Graph, class F>
inline size_t CountDirectedBalanced(Graph& DG, size_t* counts,
                                    const F& f) {
  using W = typename Graph::weight_type;
  debug(std::cout << "Starting counting"
            << "\n";);
  size_t n = DG.n;

  auto parallel_work = sequence<size_t>(n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      return DG.get_vertex(v).out_degree();
    };
    par_for(0, n, [&] (size_t i) {
      auto monoid = pbbslib::addm<size_t>();
      parallel_work[i] = DG.get_vertex(i).out_neighbors().reduce(map_f, monoid);
    });
  }
  size_t total_work = pbbslib::scan_add_inplace(parallel_work.slice());

  size_t block_size = 50000;
  size_t n_blocks = total_work/block_size + 1;
  size_t work_per_block = (total_work + n_blocks - 1) / n_blocks;
  std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
            << " work per block = " << work_per_block << "\n";

  auto run_intersection = [&](size_t start_ind, size_t end_ind) {
    uintE stk[8192];
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      auto vtx = DG.get_vertex(i);
      size_t total_ct = 0;

      uintE* nghs = (uintE*)stk;
      uintE deg = vtx.out_degree();
      if (deg > 8192) {
        nghs = pbbs::new_array_no_init<uintE>(deg);
      }
      size_t k = 0;
      auto map_seq_f = [&] (const uintE& u, const uintE& w, const W& wgh) {
        nghs[k++] = w;
      };
      vtx.out_neighbors().map(map_seq_f, false);

      auto our_seq = pbbslib::make_sequence(nghs, deg);

      auto map_f = [&](uintE u, uintE v, W wgh) {
        // Copy live neighbors of u into separate array?
        auto ngh_vtx = DG.get_vertex(v);
//        uintE ngh_stk[8192];
//        uintE ngh_deg = ngh_vtx.out_degree();
//        size_t idx = 0;
//        auto map_ngh_f = [&] (const uintE& u, const uintE& w, const W& wgh) {
//          ngh_stk[idx++] = w;
//        };
//        ngh_vtx.mapOutNgh(v, map_ngh_f, false);
//        assert(ngh_deg == idx);
//        uintE* ngh_arr = (uintE*)ngh_stk;
//        auto ngh_seq = pbbslib::make_sequence(ngh_arr, ngh_deg);
        auto iter = ngh_vtx.out_neighbors().get_iter();
        total_ct += block_vertex_ops::intersect_seq(our_seq, iter);
//        total_ct += block_vertex_ops::seq_merge(our_seq, ngh_seq);
      };
      vtx.out_neighbors().map(map_f, false);  // run map sequentially
      counts[i] = total_ct;

      if (deg > 8192) {
        pbbs::free_array(nghs);
      }
    }
  };


  par_for(0, n_blocks, 1, [&] (size_t i) {
    size_t start = i * work_per_block;
    size_t end = (i + 1) * work_per_block;
    auto less_fn = std::less<size_t>();
    size_t start_ind = pbbslib::binary_search(parallel_work, start, less_fn);
    size_t end_ind = pbbslib::binary_search(parallel_work, end, less_fn);
    run_intersection(start_ind, end_ind);
  });

  auto count_seq = pbbslib::make_sequence<size_t>(counts, DG.n);
  size_t count = pbbslib::reduce_add(count_seq);

  return count;
}

// Counts the number of triangles in the input graph.
//
// Implementation note: this converts the input graph to a directed graph in
// which we point edges from lower-degree vertices to higher-degree vertices,
// hence the function name.
//
// Arguments:
//   G
//     Graph on which we'll count triangles.
//   f: (uintE, uintE, uintE) -> void
//     Function that's run each triangle. On a triangle with vertices {u, v, w},
//     we run `f(u, v, w)`.
//
// Returns:
//   The number of triangles in `G`.
template <class Graph, class F>
inline size_t Triangle_degree_ordering(Graph& G, const F& f) {
  using W = typename Graph::weight_type;
  timer gt;
  gt.start();
  uintT n = G.n;
  auto counts = sequence<size_t>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { counts[i] = 0; });

  // 1. Rank vertices based on degree
  uintE* rank = rankNodes(G, G.n);

  // 2. Direct edges to point from lower to higher rank vertices.
  // Note that we currently only store out-neighbors for this graph to save
  // memory.
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DG = sage::filter_graph(G, pack_predicate);
  gt.stop();
  gt.reportTotal("build graph time");

  // 3. Count triangles on the digraph
  timer ct;
  ct.start();

  size_t count = CountDirectedBalanced(DG, counts.begin(), f);
  std::cout << "### Num triangles = " << count << "\n";
  DG.del();
  ct.stop();
  ct.reportTotal("count time");
  pbbslib::free_array(rank);
  return count;
}

template <class Graph, class F>
inline size_t Triangle(Graph& G, const F& f, const std::string& ordering, commandLine& P) {
  if (ordering == "degree") {
    return Triangle_degree_ordering<Graph, F>(G, f);
  } else {
    std::cerr << "Unexpected ordering: " << ordering << '\n';
    exit(1);
  }
}

}  // namespace gbbs
