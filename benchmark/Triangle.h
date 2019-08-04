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

#include <algorithm>
#include "pbbslib/sample_sort.h"
#include "pbbslib/monoid.h"
#include "ligra.h"

template <template <class W> class vertex, class W>
struct countF {
  using w_vertex = vertex<W>;
  w_vertex* V;
  size_t* counts;
  countF(w_vertex* _V, size_t* _counts) : V(_V), counts(_counts) {}

  inline bool update(uintE s, uintE d) {
    pbbslib::write_add(&counts[s], V[s].intersect(&V[d], s, d));
    return 1;
  }

  inline bool updateAtomic(uintE s, uintE d) {
    pbbslib::write_add(&counts[s], V[s].intersect(&V[d], s, d));
    return 1;
  }
  inline bool cond(uintE d) { return cond_true(d); }
};

template <class G>
inline uintE* rankNodes(G& GA, size_t n) {
  uintE* r = pbbslib::new_array_no_init<uintE>(n);
//  uintE* o = pbbslib::new_array_no_init<uintE>(n);
  sequence<uintE> o(n);

  timer t;
  t.start();
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) { o[i] = i; });
  pbbslib::sample_sort_inplace(o.slice(), [&](const uintE u, const uintE v) {
    return GA.get_vertex(u).getOutDegree() < GA.get_vertex(v).getOutDegree();
  });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { r[o[i]] = i; });
  t.stop();
  debug(t.reportTotal("Rank time"););
  return r;
}

template <class PG, class F>
inline size_t CountDirectedSlabs(PG& DG, size_t* counts, const F& f, size_t n_slabs=10) {
  using W = typename PG::weight_type;
  size_t n = DG.n;

  auto work_inefficiency = sequence<size_t>(n);
  auto parallel_work = sequence<size_t>(n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      size_t degree = DG.get_vertex(v).getOutDegree();
      assert(degree < n);
      return degree;
    };
    par_for(0, n, [&] (size_t i) {
      auto monoid = pbbslib::addm<size_t>();
      parallel_work[i] = DG.get_vertex(i).reduceOutNgh(i, map_f, monoid);
    });
  }
  size_t total_work = pbbslib::scan_add_inplace(parallel_work.slice());

  size_t block_size = 50000;
  size_t n_blocks = total_work/block_size + 1;
  size_t work_per_block = total_work / n_blocks;
  debug(std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
            << " work per block = " << work_per_block << "\n";);

  auto run_intersection = [&](size_t start_ind, size_t end_ind) {
    uintE stk[8192];
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      auto vtx = DG.get_vertex(i);
      size_t total_ct = 0;

      uintE* nghs = (uintE*)stk;
      uintE deg = vtx.getOutDegree();
      if (deg > 8192) {
        nghs = pbbs::new_array_no_init<uintE>(deg);
      }
      size_t k = 0;
      auto map_seq_f = [&] (const uintE& u, const uintE& w, const W& wgh) {
        nghs[k++] = w;
      };
      vtx.mapOutNgh(i, map_seq_f, false);
      auto our_seq = pbbslib::make_sequence(nghs, deg);

//      auto map_f = [&](uintE u, uintE v, W wgh) {
//        // Copy live neighbors of u into separate array?
//        auto ngh_vtx = DG.get_vertex(v);
////        uintE ngh_stk[8192];
////        uintE ngh_deg = ngh_vtx.getOutDegree();
////        size_t idx = 0;
////        auto map_ngh_f = [&] (const uintE& u, const uintE& w, const W& wgh) {
////          ngh_stk[idx++] = w;
////        };
////        ngh_vtx.mapOutNgh(v, map_ngh_f, false);
////        assert(ngh_deg == idx);
////        uintE* ngh_arr = (uintE*)ngh_stk;
////        auto ngh_seq = pbbslib::make_sequence(ngh_arr, ngh_deg);
//        auto iter = ngh_vtx.getOutIter();
//        total_ct += block_vertex_ops::intersect_seq(our_seq, iter);
////        total_ct += block_vertex_ops::seq_merge(our_seq, ngh_seq);
//      };
//      vtx.mapOutNgh(i, map_f, false);  // run map sequentially
//      counts[i] = total_ct;

      if (deg > 8192) {
        pbbs::free_array(nghs);
      }
    }
  };

  // How to calculate the work-inefficiency?


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




template <class PG, class F>
inline size_t CountDirectedBalanced(PG& DG, size_t* counts,
                                    const F& f) {
  using W = typename PG::weight_type;
  debug(std::cout << "Starting counting"
            << "\n";);
  size_t n = DG.n;

  auto work_inefficiency = sequence<size_t>(n);
  auto parallel_work = sequence<size_t>(n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      size_t degree = DG.get_vertex(v).getOutDegree();
      assert(degree < n);
      return degree;
    };
//    auto map_2_f = [&](uintE u, uintE v, W wgh) -> size_t {
//      auto ngh_vtx = DG.get_vertex(v);
//      size_t num_out_blocks = ngh_vtx.getNumOutBlocks();
//      if (num_out_blocks > 0) {
//        size_t last_block_id = num_out_blocks-1;
//        size_t last_block_size = ngh_vtx.out_block_degree(last_block_id);
//        if (num_out_blocks > 1) {
//          return (num_out_blocks - 1)*ngh_vtx.block_manager.edges_per_block + last_block_size;
//        } else {
//          return last_block_size;
//        }
//      }
//      return static_cast<size_t>(0);
//    };
    par_for(0, n, [&] (size_t i) {
      auto monoid = pbbslib::addm<size_t>();
      parallel_work[i] = DG.get_vertex(i).reduceOutNgh(i, map_f, monoid);
//      work_inefficiency[i] = DG.get_vertex(i).reduceOutNgh(i, map_2_f, monoid);
    });
  }
  size_t total_work = pbbslib::scan_add_inplace(parallel_work.slice());

//  size_t total_work_inefficiency = pbbslib::scan_add_inplace(work_inefficiency.slice());
//  cout << "total_work_inefficiency = " << total_work_inefficiency << endl;


  // TODO Can get even better load balance by using the actual blocks
  // corresponding to a vertex. (similar to edgeMapBlocked)

  size_t block_size = 50000;
  size_t n_blocks = total_work/block_size + 1;
//  size_t n_blocks = num_workers() * 8 + 1;
  size_t work_per_block = total_work / n_blocks;
  debug(std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
            << " work per block = " << work_per_block << "\n";);

  auto run_intersection = [&](size_t start_ind, size_t end_ind) {
    uintE stk[8192];
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      auto vtx = DG.get_vertex(i);
      size_t total_ct = 0;

      uintE* nghs = (uintE*)stk;
      uintE deg = vtx.getOutDegree();
      if (deg > 8192) {
        nghs = pbbs::new_array_no_init<uintE>(deg);
      }
      size_t k = 0;
      auto map_seq_f = [&] (const uintE& u, const uintE& w, const W& wgh) {
        nghs[k++] = w;
      };
      vtx.mapOutNgh(i, map_seq_f, false);

      auto our_seq = pbbslib::make_sequence(nghs, deg);

      auto map_f = [&](uintE u, uintE v, W wgh) {
        // Copy live neighbors of u into separate array?
        auto ngh_vtx = DG.get_vertex(v);
//        uintE ngh_stk[8192];
//        uintE ngh_deg = ngh_vtx.getOutDegree();
//        size_t idx = 0;
//        auto map_ngh_f = [&] (const uintE& u, const uintE& w, const W& wgh) {
//          ngh_stk[idx++] = w;
//        };
//        ngh_vtx.mapOutNgh(v, map_ngh_f, false);
//        assert(ngh_deg == idx);
//        uintE* ngh_arr = (uintE*)ngh_stk;
//        auto ngh_seq = pbbslib::make_sequence(ngh_arr, ngh_deg);
        auto iter = ngh_vtx.getOutIter();
        total_ct += block_vertex_ops::intersect_seq(our_seq, iter);
//        total_ct += block_vertex_ops::seq_merge(our_seq, ngh_seq);
      };
      vtx.mapOutNgh(i, map_f, false);  // run map sequentially
      counts[i] = total_ct;

      if (deg > 8192) {
        pbbs::free_array(nghs);
      }
    }
  };

  // How to calculate the work-inefficiency?


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



template <template <class W> class vertex, class W, class F>
inline size_t Triangle(symmetric_graph<vertex, W>& GA, const F& f) {
  timer gt;
  gt.start();
  uintT n = GA.n;
  auto counts = sequence<size_t>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { counts[i] = 0; });

  // 1. Rank vertices based on degree
  uintE* rank = rankNodes(GA, GA.n);

//  size_t xorr = 0;
//  for (size_t i=0; i<n; i++) {
//    xorr ^= GA.get_vertex(i).getOutDegree();
//  }
//  cout << "xorr = " << xorr << endl;
//  auto vtx_xors = pbbs::sequence<size_t>(n);
//  parallel_for(0, n, [&] (size_t v) {
//    auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
//      return v;
//    };
//    size_t id = 0;
//    auto reduce_f = [&] (const size_t& l, const size_t& r) { return l ^ r; };
//    auto reduce_m = pbbs::make_monoid(reduce_f, id);
//    vtx_xors[v] = GA.get_vertex(v).reduceOutNgh(v, map_f, reduce_m);
//  });
//  cout << "graph xor = " << pbbslib::reduce_xor(vtx_xors) << endl;

  // 2. Direct edges to point from lower to higher rank vertices.
  // Note that we currently only store out-neighbors for this graph to save
  // memory.
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DG = filter_graph<vertex, W>(GA, pack_predicate);
  gt.stop();
  gt.reportTotal("build graph time");

  // 3. Count triangles on the digraph
  timer ct;
  ct.start();

  auto get_degree = [&] (size_t i) { return DG.get_vertex(i).getOutDegree(); };
  auto degrees = pbbslib::make_sequence<size_t>(n, get_degree);
  cout << "max_degree = " << pbbslib::reduce_max(degrees) << endl;

  size_t count = CountDirectedSlabs(DG, counts.begin(), f);
  std::cout << "### Num triangles = " << count << "\n";
  DG.del();
  ct.stop();
  debug(ct.reportTotal("count time"););
  pbbslib::free_array(rank);
  return count;
}
