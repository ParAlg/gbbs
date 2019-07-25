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
inline size_t CountDirectedBalanced(PG& DG, size_t* counts,
                                    const F& f) {
  using W = typename PG::weight_type;
  debug(std::cout << "Starting counting"
            << "\n";);
  size_t n = DG.n;

  auto parallel_work = sequence<size_t>(n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      return DG.get_vertex(v).getOutDegree();
    };
    par_for(0, n, [&] (size_t i) {
      auto monoid = pbbslib::addm<size_t>();
      parallel_work[i] = DG.get_vertex(i).reduceOutNgh(i, map_f, monoid);
    });
  }
  size_t total_work = pbbslib::scan_add_inplace(parallel_work.slice());

  size_t block_size = 50000;
  size_t n_blocks = total_work/block_size + 1;
//  size_t n_blocks = num_workers() * 8 + 1;
  size_t work_per_block = total_work / n_blocks;
  debug(std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
            << " work per block = " << work_per_block << "\n";);

  auto run_intersection = [&](size_t start_ind, size_t end_ind) {
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      auto vtx = DG.get_vertex(i);
      size_t total_ct = 0;
      auto map_f = [&](uintE u, uintE v, W wgh) {
        auto ngh_vtx = DG.get_vertex(v);
        total_ct += vtx.intersect(ngh_vtx);
//        total_ct += vtx.intersect_f_par(&V[v], u, v, f);
      };
      vtx.mapOutNgh(i, map_f, false);  // run map sequentially
      counts[i] = total_ct;
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

template <template <class W> class vertex, class W, class F>
inline size_t Triangle(graph<vertex, W>& GA, const F& f) {
  timer gt;
  gt.start();
  uintT n = GA.n;
  auto counts = sequence<size_t>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { counts[i] = 0; });

  // 1. Rank vertices based on degree
  uintE* rank = rankNodes(GA, GA.n);

  // 2. Direct edges to point from lower to higher rank vertices.
  // Note that we currently only store out-neighbors for this graph to save
  // memory.
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DG = filter_graph<vertex, W>(GA, pack_predicate);
  gt.stop();
  debug(gt.reportTotal("build graph time"););

//  auto it = DG.get_vertex(10012).getOutIter();
//  cout << std::get<0>(it.cur()) << endl;
//  size_t kk = 1;
//  while (it.has_next()) {
//    cout << std::get<0>(it.next()) << endl;
//    kk++;
//  }
//  cout << "kk == " << kk << " vtx degree = " << it.degree() << endl;
//
//  cout << "DG.degree = " << DG.get_vertex(10012).getOutDegree() << endl;
//  auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
//    cout << "v = " << v << endl;
//  };
//  DG.get_vertex(10012).mapOutNgh(10012, map_f, false);
//  exit(0);

//  auto v_10012 = DG.get_vertex(10012);
//  auto v_10013 = DG.get_vertex(10013);
//  cout << "intersect = " << v_10012.intersect(v_10013) << endl;
//  exit(0);

  // 3. Count triangles on the digraph
  timer ct;
  ct.start();

  size_t count = CountDirectedBalanced(DG, counts.begin(), f);
  std::cout << "### Num triangles = " << count << "\n";
  DG.del();
  ct.stop();
  debug(ct.reportTotal("count time"););
  pbbslib::free_array(rank);
  return count;
}
