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
#include "lib/sample_sort.h"
#include "ligra.h"

template <template <class W> class vertex, class W>
struct countF {
  using w_vertex = vertex<W>;
  w_vertex* V;
  size_t* counts;
  countF(w_vertex* _V, size_t* _counts) : V(_V), counts(_counts) {}

  inline bool update(uintE s, uintE d) {
    writeAdd(&counts[s], V[s].intersect(&V[d], s, d));
    return 1;
  }

  inline bool updateAtomic(uintE s, uintE d) {
    writeAdd(&counts[s], V[s].intersect(&V[d], s, d));
    return 1;
  }
  inline bool cond(uintE d) { return cond_true(d); }
};

template <class vertex>
inline uintE* rankNodes(vertex* V, size_t n) {
  uintE* r = pbbs::new_array_no_init<uintE>(n);
  uintE* o = pbbs::new_array_no_init<uintE>(n);

  timer t;
  t.start();
  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i) { o[i] = i; });
  pbbs::sample_sort(o, n, [&](const uintE u, const uintE v) {
    return V[u].getOutDegree() < V[v].getOutDegree();
  });
  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i)
                  { r[o[i]] = i; });
  t.stop();
  t.reportTotal("Rank time");
  pbbs::free_array(o);
  return r;
}

// Directly call edgemap dense-forward.
template <class vertex, class VS, class F>
inline vertexSubset emdf(graph<vertex> GA, VS& vs, F f, const flags& fl = 0) {
  return edgeMapDenseForward<pbbs::empty>(GA, vs, f, fl);
}

template <template <class W> class vertex, class W>
inline size_t CountDirected(graph<vertex<W>>& DG, size_t* counts,
                            vertexSubset& Frontier) {
  emdf(DG, Frontier, wrap_em_f<W>(countF<vertex, W>(DG.V, counts)), no_output);
  auto count_seq = sequence<size_t>(counts, DG.n);
  size_t count = pbbs::reduce_add(count_seq);
  return count;
}

template <template <class W> class vertex, class W, class F>
inline size_t CountDirectedBalanced(graph<vertex<W>>& DG, size_t* counts,
                                    const F& f) {
  std::cout << "Starting counting "
            << "\n";
  size_t n = DG.n;

  auto parallel_work = sequence<size_t>(n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      return DG.V[v].getOutDegree();
    };
    auto reduce_f = [&](const size_t& l, const size_t& r) { return l + r; };
    size_t id = 0;
    par_for(0, n, [&] (size_t i) {
      parallel_work[i] = DG.V[i].reduceOutNgh(i, id, map_f, reduce_f);
    });
  }
  size_t total_work = pbbs::scan_add(parallel_work, parallel_work);

  size_t n_blocks = nworkers() * 8 + 1;
  size_t work_per_block = total_work / n_blocks;
  std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
            << " work per block = " << work_per_block << "\n";

  auto V = DG.V;
  auto run_intersection = [&](size_t start_ind, size_t end_ind) {
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      auto vtx = V[i];
      size_t total_ct = 0;
      auto map_f = [&](uintE u, uintE v, W wgh) {
        total_ct += vtx.intersect_f(&V[v], u, v, f);
      };
      vtx.mapOutNgh(i, map_f, false);  // run map sequentially
      counts[i] = total_ct;
    }
  };

  par_for(0, n_blocks, 1, [&] (size_t i) {
    size_t start = i * work_per_block;
    size_t end = (i + 1) * work_per_block;
    auto less_fn = std::less<size_t>();
    size_t start_ind = pbbs::binary_search(parallel_work, start, less_fn);
    size_t end_ind = pbbs::binary_search(parallel_work, end, less_fn);
    run_intersection(start_ind, end_ind);
  });

  auto count_seq = sequence<size_t>(counts, DG.n);
  size_t count = pbbs::reduce_add(count_seq);

  return count;
}

template <template <class W> class vertex, class W, class F>
inline size_t Triangle(graph<vertex<W>>& GA, const F& f) {
  timer gt;
  gt.start();
  uintT n = GA.n;
  auto counts = sequence<size_t>(n);
  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i)
                  { counts[i] = 0; });

  // 1. Rank vertices based on degree
  uintE* rank = rankNodes(GA.V, GA.n);

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

  size_t count = CountDirectedBalanced(DG, counts.start(), f);
  std::cout << "Num triangles = " << count << "\n";
  DG.clear();
  ct.stop();
  ct.reportTotal("count time");
  pbbs::free_array(rank);
  return count;
}
