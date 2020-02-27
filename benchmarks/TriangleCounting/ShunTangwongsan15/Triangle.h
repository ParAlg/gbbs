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
#include "ligra/ligra.h"

#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

template <class Graph>
struct countF {
  Graph& G;
  size_t* counts;
  countF(Graph& G, size_t* _counts) : G(G), counts(_counts) {}

  inline bool update(uintE s, uintE d) {
    auto d_vertex = G.get_vertex(d);
    pbbslib::write_add(&counts[s], G.get_vertex(s).intersect(&d_vertex, s, d));
    return 1;
  }

  inline bool updateAtomic(uintE s, uintE d) {
    auto d_vertex = G.get_vertex(d);
    pbbslib::write_add(&counts[s], G.get_vertex(s).intersect(&d_vertex, s, d));
    return 1;
  }
  inline bool cond(uintE d) { return cond_true(d); }
};

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
  t.reportTotal("Rank time");
  return r;
}

// Directly call edgemap dense-forward.
template <class Graph, class VS, class F>
inline vertexSubset emdf(Graph& G, VS& vs, F f, const flags& fl = 0) {
  return edgeMapDenseForward<pbbslib::empty>(G, vs, f, fl);
}

template <class Graph>
inline size_t CountDirected(Graph& DG, size_t* counts,
                            vertexSubset& Frontier) {
  using W = typename Graph::weight_type;
  emdf(DG, Frontier, wrap_em_f<W>(countF<Graph>(DG, counts)), no_output);
  auto count_seq = sequence<size_t>(counts, DG.n);
  size_t count = pbbslib::reduce_add(count_seq);
  return count;
}

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
      return DG.get_vertex(v).getOutDegree();
    };
    par_for(0, n, [&] (size_t i) {
      auto monoid = pbbslib::addm<size_t>();
      parallel_work[i] = DG.get_vertex(i).template reduceOutNgh<size_t>(i, map_f, monoid);
    });
  }
  size_t total_work = pbbslib::scan_add_inplace(parallel_work.slice());

  size_t block_size = 50000;
  size_t n_blocks = total_work/block_size + 1;
  size_t work_per_block = total_work / n_blocks;
  n_blocks = (total_work/work_per_block) + 1;
  std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
            << " work per block = " << work_per_block << "\n";

  auto run_intersection = [&](size_t start_ind, size_t end_ind) {
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      auto vtx = DG.get_vertex(i);
      size_t total_ct = 0;
      auto map_f = [&](uintE u, uintE v, W wgh) {
        auto v_vtx = DG.get_vertex(v);
        total_ct += vtx.intersect_f_par(&v_vtx, u, v, f);
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
  auto DG = filter_graph(G, pack_predicate);
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

template <class Graph, class F, class O>
inline size_t Triangle_degeneracy_ordering(Graph& G, const F& f, O ordering_fn) {
  using W = typename Graph::weight_type;
  timer gt;
  gt.start();
  uintT n = G.n;
  auto counts = sequence<size_t>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { counts[i] = 0; });

  timer rt; rt.start();
  auto ordering = ordering_fn(G);
  rt.stop(); rt.reportTotal("rank time");
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return (ordering[u] < ordering[v]);
  };

  auto DG = filter_graph(G, pack_predicate);
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
  return count;
}

template <class Graph, class F>
inline size_t Triangle(Graph& G, const F& f, const std::string& ordering, commandLine& P) {
  if (ordering == "degree") {
    return Triangle_degree_ordering<Graph, F>(G, f);
  } else if (ordering == "goodrich") {
    auto eps = P.getOptionDoubleValue("-e", 0.1);
    auto ff = [&] (Graph& graph) -> pbbs::sequence<uintE> {
      return goodrichpszona_degen::DegeneracyOrder_intsort(graph, eps);
    };
    return Triangle_degeneracy_ordering<Graph, F>(G, f, ff);
  } else if (ordering == "kcore") {
    auto ff = [&] (Graph& graph) -> pbbs::sequence<uintE> {
      auto dyn_arr = DegeneracyOrder(graph);
      auto arr = dyn_arr.A; dyn_arr.A = nullptr;
      dyn_arr.alloc = false;
      return pbbs::sequence<uintE>(arr, graph.n);
    };
    return Triangle_degeneracy_ordering<Graph, F>(G, f, ff);
  } else if (ordering == "barenboimelkin") {
    auto ff = [&] (Graph& graph) -> pbbs::sequence<uintE> {
      return barenboimelkin_degen::DegeneracyOrder(graph);
    };
    return Triangle_degeneracy_ordering<Graph, F>(G, f, ff);
  } else {
    std::cerr << "Unexpected ordering: " << ordering << '\n';
    exit(1);
  }
}
