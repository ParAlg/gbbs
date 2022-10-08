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

#include "gbbs/gbbs.h"

#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

namespace gbbs {

template <class Graph>
struct countF {
  Graph& G;
  size_t* counts;
  countF(Graph& G, size_t* _counts) : G(G), counts(_counts) {}

  inline bool update(uintE s, uintE d) {
    auto d_neighbors = G.get_vertex(d).out_neighbors();
    gbbs::write_add(&counts[s], G.get_vertex(s).out_neighbors().intersect(
                                    &d_neighbors, s, d));
    return 1;
  }

  inline bool updateAtomic(uintE s, uintE d) {
    auto d_neighbors = G.get_vertex(d).out_neighbors();
    gbbs::write_add(&counts[s], G.get_vertex(s).out_neighbors().intersect(
                                    &d_neighbors, s, d));
    return 1;
  }
  inline bool cond(uintE d) { return cond_true(d); }
};

template <class Graph>
inline uintE* rankNodes(Graph& G, size_t n) {
  uintE* r = gbbs::new_array_no_init<uintE>(n);
  sequence<uintE> o = sequence<uintE>::uninitialized(n);

  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { o[i] = i; });
  parlay::sample_sort_inplace(make_slice(o), [&](const uintE u, const uintE v) {
    return G.get_vertex(u).out_degree() < G.get_vertex(v).out_degree();
  });
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { r[o[i]] = i; });
  return r;
}

// Directly call edgemap dense-forward.
template <class Graph, class VS, class F>
inline vertexSubset emdf(Graph& G, VS& vs, F f, const flags& fl = 0) {
  return edgeMapDenseForward<gbbs::empty>(G, vs, f, fl);
}

template <class Graph>
inline size_t CountDirected(Graph& DG, size_t* counts, vertexSubset& Frontier) {
  using W = typename Graph::weight_type;
  emdf(DG, Frontier, wrap_em_f<W>(countF<Graph>(DG, counts)), no_output);
  auto count_seq =
      parlay::delayed_seq<size_t>(DG.n, [&](size_t i) { return counts[i]; });
  size_t count = parlay::reduce(count_seq);
  return count;
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
inline size_t CountDirectedBalanced(Graph& DG, size_t* counts, const F& f) {
  using W = typename Graph::weight_type;
  debug(std::cout << "Starting counting"
                  << "\n";);
  size_t n = DG.n;

  auto parallel_work = sequence<size_t>::uninitialized(n);
  {
    auto map_f = [&](uintE u, uintE v, W wgh) -> size_t {
      return DG.get_vertex(v).out_degree();
    };
    parallel_for(0, n, [&](size_t i) {
      auto monoid = parlay::addm<size_t>();
      parallel_work[i] = DG.get_vertex(i).out_neighbors().reduce(map_f, monoid);
    });
  }
  size_t total_work = parlay::scan_inplace(make_slice(parallel_work));

  size_t block_size = 50000;
  size_t n_blocks = total_work / block_size + 1;
  size_t work_per_block = (total_work + n_blocks - 1) / n_blocks;
  std::cout << "Total work = " << total_work << " nblocks = " << n_blocks
            << " work per block = " << work_per_block << "\n";

  auto run_intersection = [&](size_t start_ind, size_t end_ind) {
    for (size_t i = start_ind; i < end_ind; i++) {  // check LEQ
      auto our_neighbors = DG.get_vertex(i).out_neighbors();
      size_t total_ct = 0;
      auto map_f = [&](uintE u, uintE v, W wgh) {
        auto their_neighbors = DG.get_vertex(v).out_neighbors();
        total_ct += our_neighbors.intersect_f_par(&their_neighbors, f);
      };
      our_neighbors.map(map_f, false);  // run map sequentially
      counts[i] = total_ct;
    }
  };

  parallel_for(0, n_blocks, 1, [&](size_t i) {
    size_t start = i * work_per_block;
    size_t end = (i + 1) * work_per_block;
    auto less_fn = std::less<size_t>();
    size_t start_ind = parlay::binary_search(parallel_work, start, less_fn);
    size_t end_ind = parlay::binary_search(parallel_work, end, less_fn);
    run_intersection(start_ind, end_ind);
  });

  auto count_seq = gbbs::make_slice<size_t>(counts, DG.n);
  size_t count = parlay::reduce(count_seq);

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
  auto counts = sequence<size_t>::uninitialized(n);
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { counts[i] = 0; });

  // 1. Rank vertices based on degree
  timer rt;
  rt.start();
  uintE* rank = rankNodes(G, G.n);
  rt.stop();
  rt.next("rank time");

  // 2. Direct edges to point from lower to higher rank vertices.
  // Note that we currently only store out-neighbors for this graph to save
  // memory.
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DG = filterGraph(G, pack_predicate);
  // auto DG = Graph::filterGraph(G, pack_predicate);
  gt.stop();
  gt.next("build graph time");

  // 3. Count triangles on the digraph
  timer ct;
  ct.start();

  size_t count = CountDirectedBalanced(DG, counts.begin(), f);
  std::cout << "### Num triangles = " << count << "\n";
  ct.stop();
  ct.next("count time");
  gbbs::free_array(rank, G.n);
  return count;
}

template <class Graph, class F, class O>
inline size_t Triangle_degeneracy_ordering(Graph& G, const F& f,
                                           O ordering_fn) {
  using W = typename Graph::weight_type;
  timer gt;
  gt.start();
  uintT n = G.n;
  auto counts = sequence<size_t>::uninitialized(n);
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { counts[i] = 0; });

  timer rt;
  rt.start();
  auto ordering = ordering_fn(G);
  rt.stop();
  rt.next("rank time");
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return (ordering[u] < ordering[v]);
  };

  auto DG = filterGraph(G, pack_predicate);
  // auto DG = Graph::filterGraph(G, pack_predicate);
  gt.stop();
  gt.next("build graph time");

  // 3. Count triangles on the digraph
  timer ct;
  ct.start();

  size_t count = CountDirectedBalanced(DG, counts.begin(), f);
  std::cout << "### Num triangles = " << count << "\n";
  ct.stop();
  ct.next("count time");
  return count;
}

template <class Graph, class F>
inline size_t Triangle(Graph& G, const F& f, const std::string& ordering,
                       commandLine& P) {
  if (ordering == "degree") {
    return Triangle_degree_ordering<Graph, F>(G, f);
  } else if (ordering == "goodrich") {
    auto eps = P.getOptionDoubleValue("-e", 0.1);
    auto ff = [&](Graph& graph) -> sequence<uintE> {
      return goodrichpszona_degen::DegeneracyOrder_intsort(graph, eps);
    };
    return Triangle_degeneracy_ordering<Graph, F>(G, f, ff);
  } else if (ordering == "kcore") {
    auto ff = [&](Graph& graph) -> sequence<uintE> {
      auto D = DegeneracyOrder(graph);
      auto ret = sequence<uintE>::from_function(
          graph.n, [&](size_t i) { return D[i]; });
      return ret;
    };
    return Triangle_degeneracy_ordering<Graph, F>(G, f, ff);
  } else if (ordering == "barenboimelkin") {
    auto ff = [&](Graph& graph) -> sequence<uintE> {
      return barenboimelkin_degen::DegeneracyOrder(graph);
    };
    return Triangle_degeneracy_ordering<Graph, F>(G, f, ff);
  } else {
    std::cerr << "Unexpected ordering: " << ordering << '\n';
    exit(1);
  }
}

}  // namespace gbbs
