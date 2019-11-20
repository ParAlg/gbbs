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

#include "ligra/edge_map_reduce.h"
#include "ligra/bucket.h"
#include "ligra/ligra.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "pbbslib/list_allocator.h"
#include "pbbslib/integer_sort.h"

#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

enum OrderT {
  GoodrichPszonaOrder,
  BarenboimElkinOrder,
  KCoreOrder
};

template <class Graph, OrderT order>
auto get_ordering(Graph& GA, commandLine& P) -> pbbs::sequence<uintE> {
  if constexpr (order == GoodrichPszonaOrder) {
    return goodrichpszona_degen::DegeneracyOrder_intsort(GA, P);
  } else if constexpr (order == BarenboimElkinOrder) {
    return barenboimelkin_degen::DegeneracyOrder(GA, P);
  } else if constexpr (order == KCoreOrder) {
    auto dyn_arr = DegeneracyOrder(GA);
    auto ord = pbbs::sequence<uintE>(GA.n);
    parallel_for(0, GA.n, [&] (size_t i) {
      uintE u = dyn_arr.A[i];
      ord[u] = i;
    });
    dyn_arr.del();
    return ord;
  } else {
    assert(false);
  }
}



template <class Graph>
size_t CountDirectedBalanced(Graph& G, size_t k) {

}

template <class Graph, class W>
struct countF {
  Graph& G;
  size_t* counts;
  size_t k;
  countF(Graph& G, size_t* _counts, size_t k) : G(G), counts(_counts), k(k) {}

  inline size_t count_four_cliques_seq(uintE s, uintE d) {
    auto s_vtx = G.get_vertex(s);
    auto d_vtx = G.get_vertex(d);
    uintE stk[3000];
    if (s_vtx.getOutDegree() > 3000) {
      std::cout << "degree = " << s_vtx.getOutDegree() << std::endl;
      exit(-1);
    }
    uintE* ret = (uintE*)stk;
    size_t retct = s_vtx.intersect_ret(&d_vtx, s, d, ret);

    // ret now contains L = N_s \cap N_d
    // a. for each w \in L, intersect L and N_w

    size_t ct = 0;
    for (size_t k=0; k<retct; k++) {
      uintE w = ret[k];
      auto w_vtx = G.get_vertex(w);
      auto out_nghs = w_vtx.getOutNeighbors();
      auto left_seq = pbbslib::make_sequence<uintE>(w_vtx.getOutDegree(), [&] (size_t i) { return std::get<0>(out_nghs[i]); });
      auto right_seq = pbbslib::make_sequence<uintE>(retct, [&] (size_t i) { return ret[i]; });
      auto f = [&] (const uintE& u) {};
      ct += intersection::merge(left_seq, right_seq, f);
    }

    pbbslib::write_add(&counts[s], ct);
  }

  inline bool update(const uintE& s, const uintE& d, const W& w) {
    count_four_cliques_seq(s, d);
    return 1;
  }

  inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) {
    update(s, d, w);
    return 1;
  }
  inline bool cond(uintE d) { return cond_true(d); }
};

// Directly call edgemap dense-forward.
template <class Graph, class VS, class F>
inline vertexSubset emdf(Graph& G, VS& vs, F f, const flags& fl = 0) {
  return edgeMapDenseForward<pbbslib::empty>(G, vs, f, fl);
}

template <class Graph>
size_t CountDirected(Graph& G, size_t k) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  auto counts = pbbs::sequence<size_t>(n, (size_t)0);
  bool* frontier = pbbs::new_array_no_init<bool>(n);
  parallel_for(0, n, [&] (size_t i) { frontier[i] = true; });
  auto Frontier = vertexSubset(n, n, frontier);
  emdf(G, Frontier, countF<Graph, W>(G, counts.begin(), k), no_output);
  size_t count = pbbslib::reduce_add(counts);
  return count;
}

template <class Graph>
void ordering_stats(Graph& GA, pbbs::sequence<uintE>& order) {
  size_t n = GA.n;
  auto degree_seq = pbbslib::make_sequence<uintE>(n, [&] (size_t i) { return GA.get_vertex(i).getOutDegree(); });
  uintE max_degree = pbbslib::reduce_max(degree_seq);
  double avg_degree = static_cast<double>(GA.m)/GA.n;
  // median degree
  uintE med_degree = pbbs::approximate_kth_smallest(degree_seq, n/2, std::less<uintE>());

  std::cout << "avg_degree = " << avg_degree << " max_degree = " << max_degree << " med_degree = " << med_degree << std::endl;
}

template <class Graph, OrderT order>
size_t KClique(Graph& GA, commandLine& P, size_t k) {
  using W = typename Graph::weight_type;
  timer order_t; order_t.start();
  auto ordering = get_ordering<Graph, order>(GA, P);
  order_t.stop(); order_t.reportTotal("ordering time");

  timer t_filter; t_filter.start();

  // Orient edges s.t. they point upward in the order. Predicate returning true
  // => edge is preserved in the resulting directed graph.
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return (GA.get_vertex(u).getOutDegree() >= k-1) && // vertices with deg < k-1 can't participate in a k-clique: prune
           (GA.get_vertex(v).getOutDegree() >= k-1) && // ditto
           (ordering[u] < ordering[v]);
  };
  auto DG = filter_graph(GA, pack_predicate);
  double tt_filter = t_filter.stop();
  std::cout << "### Filter Graph Running Time: " << tt_filter << std::endl;
  ordering_stats(DG, ordering);

  timer count_t; count_t.start();
  size_t count = CountDirected(DG, k);
  count_t.stop(); count_t.reportTotal("count time");
  return count;
}
