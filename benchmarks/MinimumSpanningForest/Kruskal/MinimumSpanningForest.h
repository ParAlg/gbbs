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

#include <cassert>
#include "gbbs/gbbs.h"

#include "benchmarks/Connectivity/UnionFind/union_find_rules.h"
#include "benchmarks/Connectivity/common.h"

namespace gbbs {
namespace MinimumSpanningForest_kruskal {

template <template <class W> class vertex, class W,
          typename std::enable_if<!std::is_same<W, gbbs::empty>::value,
                                  int>::type = 0>
inline sequence<std::tuple<uintE, uintE, W>> MinimumSpanningForest(
    symmetric_graph<vertex, W>& GA) {
  using edge = std::tuple<uintE, uintE, W>;

  size_t n = GA.n;
  sequence<edge> edges = GA.edges();

  timer kt;
  kt.start();
  auto weight_seq = parlay::delayed_seq<W>(
      edges.size(), [&](size_t i) { return std::get<2>(edges[i]); });
  std::cout << weight_seq[0] << std::endl;
  auto max_weight = parlay::reduce_max(weight_seq);
  std::cout << "max_weight = " << max_weight << std::endl;
  timer st;
  st.start();
  auto comp = [&](const edge& l, const edge& r) {
    return std::get<2>(l) < std::get<2>(r);
  };
  parlay::sample_sort_inplace(make_slice(edges), comp);
  st.stop();
  st.next("sort time");

  auto components =
      sequence<uintE>::from_function(n, [&](size_t i) { return i; });
  constexpr auto find{find_variants::find_compress};
  auto unite{unite_variants::Unite<decltype(find)>{find}};

  W weight = 0;
  auto MST = sequence<edge>(n);
  size_t k = 0;
  for (size_t i = 0; i < edges.size(); i++) {
    auto[u, v, w] = edges[i];
    if (find(u, components) != find(v, components)) {
      unite(u, v, components);
      weight += w;
      MST[k++] = {u, v, w};
    }
  }
  kt.stop();
  kt.next(
      "kruskal time (excluding G.get_edges() to convert graph to edge-list "
      "format)");
  std::cout << "MST weight = " << weight << std::endl;
  MST.resize(k);
  return MST;
}

template <
    template <class W> class vertex, class W,
    typename std::enable_if<std::is_same<W, gbbs::empty>::value, int>::type = 0>
inline sequence<std::pair<uintE, uintE>> MinimumSpanningForest(
    symmetric_graph<vertex, W>& GA) {
  return sequence<std::pair<uintE, uintE>>();
}

}  // namespace MinimumSpanningForest_kruskal
}  // namespace gbbs
