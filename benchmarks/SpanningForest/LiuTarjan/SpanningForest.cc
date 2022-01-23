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

#include "SpanningForest.h"
#include "benchmarks/SpanningForest/BFSSF/SpanningForest.h"
#include "benchmarks/SpanningForest/check.h"
#include "gbbs/gbbs.h"

namespace gbbs {
template <class Graph>
double SF_runner(Graph& G, commandLine P) {
  std::cout << "### Application: SpanningForest (LiuTarjan-based)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t;
  t.start();
  auto n = G.n;
  auto parents =
      sequence<parent>::from_function(n, [&](size_t i) { return i; });
  auto edges = sequence<edge>(n, empty_edge);

  auto opt_connect = lt::template get_connect_function<parent_connect>();
  auto opt_update = lt::template get_update_function<root_update>();
  auto opt_shortcut = lt::template get_shortcut_function<full_shortcut>();

  auto alg =
      lt::LiuTarjanAlgorithm<decltype(opt_connect), parent_connect,
                             decltype(opt_update), root_update,
                             decltype(opt_shortcut), full_shortcut, Graph>(
          G, n, opt_connect, opt_update, opt_shortcut);

  alg.initialize(parents, edges);

  alg.template compute_spanning_forest<no_sampling>(parents, edges);

  auto filtered_edges = parlay::filter(
      make_slice(edges), [&](const edge& e) { return e != empty_edge; });

  std::cout << "sf has: " << filtered_edges.size() << " many edges"
            << std::endl;
  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;

  if (P.getOptionValue("-check")) {
    auto bfs_edges = bfs_sf::SpanningForestDet(G);
    spanning_forest::check_spanning_forest(G.n, bfs_edges, filtered_edges);
  }

  return tt;
}
}  // namespace gbbs

generate_main(gbbs::SF_runner, false);
