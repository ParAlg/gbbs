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

// Usage:
// numactl -i all ./SpanningForest -rounds 3 -s -m twitter_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc

#include "SpanningForest.h"
#include "benchmarks/SpanningForest/BFSSF/SpanningForest.h"
#include "benchmarks/SpanningForest/check.h"

namespace gbbs {
template <class Graph>
double SpanningForest_runner(Graph& G, commandLine P) {
  auto beta = P.getOptionDoubleValue("-beta", 0.2);
  std::cout << "### Application: SpanningForest" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -beta = " << beta << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  auto pack = P.getOption("-pack");
  assert(P.getOption("-s"));
  assert(!pack);  // discouraged for now. Using the optimized contraction method
                  // is faster.
  timer t;
  t.start();
  auto edges = workefficient_sf::SpanningForest(G, beta, pack,
                                                P.getOptionValue("-permute"));
  std::cout << "n = " << G.n << " #edges = " << edges.size() << std::endl;
  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;

  std::cout << "vtx 0 has degree: " << G.get_vertex(0).out_degree()
            << std::endl;

  if (P.getOptionValue("-check")) {
    auto bfs_edges = bfs_sf::SpanningForestDet(G);
    spanning_forest::check_spanning_forest(G.n, bfs_edges, edges);
  }

  if (pack) {
    // packing mutates the graph, packing out all intra-cluster edges, and can
    // only be run once unless the input graph is copied.
    exit(0);
  }
  return tt;
}
}  // namespace gbbs

generate_main(gbbs::SpanningForest_runner, false);
