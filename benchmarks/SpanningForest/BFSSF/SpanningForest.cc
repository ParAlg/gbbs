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
#include "ligra/ligra.h"
//#include "benchmarks/SpanningForest/common.h"
#include "benchmarks/SpanningForest/check.h"

template <class Graph>
double SpanningForest_runner(Graph& G, commandLine P) {
  std::cout << "### Application: SpanningForest (BFS-based)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### ------------------------------------" << endl;
  assert(P.getOption("-s"));
  timer t;
  t.start();
  auto edges_nd = bfs_sf::SpanningForest(G);
  auto edges_det = bfs_sf::SpanningForestDet(G);
  double tt = t.stop();

  size_t cor_xor = 0; size_t check_xor = 0;
  for (size_t i=0; i<edges_nd.size(); i++) {
    auto [u, v] = edges_nd[i];
    auto [up, vp] = edges_det[i];
    cor_xor ^= ((static_cast<size_t>(u) << 32L) + static_cast<size_t>(v));
    check_xor ^= ((static_cast<size_t>(up) << 32L) + static_cast<size_t>(vp));
  }
  std::cout << "correct xor = " << cor_xor << " check_xor = " << check_xor << std::endl;

  spanning_forest::check_spanning_forest(G.n, edges_nd, edges_det);
//  auto edges = P.getOption("-det") ? bfs_sf::SpanningForestDet(G) : bfs_sf::SpanningForest(G);
//  std::cout << "### Running Time: " << tt << std::endl;
//
//  if (P.getOptionValue("-check")) {
//    spanning_forest::sf_compute_and_check(G, edges);
//  }

  return tt;
}

generate_main(SpanningForest_runner, false);
