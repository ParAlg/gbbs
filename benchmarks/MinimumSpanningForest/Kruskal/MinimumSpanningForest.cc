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
// numactl -i all ./BellmanFord -src 10012 -s -m -rounds 3 twitter_wgh_SJ
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//     -w: indicate that the graph is weighted
//   optional:
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd

#include "MinimumSpanningForest.h"

namespace gbbs {

template <template <class W> class vertex, class W>
double MinimumSpanningForest_runner(symmetric_graph<vertex, W>& GA,
                                    commandLine P) {
  std::cout << "### Application: MinimumSpanningForest (Kruskal)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer mst_t;
  mst_t.start();
  MinimumSpanningForest_kruskal::MinimumSpanningForest(GA);
  double tt = mst_t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

generate_symmetric_weighted_main(gbbs::MinimumSpanningForest_runner, true);
