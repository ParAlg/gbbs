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
// numactl -i all ./SSBetweennessCentrality -src 10012 -s -m -rounds 3
// twitter_SJ
// flags:
//   required:
//     -src: the source to compute centrality contributions from
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#include "SSBetweennessCentrality.h"

namespace gbbs {

template <class Graph>
double SSBetweennessCentrality_runner(Graph& G, commandLine P) {
  uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
  std::cout << "### Application: SSBetweennessCentrality" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -src = " << src << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t;
  t.start();
  if (P.getOptionValue("-fa")) {
    auto scores = bc::SSBetweennessCentrality_EM(G, src);
    for (size_t i = 0; i < 100; i++) {
      std::cout << scores[i] << std::endl;
    }
  } else if (P.getOptionValue("-ligra")) {
    auto scores = bc::SSBetweennessCentrality(G, src);
    for (size_t i = 0; i < 100; i++) {
      std::cout << scores[i] << std::endl;
    }
  } else {
    /* Default: no contention --- reduceNgh technique */
    auto scores = bc_bfs::SSBetweennessCentrality_BFS(G, src);
    for (size_t i = 0; i < 100; i++) {
      std::cout << scores[i] << std::endl;
    }
  }

  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;

  return tt;
}

}  // namespace gbbs

generate_main(gbbs::SSBetweennessCentrality_runner, false);
