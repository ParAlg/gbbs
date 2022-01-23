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
// > numactl -i all ./Coloring -s -m clueweb_sym.bytepda
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//   optional:
//     -c : indicate that the graph should be mmap'd
//     -m : indicate that the graph is compressed
//     -lf : use the LF (largest degree first) herustic
//     -stats : output statistics on the resulting coloring
//     -verify : verify that the algorithm produced a valid coloring
//
// The default heuristic used is LLF, which provably achieves polynomial
// parallelism (see "Ordering Heuristics for Parallel Graph Coloring" by
// Hasenplaugh et al.)

#include "GraphColoring.h"

#include <fstream>
#include <iostream>

namespace gbbs {
template <class Graph>
double Coloring_runner(Graph& G, commandLine P) {
  bool runLF = P.getOption("-lf");
  std::cout << "### Application: Coloring" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -lf = " << runLF << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t;
  t.start();
  auto colors = Coloring(G, runLF);
  double tt = t.stop();
  if (P.getOption("-stats")) {
    std::cout << "num_colors = " << parlay::reduce_max(colors) << "\n";
  }
  if (P.getOption("-verify)")) {
    verify_coloring(G, colors);
  }

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}
}  // namespace gbbs

generate_main(gbbs::Coloring_runner, false);
