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

// Usage Example:
// numactl -i all ./MaximalMatching -s -c -m clueweb_sym.bytepda
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed

#include "MaximalMatching.h"

#include "gbbs/gbbs.h"

#include <fstream>
#include <iostream>

namespace gbbs {

template <template <class W> class vertex, class W>
double MaximalMatching_runner(symmetric_graph<vertex, W>& G, commandLine P) {
  std::cout << "### Application: Maximal Matching" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: (n/a)" << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  assert(P.getOption("-s"));  // input graph must be symmetric

  timer t;
  t.start();
  auto matching = MaximalMatching(G);
  double tt = t.stop();

  verify_matching(G, matching);

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

generate_symmetric_main(gbbs::MaximalMatching_runner, false);
