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
// numactl -i all ./LDD -beta 0.2 -rounds 3 -s -m com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -fa : run the fetch-and-add implementation of k-core
//     -nb : the number of buckets to use in the bucketing implementation

#include "LowDiameterDecomposition.h"

namespace gbbs {

template <class Graph>
double LDD_runner(Graph& G, commandLine P) {
  double beta = P.getOptionDoubleValue("-beta", 0.2);
  bool permute = P.getOption("-permute");
  std::cout << "### Application: LDD (Low-Diameter Decomposition)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -beta = " << beta << " -permute = " << permute
            << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  assert(P.getOption("-s"));
  timer t;
  t.start();
  auto ldd = LDD(G, beta, permute);
  double tt = t.stop();
  if (P.getOption("-stats")) {
    ldd_utils::num_clusters(ldd);
    ldd_utils::cluster_sizes(ldd);
    ldd_utils::num_intercluster_edges(G, ldd);
  }

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

generate_main(gbbs::LDD_runner, false);
