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
// numactl -i all ./KCore -rounds 3 -s -m com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -fa : run the fetch-and-add implementation of k-core
//     -nb : the number of buckets to use in the bucketing implementation

#include "TwoThreeNucleusDecomposition.h"
#include <math.h>
#include <fstream>
namespace gbbs {
template <class Graph>
double AppTwoThreeNucleusDecomposition_runner(Graph& GA, commandLine P) {
  long r = P.getOptionLongValue("-r", 2); // k as in k-cliques
  long ss = P.getOptionLongValue("-ss", 3); // k as in k-cliques
  long table_type = P.getOptionLongValue("-tt", 3); // 1 = 1 lvl, 2 = 2 lvls, 3 = multi
  long num_levels = P.getOptionLongValue("-nl", 2); // only for multi, # levels
  bool relabel = P.getOptionValue("-relabel"); // for true, relabel graph
  bool contiguous_space = P.getOptionValue("-contig"); // for true, contiguous space
  bool verify = P.getOptionValue("-verify"); // for testing
  long efficient = P.getOptionLongValue("-efficient", 1); // for list buffer
  bool use_compress = P.getOptionValue("-compress"); //only for 2, 3

  std::cout << "### Application: Nucleus Decomposition" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  assert(P.getOption("-s"));

  timer t; t.start();

  NucleusDecomposition(GA, r, ss, table_type, num_levels, relabel, contiguous_space, verify, efficient, use_compress);

  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;

  return tt;
}
}
generate_symmetric_main(AppTwoThreeNucleusDecomposition_runner, false);
