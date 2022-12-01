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

#include "NucleusDecomposition.h"
#include <math.h>
#include <fstream>
namespace gbbs {
template <class Graph>
double AppNucleusDecompositionVerification_runner(Graph& GA, commandLine P) {
  long r = P.getOptionLongValue("-r", 3); // k as in k-cliques
  long ss = P.getOptionLongValue("-ss", 4); // k as in k-cliques

  long table_type = P.getOptionLongValue("-tt", 3); // 1 = 1 lvl, 2 = 2 lvls, 3 = multi
  long num_levels = P.getOptionLongValue("-nl", 2); // only for multi, # levels
  bool relabel = P.getOptionValue("-relabel"); // for true, relabel graph
  bool contiguous_space = P.getOptionValue("-contig"); // for true, contiguous space

  long table_type_verify = P.getOptionLongValue("-ttv", 3); // 1 = 1 lvl, 2 = 2 lvls, 3 = multi
  long num_levels_verify = P.getOptionLongValue("-nlv", 2); // only for multi, # levels
  bool relabel_verify = P.getOptionValue("-relabelv"); // for true, relabel graph
  bool contiguous_space_verify = P.getOptionValue("-contigv"); // for true, contiguous space


  std::cout << "### Application: Nucleus Decomposition Verification" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  assert(P.getOption("-s"));

  timer t; t.start();
  auto count = NucleusDecomposition(GA, r, ss, table_type, num_levels, relabel, contiguous_space);
  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;

  timer t2; t2.start();
  auto count2 = NucleusDecomposition(GA, r, ss, table_type_verify,
    num_levels_verify, relabel_verify, contiguous_space_verify);
  double tt2 = t2.stop();
  std::cout << "### Running Time: " << tt2 << std::endl;

  return tt;
}
}
generate_symmetric_main(AppNucleusDecompositionVerification_runner, false);
