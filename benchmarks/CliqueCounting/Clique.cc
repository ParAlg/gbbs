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

#include "Clique.h"
#include <math.h>
#include <fstream>

//#include "kClistNodeParallel.c"

// -i 0 (simple gbbs intersect), -i 2 (simd intersect), -i 1 (set intersect)
// -space 0 = induced
// -subspace 0 = dyn, 1 = alloc, 2 = stack
// -gen
// -k clique size
// -o 0 (approx goodrich), 1 (densest using work efficient densest subgraph, exact), 2 (densest using approx densest subgraph)

// count in total, count per vert
// if count per vert, peel or no

// right now, -b = count per vert and peel; no -b means count in total
// if -b, then -f (did not do per vert relabeling for relabeled graph

template <class Graph>
double AppKCore_runner(Graph& GA, commandLine P) {
  double epsilon = P.getOptionDoubleValue("-e", 0.1); // epsilon for rank 0, 1
  long space = P.getOptionLongValue("-space", 5); // just use -space 5 forget about everything else
  long k = P.getOptionLongValue("-k", 3); // k-cliques
  long order = P.getOptionLongValue("-o", 0); //  ranking 0--4
  long recursive_level = P.getOptionLongValue("-r", 0); // -r 1 means kick off by edge, -r 2 and up is basically ineffective but same idea
  bool label = P.getOptionValue("-l"); // for -space 5, use a label intersect O(n) when setting up graph
  bool filter = P.getOptionValue("-f"); // filter only -- needed for peeeling
  bool use_base = P.getOptionValue("-b"); // this flag means count per vert + peel
  bool par_serial = P.getOptionValue("-p"); // this flag means switch to serial peeling after 200 in active
  bool sparsify = P.getOptionValue("--sparse");
  long sparsify_denom = P.getOptionLongValue("--colors", 0);
  bool approx_peel = P.getOptionValue("--approxpeel");
  double approx_eps = P.getOptionDoubleValue("--approxeps", 0.1);
  std::cout << "### Application: AppKCore" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: -k = " << k << " -e (epsilon) = " << epsilon << std::endl;
  std::cout << "### ------------------------------------" << endl;
  assert(P.getOption("-s"));

  timer t; t.start();
  size_t count = 0;
  if (sparsify) {
    auto GA_sparse = clr_sparsify_graph(GA, sparsify_denom, 7398234);
    count = Clique(GA_sparse, k, order, epsilon, space, label, filter, use_base, recursive_level, par_serial, approx_peel, approx_eps);
    std::cout << "sparse count: " << count << std::endl;
    count = count * pow(sparsify_denom,k-1);
  } else {
    count = Clique(GA, k, order, epsilon, space, label, filter, use_base, recursive_level, par_serial, approx_peel, approx_eps);
  }
  double tt = t.stop();
  std::cout << "count: " << count << std::endl;
  std::cout << "### Running Time: " << tt << std::endl;

  return tt;
}

generate_symmetric_main(AppKCore_runner, false);
