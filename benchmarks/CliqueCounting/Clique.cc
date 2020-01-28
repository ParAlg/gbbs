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
  double epsilon = P.getOptionDoubleValue("-e", 0.1);
  long space = P.getOptionLongValue("-space", 2);
  long k = P.getOptionLongValue("-k", 3);
  long order = P.getOptionLongValue("-o", 0);
  bool label = P.getOptionValue("-l");
  bool filter = P.getOptionValue("-f");
  bool use_base = P.getOptionValue("-b");
  std::cout << "### Application: AppKCore" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: -k = " << k << " -e (epsilon) = " << epsilon << std::endl;
  std::cout << "### ------------------------------------" << endl;
  assert(P.getOption("-s"));

  /*timer tclist; tclist.start();
  std::string file = P.getOptionValue("-file", "");
  auto countcount = kClist(k, file.c_str());
  double ttclist = tclist.stop();
  std::cout << "count: " << countcount << std::endl;
  std::cout << "### Running Time: " << ttclist << std::endl;
  return ttclist;*/

  /*std::string rankfile = P.getOptionValue("-rankfile", "");
  uintE* r = nullptr;
  if (rankfile != "") {
    std::ifstream infile(rankfile);
    r = pbbslib::new_array_no_init<uintE>(GA.n);
    int a;
    size_t idx = 0;
    while (infile >> a) {
      r[idx++] = a;
    }
  }*/


  timer t; t.start();
  uintE* per_vert = use_base ? (uintE*) malloc(sizeof(uintE)*GA.n) : nullptr;
  size_t count = KClique(GA, k, order, epsilon, space, label, filter, use_base, per_vert);
  double tt = t.stop();
  std::cout << "count: " << count << std::endl;
  std::cout << "### Running Time: " << tt << std::endl;

  //TODO peeling
  if (!use_base) return tt;
  assert(filter);
  size_t allcount=0;
  for (size_t j=0; j < GA.n; j++) {
    allcount += per_vert[j];
  }
  assert(allcount == count);
  timer t2; t2.start();
  sequence<uintE> cores = Peel(GA, k-1, per_vert, label);
  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;
  free(per_vert);

  return tt;
}

generate_symmetric_main(AppKCore_runner, false);
