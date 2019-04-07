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
// numactl -i all ./DensestSubgraph -rounds 3 -s -m com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -fa : run the fetch-and-add implementation of k-core
//     -nb : the number of buckets to use in the bucketing implementation

#include "DensestSubgraph.h"
#include "ligra.h"

template <class vertex>
double DensestSubgraph_runner(graph<vertex>& GA, commandLine P) {
  size_t num_buckets = P.getOptionLongValue("-nb", 16);
  std::cout << "### Application: DensestSubgraph" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: n/a = " << std::endl;
  if (num_buckets != (1 << pbbslib::log2_up(num_buckets))) {
    std::cout << "Number of buckets must be a power of two."
              << "\n";
    exit(-1);
  }
  assert(P.getOption("-s"));

  timer t; t.start();
  WorkInefficientDensestSubgraph(GA, num_buckets);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

generate_main(DensestSubgraph_runner, false);
