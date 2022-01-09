// Copyright (c) 2019 Laxman Dhulipala, Guy Blelloch, and Julian Shun
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
// numactl -i all ./KTruss -rounds 3 -s -m com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -no_buckets : run an implementation that does not use bucketing
//                   (using bucketing is default)
//     -nb : the number of buckets to use in the bucketing implementation

#include "KTruss.h"

namespace gbbs {

template <class Graph>
double KTruss_runner(Graph& G, commandLine P) {
  size_t num_buckets = P.getOptionLongValue("-nb", 16);
  bool no_buckets = P.getOption("-no_buckets");
  std::cout << "### Application: KTruss" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -nb (num_buckets) = " << num_buckets
            << " -no_buckets = " << no_buckets << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  if (num_buckets != static_cast<size_t>((1 << parlay::log2_up(num_buckets)))) {
    std::cout << "Number of buckets must be a power of two."
              << "\n";
    exit(-1);
  }
  assert(P.getOption("-s"));

  // runs the fetch-and-add based implementation if set.
  timer t;
  t.start();
  // auto trusses = (!no_buckets) ? KTruss(G, num_buckets) :
  // KTruss_no_bucket(G);
  KTruss_ht(G, num_buckets);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;

  exit(0);  // Note that running k-truss mutates the input graph to enable the
            // packing optimization.
  return tt;
}

}  // namespace gbbs

generate_symmetric_main(gbbs::KTruss_runner, false);
