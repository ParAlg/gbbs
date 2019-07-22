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
// numactl -i all ./Spanner -rounds 3 -s -m twitter_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc

#include <math.h>

#include "Spanner.h"
#include "ligra.h"

// Beta should be set to log n/2k. See Corollary 3.1 and Lemma 3.2 in MPVX'15.
template <class G>
double Spanner_runner(G& GA, commandLine P) {
  size_t n = GA.n;
  size_t k = P.getOptionLongValue("-k", 4);
  double beta = log(n)/(2*k);
  std::cout << "### Application: Spanner (O(k)-spanner from MPXV)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: -k = " << k << " => \\beta = \\log n/2k = " << beta << std::endl;
  std::cout << "### ------------------------------------" << endl;

  assert(P.getOption("-s"));
  timer t;
  t.start();
  auto spanner = spanner::Spanner(GA, beta);
  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;
  std::cout << "### ------------------------------------" << endl;
  return tt;
}

generate_main(Spanner_runner, false);
