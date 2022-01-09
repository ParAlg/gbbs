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
// numactl -i all ./MaximalMatching -s -c -m clueweb_sym.bytepda
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -stats : print the #ccs, and the #vertices in the largest cc

#include "MaximalMatching.h"

#include <fstream>
#include <iomanip>
#include <iostream>

#include "gbbs/gbbs.h"

namespace gbbs {

void print_stats(commandLine& P, size_t query_cutoff, size_t max_query_length,
                 size_t total_work, double fraction_covered, double rt) {
  std::cout << "{" << std::endl;
  std::cout << "  \"test_type\": \"Yoshida matching result\"," << std::endl;
  std::cout << "  \"graph\" : \"" << P.getArgument(0) << "\"," << std::endl;
  std::cout << "  \"query_cutoff\" : " << query_cutoff << "," << std::endl;
  std::cout << "  \"time\" : " << std::setprecision(10) << rt << ","
            << std::endl;
  std::cout << "  \"max_query_length\" : " << max_query_length << ","
            << std::endl;
  std::cout << "  \"total_work\" : " << total_work << "," << std::endl;
  std::cout << "  \"fraction_covered\" : " << fraction_covered << std::endl;
  std::cout << "}" << std::endl;
}

template <class Graph>
double MaximalMatching_runner(Graph& G, commandLine P) {
  size_t query_cutoff =
      P.getOptionLongValue("-query_cutoff", std::numeric_limits<size_t>::max());
  std::cout
      << "### Application: Maximal Matching (Yoshida). Status=Experimental"
      << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -query_cutoff = " << query_cutoff << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  assert(P.getOption("-s"));  // input graph must be symmetric
  timer t;
  t.start();
  auto[max_query_length, total_work, fraction_covered] =
      MaximalMatching(G, query_cutoff);
  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;
  print_stats(P, query_cutoff, max_query_length, total_work, fraction_covered,
              tt);
  return tt;
}

}  // namespace gbbs

generate_symmetric_main(gbbs::MaximalMatching_runner, false);
