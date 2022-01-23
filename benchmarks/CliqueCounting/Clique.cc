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

namespace gbbs {

long strToDirectType(std::string order_str) {
  if (order_str == "GOODRICHPSZONA")
    return 0;
  else if (order_str == "BARENBOIMELKIN")
    return 1;
  else if (order_str == "KCORE")
    return 2;
  else if (order_str == "DEGREE")
    return 3;
  else if (order_str == "ORIGINAL")
    return 4;
  ABORT("Unexpected directed type (str): " << order_str);
}

long strToParallelType(std::string rec_str) {
  if (rec_str == "VERT")
    return 0;
  else if (rec_str == "EDGE")
    return 1;
  ABORT("Unexpected parallel type (str): " << rec_str);
}

template <class Graph>
double AppKCore_runner(Graph& GA, commandLine P) {
  // These are user-facing options
  long k = P.getOptionLongValue("-k", 3);  // k as in k-cliques
  double epsilon = P.getOptionDoubleValue(
      "-e", 0.1);  // epsilon, for Goodrich-Pszona or Barenboim-Elkin directing

  // Options to direct graph: GOODRICHPSZONA, BARENBOIMELKIN, KCORE, DEGREE, or
  // ORIGINAL
  auto order_str = P.getOptionValue("--directType", "");
  auto recursive_level_str = P.getOptionValue(
      "--parallelType", "");  // parallelism per vert (VERT) or per edge (EDGE)

  // If set, do not use linear space to kick off the first level of recursion
  // (uses O(alpha^2) space)
  bool label_str = P.getOptionValue("--saveSpace");
  bool use_base_str = P.getOptionValue("--peel");  // run vertex peeling

  bool sparsify = P.getOptionValue(
      "--sparse");  // if set, use colorful sparsification for approx counting
  long sparsify_denom = P.getOptionLongValue(
      "--colors", 0);  // number of colors for colorful sparsification

  bool approx_peel = P.getOptionValue(
      "--approxpeel");  // if set, use approximate vertex peeling
  double approx_eps = P.getOptionDoubleValue(
      "--approxeps", 0.1);  // epsilon for approximate vertex peeling

  // These are debugging options (set to correct defaults for client)
  long order = P.getOptionLongValue(
      "-o", 0);  // directing type; values are 0 - 4, reflecting the above order
  long recursive_level =
      P.getOptionLongValue("-r", 0);        // parallelism over recursive levels
  bool use_base = P.getOptionValue("-ub");  // same as use_base_str
  bool filter = P.getOptionValue(
      "-f");  // filter graph only -- required for vertex peeling

  long space = P.getOptionLongValue("-space", 5);

  std::cout << "### Application: Clique Counting" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: -k = " << k << " -e (epsilon) = " << epsilon
            << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  assert(P.getOption("-s"));

  // Set user-facing options to match backend
  if (!order_str.empty()) order = strToDirectType(order_str);
  if (!recursive_level_str.empty())
    recursive_level = strToParallelType(recursive_level_str);
  bool label = !label_str;
  use_base = use_base || use_base_str;
  if (use_base) filter = true;

  timer t;
  t.start();

  size_t count = 0;
  if (sparsify) {
    // Sparsify graph, with random seed
    clr_sparsify_graph(GA, sparsify_denom, 7398234);

    // k-clique counting
    count = Clique(GA, k, order, epsilon, space, label, filter, use_base,
                   recursive_level, approx_peel, approx_eps);
    std::cout << "sparse count: " << count << std::endl;
    count = count * pow(sparsify_denom, k - 1);
  } else {
    // k-clique counting
    count = Clique(GA, k, order, epsilon, space, label, filter, use_base,
                   recursive_level, approx_peel, approx_eps);
  }

  double tt = t.stop();

  std::cout << "count: " << count << std::endl;
  std::cout << "### Running Time: " << tt << std::endl;

  return tt;
}

}  // namespace gbbs

generate_symmetric_main(gbbs::AppKCore_runner, false);
