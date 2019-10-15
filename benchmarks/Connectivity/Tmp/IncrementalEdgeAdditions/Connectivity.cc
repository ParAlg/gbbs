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

#include "Connectivity.h"
//#include "barabasi_albert.h"
#include "rmat.h"


template <class Graph>
double CC_runner(Graph& G, commandLine P) {
  std::cout << "### Application: Incremental Edge Connectivity (edge insertions only)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### ------------------------------------" << endl;

  /* 1. Generate batches of updates (mix of inserts and queries) */
  /* options = {
   *   a) generation_model { barabasi-albert, G(n,p), RMAT, ... },
   *   b) base graph { NONE, SOME(base_graph) },
   *   c) update/query ratio = double
   * } */
  /* models = random, 3D grid, localX, rMatX  */

  std::string model = P.getOptionValue("-model", "rmat");
  size_t num_updates = P.getOptionLongValue("-num_updates", G.n*10);
  size_t batch_size = P.getOptionLongValue("-batch_size", 1000000);
  double query_update_ratio = P.getOptionDoubleValue("-ratio", 1.0); // queries:updates 1:1

  size_t n = 1UL << pbbs::log2_up(G.n);

  pbbs::sequence<std::tuple<uintE, uintE>> updates;
  if (model == "rmat") {
    updates = rmat::generate_updates(n, num_updates, 4);
  } else {
    std::cout << "Unknown model: " << model << std::endl;
  }

  timer tt; tt.start();
  /* 2. select options for unite/find, and time each batch + total */

  size_t n_batches = (num_updates + batch_size - 1) / batch_size;

  for (size_t i=0; i<n_batches; i++) {

  }

  double t_out = tt.stop();

  /* 3. calculate throughput/time statistics to report */

  return t_out; // return total time, or total time/#batches for tim
}

generate_main(CC_runner, false);
