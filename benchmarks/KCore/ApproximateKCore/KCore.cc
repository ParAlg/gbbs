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

#include "KCore.h"
#include "gbbs/dynamic_graph_io.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

namespace gbbs {
template <class Graph>
double KCore_runner(Graph& G, commandLine P) {
  size_t num_buckets = P.getOptionLongValue("-nb", 16);
  double eps = P.getOptionDoubleValue("-eps", 0.2);
  double delta = P.getOptionDoubleValue("-delta", 0.1);
  bool use_pow = P.getOptionValue("-use_pow");
  std::cout << "### Application: KCore" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -nb (num_buckets) = " << num_buckets << " epsilon = " << eps << " use_pow = " << use_pow << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  if (num_buckets != static_cast<size_t>((1 << pbbslib::log2_up(num_buckets)))) {
    std::cout << "Number of buckets must be a power of two."
              << "\n";
    exit(-1);
  }
  assert(P.getOption("-s"));

  // Option to use dynamic graph instead of static
  const std::string kInputFlag{"-i"};
  const char* const input_file{P.getOptionValue(kInputFlag)};
  long num_dynamic_edges = P.getOptionLongValue("-num_dynamic_edges", 1);

  bool use_dynamic = (input_file && input_file[0]);
  
  using W = typename Graph::weight_type;
  BatchDynamicEdges<W> batch_edge_list = use_dynamic ?
    read_batch_dynamic_edge_list<W>(input_file) : BatchDynamicEdges<W>();
  if (use_dynamic && num_dynamic_edges == 0) num_dynamic_edges = batch_edge_list.edges.size();
  symmetric_graph<symmetric_vertex, W> dynamic_graph = 
    dynamic_edge_list_to_symmetric_graph(batch_edge_list, use_dynamic ? num_dynamic_edges : 0);

  // runs the fetch-and-add based implementation if set.
  timer t; t.start();
  auto cores = use_dynamic ? approximate_kcore::KCore(dynamic_graph, num_buckets, eps, delta, use_pow) : 
    approximate_kcore::KCore(G, num_buckets, eps, delta, use_pow);
  double tt = t.stop();

  if (use_dynamic) {
    std::cout << "### Batch Running Time: " << tt << std::endl;
    std::cout << "### Batch Num: " << num_dynamic_edges << std::endl;
  }

  uintE max_core = parlay::reduce(cores, parlay::maxm<uintE>());
  std::cout << "### Coreness Estimate: " << max_core << std::endl;

  double mult_appx = (2 + 2*eps);
  if (P.getOptionValue("-stats")) {
    auto true_cores = use_dynamic ? KCore(dynamic_graph, num_buckets) : KCore(G, num_buckets);
    auto n = use_dynamic ? dynamic_graph.n : G.n;

    uintE max_true_core = parlay::reduce(true_cores, parlay::maxm<uintE>());
    std::cout << "### Coreness Exact: " << max_true_core << std::endl;

    size_t bad = 0;
    double total_error = 0.0;
    double max_error = 0.0;
    double min_error = std::numeric_limits<double>::max();
    double denominator = 0;
    for (size_t i=0; i<n; i++) {
      double true_core = true_cores[i];
      double appx_core = cores[i];
      if (appx_core > (mult_appx*true_core)) {
        //std::cout << "overappx, true_core = " << true_core << " appx_core = " << appx_core << std::endl;
        bad++;
      }
      if (appx_core < (true_core / mult_appx)) {
        //std::cout << "underappx, true_core = " << true_core << " appx_core = " << appx_core << std::endl;
        bad++;
      }
      if (true_core != 0 && appx_core != 0) {
        auto this_error = std::max(true_core, appx_core) / std::min(true_core, appx_core);
        total_error += this_error;
        max_error = std::max(max_error, this_error);
        min_error = std::min(min_error, this_error);
        denominator++;
      }
    }
    min_error = min_error == std::numeric_limits<double>::max() ? 0 : min_error;
    auto avg_error = (denominator == 0) ? 0 : (total_error / denominator);
    std::cout << "### Num Bad: " << bad << std::endl;
    std::cout << "### Per Vertex Average Coreness Error: " << avg_error << std::endl;
    std::cout << "### Per Vertex Min Coreness Error: " << min_error << std::endl;
    std::cout << "### Per Vertex Max Coreness Error: " << max_error << std::endl;
  }
  
  if (!use_dynamic) std::cout << "### Running Time: " << tt << std::endl;

  return tt;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::KCore_runner, false);
