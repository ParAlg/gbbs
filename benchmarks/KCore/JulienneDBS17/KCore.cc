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

namespace gbbs {

template <class Graph>
double KCore_runner(Graph& G, commandLine P) {
  size_t num_buckets = P.getOptionLongValue("-nb", 16);
  bool fa = P.getOption("-fa");
  std::cout << "### Application: KCore" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -nb (num_buckets) = " << num_buckets << " -fa (use fetch_and_add) = " << fa << std::endl;
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
  long num_dynamic_edges = P.getOptionLongValue("-num_dynamic_edges", 0);
  long batch_size = P.getOptionLongValue("-b", 1);
  bool use_stats = P.getOptionValue("-stats");
  long start_size = P.getOptionLongValue("-start_size", 0);
  long end_size = P.getOptionLongValue("-end_size", 0);

  bool use_dynamic = (input_file && input_file[0]);

  if (!use_dynamic) {
    timer t; t.start();
    auto cores = (fa) ? KCore_FA(G, num_buckets) : KCore(G, num_buckets);
    double tt = t.stop();
    uintE max_core = parlay::reduce(cores, parlay::maxm<uintE>());
    std::cout << "### Coreness Exact: " << max_core << std::endl;
    std::cout << "### Running Time: " << tt << std::endl;
    return tt;
  }
  
  using W = typename Graph::weight_type;
  BatchDynamicEdges<W> batch_edge_list = use_dynamic ?
    read_batch_dynamic_edge_list<W>(input_file) : BatchDynamicEdges<W>();
  
  if (num_dynamic_edges != 0) {
    symmetric_graph<symmetric_vertex, W> dynamic_graph = 
      dynamic_edge_list_to_symmetric_graph(batch_edge_list, use_dynamic ? num_dynamic_edges : 0);

    // runs the fetch-and-add based implementation if set.
    timer t; t.start();
    auto cores = (fa) ? KCore_FA(dynamic_graph, num_buckets) : KCore(dynamic_graph, num_buckets);
    double tt = t.stop();

    std::cout << "### Batch Running Time: " << tt << std::endl;
    std::cout << "### Batch Num: " << num_dynamic_edges << std::endl;
    if (use_stats) {
      uintE max_core = parlay::reduce(cores, parlay::maxm<uintE>());
      std::cout << "### Coreness Exact: " << max_core << std::endl;
    }
    return tt;
  }

  timer t1; t1.start();
  auto batch = batch_edge_list.edges;
  if (end_size == 0) end_size = batch.size();
  for (size_t i = start_size; i < end_size; i += batch_size) {
    num_dynamic_edges = std::min((size_t)end_size, i + batch_size);
    auto dynamic_graph = dynamic_edge_list_to_symmetric_graph(batch_edge_list, num_dynamic_edges);
    timer t; t.start();
    auto cores = (fa) ? KCore_FA(dynamic_graph, num_buckets) : KCore(dynamic_graph, num_buckets);
    double tt = t.stop();
    std::cout << "### Batch Running Time: " << tt << std::endl;
    std::cout << "### Batch Num: " << num_dynamic_edges << std::endl;
    if (use_stats) {
      uintE max_core = parlay::reduce(cores, parlay::maxm<uintE>());
      std::cout << "### Coreness Exact: " << max_core << std::endl;
    }
  }
  double tt1 = t1.stop();
  
  return tt1;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::KCore_runner, false);
