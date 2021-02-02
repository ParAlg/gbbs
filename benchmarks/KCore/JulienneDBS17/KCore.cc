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
    auto batch = batch_edge_list.edges;
    // Do an in-place sort for inserts and deletes in the range i to num_dynamic_edges
    auto compare_tup = [&] (const DynamicEdge<W>& l, const DynamicEdge<W>& r) { return l.insert > r.insert; };
    parlay::sort_inplace(parlay::make_slice(batch.data(), batch.data() + num_dynamic_edges), compare_tup);
    // Find the first deletion
    size_t j = 0;
    for (j = 0; j < num_dynamic_edges; j++) {
      if (!batch[j].insert) break;
    }

    timer t; t.start();
    // First do insertions
    symmetric_graph<symmetric_vertex, W> dynamic_graph_insert = dynamic_edge_list_to_symmetric_graph(batch_edge_list, j);
    auto cores_insert = (fa) ? KCore_FA(dynamic_graph_insert, num_buckets) : KCore(dynamic_graph_insert, num_buckets);

    // Then do deletions
    symmetric_graph<symmetric_vertex, W> dynamic_graph =  dynamic_edge_list_to_symmetric_graph(batch_edge_list, num_dynamic_edges);

    // runs the fetch-and-add based implementation if set.
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

  auto batch = batch_edge_list.edges;
  if (end_size == 0 || end_size > batch.size()) end_size = batch.size();
  auto end_seq = parlay::sequence<uintE>(3 + 2 * ((end_size - start_size) / batch_size), (uintE)0);
  end_seq[0] = start_size;
  for (size_t i = start_size; i < end_size; i += batch_size) {
    num_dynamic_edges = std::min((size_t)end_size, i + batch_size);

    // Do an in-place sort for inserts and deletes in the range i to num_dynamic_edges
    auto compare_tup = [&] (const DynamicEdge<W>& l, const DynamicEdge<W>& r) { return l.insert > r.insert; };
    parlay::sort_inplace(parlay::make_slice(batch.data() + i, batch.data() + num_dynamic_edges), compare_tup);

    // Find the first insert
    size_t j = i;
    for (j = i; j < num_dynamic_edges; j++) {
      if (!batch[j].insert) break;
    }
    end_seq[1 + 2 * ((i - start_size) / batch_size)] = j;
    end_seq[1 + 2 * ((i - start_size) / batch_size) + 1] = num_dynamic_edges;
  }

  timer t1; t1.start();
  for (size_t i = 0; i < end_seq.size() - 1; i++) {
    size_t start = end_seq[i];
    size_t end = end_seq[i + 1];
    if (i != 0 && (start == 0 || end == 0)) break;
  
    num_dynamic_edges = end;
  
    auto dynamic_graph = dynamic_edge_list_to_symmetric_graph(batch_edge_list, num_dynamic_edges);
    timer t; t.start();
    auto cores = (fa) ? KCore_FA(dynamic_graph, num_buckets) : KCore(dynamic_graph, num_buckets);
    double tt = t.stop();
    std::cout << "### Batch Running Time: " << tt << std::endl;
    if (num_dynamic_edges % batch_size == 0 || num_dynamic_edges == batch.size())
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
