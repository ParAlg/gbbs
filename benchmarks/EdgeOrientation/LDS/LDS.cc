// Usage (dynamic):
// numactl -i all ./LDS -rounds 3 -s -eps 0.4 -delta 3 -i <dynamic graph> -b 1000 <static graph>
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -i : file path to the dynamic graph, with a single dynamic edge on each
//          line, in the format "<+/-> u v", where + denotes insertion, -
//          denotes deletion, and (u, v) are the vertex ids denoting the edge
//          (e.g., "+ 0 1")
//     -b : batch size (output includes time per batch, from start_size to
//          end_size)
//     -eps: epsilon
//     -delta: lambda
//     -stats : indicates whether to output comparisons to exact coreness
//              values
//     -ins-opt : indicates whether to set lambda such that
//                (2 + 3 / lambda) = 1.1
// Note: The static graph is ignored if -i is specified.

#include "LDS.h"

namespace gbbs {
template <class Graph>
double LDS_runner(Graph& G, commandLine P) {
  std::cout << "### Application: LDS" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: " <<  std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  assert(P.getOption("-s"));

  // Load dynamic graph options
  const std::string kInputFlag{"-i"};
  const char* const input_file{P.getOptionValue(kInputFlag)};
  long batch_size = P.getOptionLongValue("-b", 1);
  bool compare_exact = P.getOption("-stats");
  bool optimized_insertion = P.getOption("-ins-opt");
  bool get_size = P.getOption("-size");

  // Options not exposed to general users
  const char* const init_graph_file(P.getOptionValue("-init_graph_file"));

  // Options for the approximation algorithm
  double eps = P.getOptionDoubleValue("-eps", 3);
  double delta = P.getOptionDoubleValue("-delta", 9);

  // Load the dynamic graph
  using W = typename Graph::weight_type;
  bool use_dynamic = (input_file && input_file[0]);
  BatchDynamicEdges<W> batch_edge_list = use_dynamic ?
    read_batch_dynamic_edge_list<W>(input_file) :
    BatchDynamicEdges<W>{};
  if (use_dynamic && batch_size == 0) batch_size = batch_edge_list.edges.size();

  // Prepend the initial graph if specified, and store the offset of the
  // initial graph in the resulting dynamic edges list
  size_t offset = 0;
  if (use_dynamic && init_graph_file) {
    BatchDynamicEdges<W> init_graph_list =
      read_batch_dynamic_edge_list<W>(init_graph_file);
    offset = prepend_dynamic_edge_list(batch_edge_list, init_graph_list);
  }

  // Run LDS
  timer t; t.start();
  RunLDS(G, batch_edge_list, batch_size, compare_exact, eps, delta,
         optimized_insertion, offset, get_size);
  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;

  return tt;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::LDS_runner, false);
