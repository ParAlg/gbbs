// Usage (dynamic):

#include "StreamingTriangle.h"

namespace gbbs {
template <class Graph>
double LDS_runner(Graph& G, commandLine P) {
  std::cout << "### Application: StreamingTriangle" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: " <<  std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  assert(P.getOption("-s"));

  // Load dynamic graph options
  const std::string kInputFlag{"-i"};
  size_t multiplier = P.getOptionLongValue("-mult", 5);
  size_t trials = P.getOptionLongValue("-trials", 10);
  size_t triangles = P.getOptionLongValue("-triangles", 1);
  double eps = P.getOptionDoubleValue("-eps", 0.2);
  const char* const input_file{P.getOptionValue(kInputFlag)};

  // Load the dynamic graph
  using W = typename Graph::weight_type;
  bool use_dynamic = (input_file && input_file[0]);
  BatchDynamicEdges<W> batch_edge_list = use_dynamic ?
    read_batch_dynamic_edge_list<W>(input_file) : BatchDynamicEdges<W>{};

  // Run StreamingTriangle
  ApproximateTriangles(batch_edge_list, multiplier, trials, P, triangles, eps);

  return 0;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::LDS_runner, false);
