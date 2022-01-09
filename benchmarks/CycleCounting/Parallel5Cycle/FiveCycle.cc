#include "FiveCycle.h"

namespace gbbs {
template <class Graph>
std::tuple<ulong, double> Count5Cycle_runner(Graph& G, long order_type,
                                             bool experiment, bool escape,
                                             bool no_schedule, bool serial) {
  std::cout << "### Direct Type: ";
  if (order_type == 0) {
    std::cout << "Goodrich-Pszona" << std::endl;
  } else if (order_type == 1) {
    std::cout << "Barenboim-Elkin" << std::endl;
  } else if (order_type == 2) {
    std::cout << "Degree ordering" << std::endl;
  }

  // runs the fetch-and-add based implementation if set.
  std::cout << "### Threads: " << num_workers() << std::endl;
  ulong numCycles;
  timer t;
  t.start();
  if (experiment) {
    std::cout << "### Experiment (parallel)" << std::endl;
    numCycles = Count5Cycle_experiment(G, order_type);
  } else if (escape) {
    if (serial) {
      std::cout << "### ESCAPE (Pure Serial)" << std::endl;
      numCycles = Count5Cycle_ESCAPE(G, order_type);
    } else {
      std::cout << "### ESCAPE (Parallel)" << std::endl;
      numCycles = Count5Cycle_ESCAPE_par(G, order_type);
    }
  } else {
    if (serial) {
      std::cout << "### Kowalik (Pure Serial)" << std::endl;
      numCycles = Count5Cycle_serial(G, order_type);
    } else if (no_schedule) {
      std::cout << "### Kowalik (Parallel No Scheduling)" << std::endl;
      numCycles = Count5Cycle_no_scheduling(G, order_type);
    } else {
      std::cout << "### Kowalik (Parallel)" << std::endl;
      numCycles = Count5Cycle(G, order_type);
    }
  }
  double tt = t.stop();
  std::cout << "### Number of Cycles: " << numCycles << std::endl;
  std::cout << "### Running Time: " << tt << std::endl;

  return std::make_tuple(numCycles, tt);
}

template <class Graph>
double Count5Cycle_runner(Graph& G, commandLine P) {
  std::cout << "ULONG_MAX: " << ULONG_MAX << std::endl;
  std::cout << "### Application: Cycle" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  assert(P.getOption("-s"));  // make sure input graph is symmetric

  bool sparsify = P.getOptionValue(
      "--sparse");  // if set, use colorful sparsification for approx counting
  long sparsify_denom = P.getOptionLongValue(
      "--colors", 0);  // number of colors for colorful sparsification
  long esparsify_denom = P.getOptionLongValue(
      "--edenom", 0);  // probability of keeping an edge for edge sparsification

  bool serial = P.getOptionValue("--serial");
  bool no_schedule = P.getOptionValue("--no-schedule");
  bool experiment = P.getOptionValue("--exp");
  bool escape = P.getOptionValue("--escape");
  long order_type = P.getOptionLongValue("-o", 0);

  if (sparsify) {
    // Colorful sparsify graph, with random seed
    if (sparsify_denom != 0) {
      auto G_sparse = clr_sparsify_graph(G, sparsify_denom, 7398234);
      auto ret = Count5Cycle_runner(G_sparse, order_type, experiment, escape,
                                    no_schedule, serial);
      auto count = std::get<0>(ret);
      count = count * pow(sparsify_denom, 4);  // TODO: fix
      std::cout << "un-sparse count: " << count << std::endl;
      return std::get<1>(ret);
    } else if (esparsify_denom != 0) {
      // edge sparsify
      auto G_sparse = edge_sparsify_graph(G, esparsify_denom, 7398234);
      auto ret = Count5Cycle_runner(G_sparse, order_type, experiment, escape,
                                    no_schedule, serial);
      auto count = std::get<0>(ret);
      count = count * pow(esparsify_denom, 5);  // TODO: fix
      std::cout << "un-sparse count: " << count << std::endl;
      return std::get<1>(ret);
    }
  }

  return std::get<1>(Count5Cycle_runner(G, order_type, experiment, escape,
                                        no_schedule, serial));
}
}  // namespace gbbs

generate_symmetric_main(gbbs::Count5Cycle_runner, false);