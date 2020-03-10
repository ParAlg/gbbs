#include "FiveCycle.h"

template <class Graph>
double Count5Cycle_runner(Graph& G, commandLine P) {
  std::cout << "ULONG_MAX: " << ULONG_MAX << std::endl;
  std::cout << "### Application: Cycle" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### ------------------------------------" << endl;
  assert(P.getOption("-s"));

  // runs the fetch-and-add based implementation if set.
  timer t; t.start();
  ulong numCycles = Count5Cycle(G);
  double tt = t.stop();
  std::cout << "### Number of Cycles: " << numCycles << std::endl;
  std::cout << "### Running Time: " << tt << std::endl;

  return tt;
}

generate_symmetric_main(Count5Cycle_runner, false);