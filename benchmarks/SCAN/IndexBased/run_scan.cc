// Runs SCAN to cluster a graph and times how long the computation takes.
//
// Usage example:
//     bazel run //benchmarks/SCAN/IndexBased:SCAN_main -- -s <path to graph>
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
#include "benchmarks/SCAN/IndexBased/scan.h"

#include "ligra/ligra.h"

// Executes SCAN on the input graph and reports stats on the execution.
template <class Graph>
double RunScan(Graph& graph, commandLine parameters) {
  timer timer{};
  timer.start();
  const indexed_scan::Index scan_index{&graph};

  constexpr uint64_t kMu{5};
  constexpr float kEpsilon{0.6};
  const indexed_scan::Clustering clustering{scan_index.Cluster(kMu, kEpsilon)};
  const double running_time{timer.stop()};
  std::cout << "Running Time: " << running_time << std::endl;
  return running_time;
}

static constexpr bool kMutatesGraph{false};
generate_symmetric_main(RunScan, kMutatesGraph);
