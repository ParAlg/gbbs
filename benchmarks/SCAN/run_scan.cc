// Runs SCAN to cluster a graph and times how long the computation takes.
//
// Usage example:
// numactl -i all ./scan -rounds 3 -s -m twitter_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc

#include "benchmarks/SCAN/scan.h"

#include "ligra/ligra.h"
#include "utils/assert.h"

// Execute SCAN on the input graph and report stats on the execution.
template <class Graph>
double RunScan(Graph& graph, commandLine parameters) {
  ASSERT(parameters.getOption("-s"), "Input graph must be symmetric");

  timer timer{};
  timer.start();
  const scan::ScanIndex scan_index{&graph};

  constexpr uint64_t kMu{5};
  constexpr float kEpsilon{0.6};
  const scan::Clustering clustering{scan_index.Cluster(kMu, kEpsilon)};
  const double running_time{timer.stop()};
  std::cout << "Running Time: " << running_time << std::endl;
  return running_time;
}

static const bool kMutates{false};
generate_symmetric_main(RunScan, kMutates);
