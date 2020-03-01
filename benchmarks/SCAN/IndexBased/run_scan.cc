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
  timer index_construction_timer{};
  index_construction_timer.start();
  const indexed_scan::Index scan_index{&graph};
  const double index_construction_time{index_construction_timer.stop()};

  constexpr uint64_t kMu{5};
  constexpr float kEpsilon{0.6};
  timer cluster_timer{};
  cluster_timer.start();
  const indexed_scan::Clustering clustering{scan_index.Cluster(kMu, kEpsilon)};
  const double cluster_time{cluster_timer.stop()};

  std::cout << "Index construction time: " << index_construction_time << '\n';;
  std::cout << "Clustering time: " << cluster_time << '\n';;
  const double running_time{index_construction_time + cluster_time};
  std::cout << "Total running time: " << running_time << '\n';
  return running_time;
}

static constexpr bool kMutatesGraph{false};
generate_symmetric_main(RunScan, kMutatesGraph);
