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
//     -rounds : the number of times to run the index construction + clustering
//     -cluster-rounds : the number of times to cluster the graph per index
//       construction
//     -mu : SCAN parameter mu
//     -epsilon : SCAN parameter epsilon
#include <string>

#include "benchmarks/SCAN/IndexBased/scan.h"
#include "benchmarks/SCAN/IndexBased/similarity_measure.h"
#include "benchmarks/SCAN/IndexBased/utils.h"
#include "gbbs/gbbs.h"

namespace gbbs {
// Executes SCAN on the input graph and reports stats on the execution.
template <class Graph>
double RunScan(Graph& graph, commandLine parameters) {
  const size_t cluster_rounds{
      parameters.getOptionLongValue("-cluster-rounds", 1)};
  const uint64_t mu{parameters.getOptionLongValue("-mu", 5)};
  const float epsilon{
      static_cast<float>(parameters.getOptionDoubleValue("-epsilon", 0.6))};
  std::cout << "Scan parameters: mu = " << mu << ", epsilon = " << epsilon
            << '\n';

  timer index_construction_timer{"Index construction time"};
  const indexed_scan::Index scan_index{&graph, scan::CosineSimilarity{}};
  index_construction_timer.stop();

  timer cluster_timer{"Clustering time over " + std::to_string(cluster_rounds) +
                      " rounds"};
  for (size_t i = 0; i < cluster_rounds; i++) {
    cluster_timer.start();
    const indexed_scan::Clustering clustering{scan_index.Cluster(mu, epsilon)};
    cluster_timer.stop();
    std::cout << "Modularity: " << scan::Modularity(&graph, clustering) << '\n';
  }

  index_construction_timer.next("");
  cluster_timer.next("");
  const double running_time{index_construction_timer.total_time() +
                            cluster_timer.total_time()};
  std::cout << "Total SCAN running time: " << running_time << '\n';
  return running_time;
}

static constexpr bool kMutatesGraph{false};

}  // namespace gbbs

generate_symmetric_main(gbbs::RunScan, gbbs::kMutatesGraph);
