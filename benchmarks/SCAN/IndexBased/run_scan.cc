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
  timer index_construction_timer{"Index construction time"};
  const indexed_scan::Index scan_index{&graph, scan::CosineSimilarity{}};
  index_construction_timer.reportTotal("");

  constexpr uint32_t kMu{3};
  constexpr bool kDeterministic{true};
  pbbs::sequence<float> epsilons{
      99, [&](const size_t i) { return (i + 1) * .01;}};
  pbbs::sequence<scan::Clustering> clusterings{
    scan_index.Cluster(kMu, epsilons, kDeterministic)};
  bool good{true};
  for (size_t i = 0; i < 99; i++) {
    scan::Clustering clustering{
      scan_index.Cluster(kMu, epsilons[i], kDeterministic)};
    bool check{clustering == clusterings[i]};
    // std::cerr << epsilons[i] << ": " << check << '\n';
    if (!check) {
      good = false;
    }
  }
  if (!good) {
    std::cerr << "######## FAILURE" << std::endl;
  } else {
    std::cerr << "!! Good!\n";
  }

  for (int i = 0; i < 3; i++) {
    timer cluster_sweep_timer{"Clustering sweep time"};
    scan_index.Cluster(kMu, epsilons);
    cluster_sweep_timer.reportTotal("");

    timer cluster_timer{"Clustering individual time"};
    for (auto e : epsilons) {
      scan_index.Cluster(kMu, e);
    }
    cluster_timer.reportTotal("");
  }

  return 0.0;
}

static constexpr bool kMutatesGraph{false};

}  // namespace gbbs

generate_symmetric_once_main(gbbs::RunScan, gbbs::kMutatesGraph);
