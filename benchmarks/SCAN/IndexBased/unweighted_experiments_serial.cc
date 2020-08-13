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
#include <vector>

#include "benchmarks/SCAN/IndexBased/scan.h"
#include "benchmarks/SCAN/IndexBased/similarity_measure.h"
#include "benchmarks/SCAN/IndexBased/utils.h"
#include "gbbs/gbbs.h"
#include "pbbslib/assert.h"

namespace gbbs {

namespace {

template <typename T>
double Median(std::vector<T> v) {
  const size_t size = v.size();
  if (size == 0) {
    return 0;
  }
  sort(v.begin(), v.end());
  if (size % 2 == 0) {
    return (v[size / 2 - 1] + v[size / 2]) / 2;
  } else {
    return v[size / 2];
  }
}

std::string ParametersToString(const size_t mu, const float epsilon) {
  return std::to_string(mu) + "-" + std::to_string(epsilon);
}

}  // namespace

// Executes SCAN on the input graph and reports stats on the execution.
template <class Graph>
double RunScan(Graph& graph, commandLine parameters) {
  const size_t max_degree{pbbslib::reduce_max(pbbslib::make_sequence<size_t>(
    graph.n,
    [&](const size_t i) {
      return graph.get_vertex(i).getOutDegree();
    }))};
  {
    // get the graph into cache so timings are consistent
    indexed_scan::Index{&graph, scan::ApproxCosineSimilarity{1, 0}};
  }

  constexpr size_t index_rounds{5};
  {
    std::cerr << "**********************\n";
    std::cerr << "** CosineSimilarity **\n";
    std::cerr << "**********************\n";

    std::vector<double> exact_index_times;
    std::cerr << "Index construction:";
    for (size_t i{0}; i < index_rounds - 1; i++) {
      timer exact_index_timer{"Index construction"};
      const indexed_scan::Index exact_index{&graph, scan::CosineSimilarity{}};
      exact_index_times.emplace_back(exact_index_timer.stop());
      std::cerr << ' ' << exact_index_times.back();
    }
    // on last trial, keep the index
    timer exact_index_timer{"Index construction"};
    const indexed_scan::Index exact_index{&graph, scan::CosineSimilarity{}};
    exact_index_times.emplace_back(exact_index_timer.stop());
    std::cerr << ' ' << exact_index_times.back();
    std::cerr << "\n** Index construction median: " << Median(exact_index_times) << "\n\n";

    // timing trials for querying
    constexpr size_t cluster_rounds{5};
    {
      for (size_t mu{2}; mu < max_degree; mu *= 2) {
        std::vector<double> query_times;
        constexpr double kEpsilon{0.6};
        std::cerr << "query " << ParametersToString(mu, kEpsilon) << ":";
        for (size_t i{0}; i < cluster_rounds; i++) {
          timer query_timer{};
          // unused variable `clusters`, but keep it anyway so that it doesn't get
          // destructed within the timer
          const scan::Clustering clusters{exact_index.Cluster(mu, kEpsilon)};
          query_times.emplace_back(query_timer.stop());
          std::cerr << ' ' << query_times.back();
        }
        std::cerr << "\n** Cluster median " << ParametersToString(mu, kEpsilon) << ": " << Median(query_times) << "\n";
      }
    }
    {
      for (double epsilon{0.01}; epsilon < 0.995; epsilon += 0.01) {
        std::vector<double> query_times;
        constexpr double kMu{5};
        std::cerr << "query " << ParametersToString(kMu, epsilon) << ":";
        for (size_t i{0}; i < cluster_rounds; i++) {
          timer query_timer{"query"};
          const scan::Clustering clusters{exact_index.Cluster(kMu, epsilon)};
          query_times.emplace_back(query_timer.stop());
          std::cerr << ' ' << query_times.back();
        }
        std::cerr << "\n** Cluster median " << ParametersToString(kMu, epsilon) << ": " << Median(query_times) << "\n";
      }
    }

    std::cerr << "****************************\n";
    std::cerr << "** ApproxCosineSimilarity **\n";
    std::cerr << "****************************\n";

    for (uint32_t num_samples{64}; 2 * num_samples < max_degree; num_samples *= 2) {
      std::cerr << "Samples=" << num_samples << '\n';
      std::cerr << "----------\n";
      std::vector<double> approx_index_times;
      std::cerr << "Index construction:";
      for (size_t i{0}; i < index_rounds; i++) {
        timer approx_index_timer{"Index construction"};
        const indexed_scan::Index approx_index{&graph, scan::ApproxCosineSimilarity{num_samples, i}};
        approx_index_times.emplace_back(approx_index_timer.stop());
        std::cerr << ' ' << approx_index_times.back();
      }
      std::cerr << "\n** Index construction median: " << Median(approx_index_times) << "\n\n";
    }
  }

  {
    std::cerr << "*****************************\n";
    std::cerr << "** ApproxJaccardSimilarity **\n";
    std::cerr << "*****************************\n";
    for (uint32_t num_samples{64}; 2 * num_samples < max_degree; num_samples *= 2) {
      std::cerr << "Samples=" << num_samples << '\n';
      std::cerr << "----------\n";
      std::vector<double> approx_index_times;
      std::cerr << "Index construction:";
      for (size_t i{0}; i < index_rounds; i++) {
        timer approx_index_timer{"Index construction"};
        const indexed_scan::Index approx_index{&graph, scan::ApproxJaccardSimilarity{num_samples, i}};
        approx_index_times.emplace_back(approx_index_timer.stop());
        std::cerr << ' ' << approx_index_times.back();
      }
      std::cerr << "\n** Index construction median: " << Median(approx_index_times) << "\n\n";
    }
  }

  return 0.0;
}

static constexpr bool kMutatesGraph{false};

}  // namespace gbbs

generate_symmetric_once_main(gbbs::RunScan, gbbs::kMutatesGraph);
