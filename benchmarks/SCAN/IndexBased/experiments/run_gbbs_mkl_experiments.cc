// Runs SCAN to cluster a graph. Times how long the computation takes and
// measures clustering quality.
//
// Usage example:
//     bazel run //benchmarks/SCAN/IndexBased/experiments:run_gbbs_mkl_experiments -- -s <path to graph>
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//   optional:
//     -wf : indicate that the graph is weighted with floating-point numbers
//     -m : indicate that the graph should be mmap'd
#include <algorithm>
#include <vector>

#include "benchmarks/SCAN/IndexBased/scan.h"
#include "benchmarks/SCAN/IndexBased/similarity_measure.h"
#include "benchmarks/SCAN/IndexBased/utils.h"
#include "gbbs/gbbs.h"
#include "pbbslib/assert.h"

namespace gbbs {

namespace {

constexpr bool verbose{false};

struct CsvEntry {
  std::string similarity_measure{};
  int approximation_samples{0};
  std::string operation{};
  uint64_t mu{0};
  float epsilon{-1};
  double median_time{-1};
  double quality{-1};

  static void OutputHeader() {
    std::cerr
      << "similarity measure,"
      << "approximation samples,"
      << "operation,"
      << "mu,"
      << "epsilon,"
      << "median_time,"
      << "quality\n";
  }

  void Output() const {
    std::cerr << similarity_measure;
    std::cerr << ',';
    if (approximation_samples) std::cerr << approximation_samples;
    std::cerr << ',';
    std::cerr << operation;
    std::cerr << ',';
    if (mu) std::cerr << mu;
    std::cerr << ',';
    if (epsilon >= 0) std::cerr << epsilon;
    std::cerr << ',';
    if (median_time >= 0) std::cerr << median_time;
    std::cerr << ',';
    if (quality >= 0) std::cerr << quality;
    std::cerr << '\n';
  }
};

template <typename T>
std::string SimilarityToString() {
  if constexpr (std::is_same_v<T, scan::CosineSimilarity>) {
    return "cosine similarity";
  } else if constexpr (std::is_same_v<T, scan::ApproxCosineSimilarity>) {
    return "approximate cosine similarity";
  } else if constexpr (std::is_same_v<T, scan::JaccardSimilarity>) {
    return "jaccard similarity";
  } else if constexpr (std::is_same_v<T, scan::ApproxJaccardSimilarity>) {
    return "approximate jaccard similarity";
  } else if constexpr (std::is_same_v<T, scan::DenseCosineSimilarity>) {
    return "cosine similarity";
  }
}

template <typename T>
double Median(std::vector<T> v) {
  const size_t size = v.size();
  if (size == 0) {
    return 0;
  }
  std::sort(v.begin(), v.end());
  if (size % 2 == 0) {
    return (v[size / 2 - 1] + v[size / 2]) / 2;
  } else {
    return v[size / 2];
  }
}

void PrintClock() {
  if (!verbose) {
    return;
  }
  std::time_t time{std::chrono::system_clock::to_time_t(
      std::chrono::system_clock::now())};
  std::cerr << "[Clock] " << std::ctime(&time);
}

// Output index construction time as the median of several trials.
template <class Graph, class SimilarityMeasure>
indexed_scan::Index BuildIndexAndOutputTimes(
    Graph* graph,
    const SimilarityMeasure& similarity_measure,
    const size_t num_rounds) {
  std::vector<double> times;
  if (verbose) { std::cerr << "Index construction:"; }
  for (size_t i{0}; i < num_rounds - 1; i++) {
    timer timer{"Index construction"};
    const indexed_scan::Index index{graph, similarity_measure};
    times.emplace_back(timer.stop());
    if (verbose) { std::cerr << ' ' << times.back(); }
  }
  // on last trial, keep the index
  timer timer{"Index construction"};
  indexed_scan::Index index{graph, scan::CosineSimilarity{}};
  times.emplace_back(timer.stop());
  if (verbose) { std::cerr << ' ' << times.back(); }

  CsvEntry csv;
  csv.similarity_measure = SimilarityToString<SimilarityMeasure>();
  csv.operation = "index construction";
  csv.median_time = Median(times);
  csv.Output();

  return index;
}

// Executes SCAN on the input graph and reports stats on the execution.
template <class Graph>
void RunScan(Graph& graph, const commandLine& params) {
  {
    // get the graph into cache so timings are consistent
    indexed_scan::Index{&graph, scan::ApproxCosineSimilarity{1, 0}};
  }
  PrintClock();
  constexpr size_t kIndexRounds{5};
  BuildIndexAndOutputTimes(&graph, scan::DenseCosineSimilarity{}, kIndexRounds);
}

}  // namespace

}  // namespace gbbs

int main(int argc, char* argv[]) {
  gbbs::commandLine params{
    argc, argv, "-s [-wf] [--serial] <input graph file>"};
  const char* input_graph_file{params.getArgument(0)};
  const bool is_graph_symmetric{params.getOptionValue("-s")};
  ASSERT(is_graph_symmetric);
  const bool is_graph_compressed{params.getOptionValue("-c")};
  ASSERT(!is_graph_compressed);
  const bool is_graph_float_weighted{params.getOptionValue("-wf")};
  const bool should_mmap_graph{params.getOptionValue("-m")};
  const bool is_graph_binary{params.getOptionValue("-b")};
  std::cerr << "\n\nBEGIN GBBS EXPERIMENTS OUTPUT\n";
  gbbs::CsvEntry::OutputHeader();
  if (is_graph_float_weighted) {
    auto graph{gbbs::gbbs_io::read_weighted_symmetric_graph<float>(
        input_graph_file, should_mmap_graph, is_graph_binary)};
    gbbs::alloc_init(graph);
    gbbs::RunScan(graph, params);
    graph.del();
  } else {
    auto graph{gbbs::gbbs_io::read_unweighted_symmetric_graph(
        input_graph_file, should_mmap_graph, is_graph_binary)};
    gbbs::alloc_init(graph);
    gbbs::RunScan(graph, params);
    graph.del();
  }
  gbbs::alloc_finish();
}
