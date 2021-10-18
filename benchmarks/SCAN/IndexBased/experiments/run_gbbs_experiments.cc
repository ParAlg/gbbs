// Runs SCAN to cluster a graph. Times how long the computation takes and
// measures clustering quality.
//
// Usage example:
//     bazel run //benchmarks/SCAN/IndexBased/experiments:run_gbbs_experiments -- -s <path to graph>
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//   optional:
//     -wf : indicate that the graph is weighted with floating-point numbers
//     --serial: indicate that we're only gathering serial timings and not
//       measuring clustering quality
//     -m : indicate that the graph should be mmap'd
#include <algorithm>
#include <string>
#include <vector>

#define PY_SSIZE_T_CLEAN
#include <python3.6m/Python.h>

#include "benchmarks/SCAN/IndexBased/scan.h"
#include "benchmarks/SCAN/IndexBased/similarity_measure.h"
#include "benchmarks/SCAN/IndexBased/utils.h"
#include "gbbs/gbbs.h"
#include "gbbs/helpers/assert.h"

namespace gbbs {

namespace {

constexpr bool verbose{false};

// Wrapper class for computing the adjusted Rand index between two clusterings.
class ARIQuerier {
 public:
  ARIQuerier();
  ~ARIQuerier();

  double
  AdjustedRandIndex(const scan::Clustering& a, const scan::Clustering& b) const;

 private:
  PyObject* py_module_name_{nullptr};
  PyObject* py_module_{nullptr};
  PyObject* py_ari_func_{nullptr};
} ari_querier;

ARIQuerier::ARIQuerier() {
  Py_Initialize();
  py_module_name_ = PyUnicode_DecodeFSDefault("sklearn.metrics");
  py_module_ = PyImport_Import(py_module_name_);
  py_ari_func_ = PyObject_GetAttrString(py_module_, "adjusted_rand_score");
  ASSERT(py_module_name_ != nullptr);
  ASSERT(py_module_ != nullptr);
  ASSERT(py_ari_func_ != nullptr);
  ASSERT(PyCallable_Check(py_ari_func_));
}

ARIQuerier::~ARIQuerier() {
  Py_DECREF(py_ari_func_);
  Py_DECREF(py_module_);
  Py_DECREF(py_module_name_);
  Py_FinalizeEx();
}

double ARIQuerier::AdjustedRandIndex(
    const scan::Clustering& a,
    const scan::Clustering& b) const {
  PyObject* py_args{PyTuple_New(2)};
  for (const auto& [i, arg] : {std::make_pair(0, a), std::make_pair(1, b)}) {
    PyObject* py_list{PyList_New(arg.size())};
    for (size_t j{0}; j < arg.size(); j++) {
      PyObject* py_long{PyLong_FromLong(arg[j])};
      PyList_SetItem(py_list, j, py_long);
    }
    PyTuple_SetItem(py_args, i, py_list);
  }
  PyObject* py_return_value{PyObject_CallObject(py_ari_func_, py_args)};
  ASSERT(py_return_value != nullptr);
  const double adjusted_rand_index{PyFloat_AsDouble(py_return_value)};
  Py_DECREF(py_return_value);
  Py_DECREF(py_args);
  return adjusted_rand_index;
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

std::string ParametersToString(const size_t mu, const float epsilon) {
  return std::to_string(mu) + "-" + std::to_string(epsilon);
}

void PrintClock() {
  std::time_t time{std::chrono::system_clock::to_time_t(
      std::chrono::system_clock::now())};
  std::cerr << "[Clock] " << std::ctime(&time);
}

template <class Graph>
size_t GetMaxDegree(Graph* graph) {
  return parlay::reduce_max(parlay::delayed_seq<size_t>(
    graph->n,
    [&](const size_t i) {
      return graph->get_vertex(i).out_degree();
    }));
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
  std::cerr << "** Index construction median time: " << Median(times) << "\n\n";
  return index;
}

// Output query time as the median of several trials.
void OutputQueryTimes(
    const indexed_scan::Index& index,
    const size_t max_degree,
    const size_t num_rounds) {
  // Output query times with fixed epsilon and varying mu.
  for (size_t mu{2}; mu <= max_degree + 1; mu *= 2) {
    std::vector<double> query_times;
    constexpr double kEpsilon{0.6};
    if (verbose) { std::cerr << "query " << ParametersToString(mu, kEpsilon) << ":"; }
    for (size_t i{0}; i < num_rounds; i++) {
      timer query_timer{};
      // unused variable `clusters`, but keep it anyway so that it doesn't get
      // destructed within the timer
      const scan::Clustering clusters{index.Cluster(mu, kEpsilon)};
      query_times.emplace_back(query_timer.stop());
      if (verbose) { std::cerr << ' ' << query_times.back(); }
    }
    std::cerr << "** Cluster median time " << ParametersToString(mu, kEpsilon)
      << ": " << Median(query_times) << "\n";
  }

  // Output query times with fixed mu and varying epsilon.
  constexpr double kMu{5};
  const sequence<float> epsilons{.1, .2, .3, .4, .5, .6, .7, .8, .9};
  std::vector<double> total_query_times(num_rounds, 0);
  for (float epsilon : epsilons) {
    std::vector<double> query_times;
    if (verbose) { std::cerr << "query " << ParametersToString(kMu, epsilon) << ":"; }
    for (size_t i{0}; i < num_rounds; i++) {
      timer query_timer{"query"};
      const scan::Clustering clusters{index.Cluster(kMu, epsilon)};
      query_times.emplace_back(query_timer.stop());
      if (verbose) { std::cerr << ' ' << query_times.back(); }
      total_query_times[i] += query_times.back();
    }
    std::cerr << "** Cluster median time " << ParametersToString(kMu, epsilon)
      << ": " << Median(query_times) << "\n";
  }
  std::vector<double> bulk_query_times;
  constexpr auto clustering_no_op{[](scan::Clustering, size_t) {}};
  for (size_t i{0}; i < num_rounds; i++) {
    timer bulk_query_timer{"query"};
    index.Cluster(kMu, epsilons, clustering_no_op);
    bulk_query_times.emplace_back(bulk_query_timer.stop());
    if (verbose) { std::cerr << ' ' << bulk_query_times.back(); }
  }
  std::cerr << "** Bulk cluster median time mu=5: " << Median(bulk_query_times)
    << " vs. total individual queries: " << Median(total_query_times) << "\n";
}

// Output approximate index construction time as the median of several trials
// with different pseudorandom seeds.
template <class Graph, class SimilarityMeasure>
void OutputBuildApproximateIndexTimes(
    Graph* graph,
    const uint32_t num_samples,
    const size_t num_rounds) {
  std::vector<double> times;
  if (verbose) { std::cerr << "Index construction:"; }
  for (size_t i{0}; i < num_rounds; i++) {
    timer timer{"Index construction"};
    const indexed_scan::Index index{graph, SimilarityMeasure{num_samples, i}};
    times.emplace_back(timer.stop());
    if (verbose) { std::cerr << ' ' << times.back(); }
  }
  std::cerr << "** Index construction median time: " << Median(times) << "\n\n";
}

struct QueryInfo {
  uint64_t mu;
  float epsilon;
  scan::Clustering clusters;
  double modularity;
};

// Get a good clustering, where goodness is measured based on modularity,
// among parameters mu in {2, 4, 8, 16, ...} and epsilon in {.01, .02, .03, ...,
// .99}.
template <class Graph>
QueryInfo SearchForClusters(
    Graph* graph,
    const indexed_scan::Index index,
    const size_t max_degree) {
  const sequence<float> epsilons(
      99, [](const size_t i) { return (i + 1) * .01; });
  double best_modularity{-2.0};
  uint64_t best_mu{0};
  float best_epsilon{-1.0};
  scan::Clustering best_clusters{};
  for (size_t mu{2}; mu <= max_degree + 1; mu *= 2) {
    const auto update_best_clusters{
      [&](scan::Clustering clusters, size_t i) {
        const double modularity{scan::Modularity(graph, clusters)};
        if (modularity > best_modularity) {
          best_modularity = modularity;
          best_mu = mu;
          best_epsilon = epsilons[i];
          best_clusters = std::move(clusters);
        }
      }};
    constexpr bool kDeterministic{true};  // for consistency
    index.Cluster(mu, epsilons, update_best_clusters, kDeterministic);
  }
  return QueryInfo {
    .mu = best_mu,
    .epsilon = best_epsilon,
    .clusters = std::move(best_clusters),
    .modularity = best_modularity};
}

// Builds an index with an approximate similarity measure and outputs
//   - clustering modularity with approximate similarity measure at the best
//   parameters found by `SearchForClusters` relative to the exact similarity
//   measure,
//   - adjusted Rand index between approximate clustering and exact clustering
//   at the best parameters found by `SearchForClusters` relative to the exact
//   similarity measure,
//   - clustering modularity with approximate similarity measure at the best
//   parameters found by `SearchForClusters` relative to the approximate
//   measure.
// Measurements are the mean across several rounds with different pseudorandom
// seeds.
template <class Graph, class SimilarityMeasure>
void OutputApproximateQuality(
    Graph* graph,
    const QueryInfo& best_exact_clustering,
    const size_t max_degree,
    const uint32_t num_samples,
    const size_t num_rounds) {
  double total_modularity_at_exact = 0.0;
  double total_modularity_at_approx = 0.0;
  double total_ari = 0.0;
  for (size_t seed{0}; seed < num_rounds; seed++) {
    const indexed_scan::Index approx_index{
      graph, SimilarityMeasure{num_samples, seed}};
    {
      constexpr bool kDeterministic{true};  // for consistency
      const scan::Clustering clusters{approx_index.Cluster(
          best_exact_clustering.mu,
          best_exact_clustering.epsilon,
          kDeterministic)};
      const double modularity{scan::Modularity(graph, clusters)};
      const double ari{ari_querier.AdjustedRandIndex(
          best_exact_clustering.clusters, clusters)};
      total_modularity_at_exact += modularity;
      total_ari += ari;
      std::cerr << "At exact params: modularity=" << modularity
        << " ARI=" << ari << '\n';
    }
    {
      const QueryInfo best_approx_query{
        SearchForClusters(graph, approx_index, max_degree)};
      total_modularity_at_approx += best_approx_query.modularity;
      std::cerr << "At best approx params "
        << ParametersToString(best_approx_query.mu, best_approx_query.epsilon)
        << " modularity=" << best_approx_query.modularity << '\n';
    }
  }
  std::cerr << "Mean modularity @ exact params="
    << total_modularity_at_exact / num_rounds
    << " mean ari=" << total_ari / num_rounds << '\n';
  std::cerr << "Mean modularity @ approx best params="
    << total_modularity_at_approx / num_rounds << '\n';
}

// Run experiments with cosine similarity as the similarity measure.
template <class Graph>
void RunCosineExperiments(
    Graph* graph,
    const size_t num_index_rounds,
    const size_t num_cluster_rounds,
    const bool skip_approx_experiments) {
  const size_t max_degree{GetMaxDegree(graph)};
  std::cerr << "\n";
  std::cerr << "**********************\n";
  std::cerr << "** CosineSimilarity **\n";
  std::cerr << "**********************\n";

  PrintClock();
  indexed_scan::Index exact_index{BuildIndexAndOutputTimes(
      graph, scan::CosineSimilarity{}, num_index_rounds)};
  PrintClock();
  OutputQueryTimes(exact_index, max_degree, num_cluster_rounds);

  if (skip_approx_experiments) {
    return;
  }

  std::cerr << "****************************\n";
  std::cerr << "** ApproxCosineSimilarity **\n";
  std::cerr << "****************************\n";

  const QueryInfo best_exact_clustering{
    SearchForClusters(graph, exact_index, max_degree)};
  std::cerr << "Best exact params "
    << ParametersToString(
        best_exact_clustering.mu, best_exact_clustering.epsilon)
    << ", modularity " << best_exact_clustering.modularity << '\n';
  exact_index = indexed_scan::Index{};  // destruct index to save memory

  for (uint32_t num_samples{64}; 2 * num_samples < max_degree; num_samples *= 2) {
    std::cerr << "\n";
    std::cerr << "Samples=" << num_samples << '\n';
    std::cerr << "----------\n";
    PrintClock();
    OutputBuildApproximateIndexTimes<Graph, scan::ApproxCosineSimilarity>(
          graph, num_samples, num_index_rounds);
    PrintClock();
    OutputApproximateQuality<Graph, scan::ApproxCosineSimilarity>(
          graph,
          best_exact_clustering,
          max_degree,
          num_samples,
          num_index_rounds);
  }
}

// Run experiments with Jaccard similarity as the similarity measure.
template <class Graph>
void RunJaccardExperiments(
    Graph* graph,
    const size_t num_index_rounds,
    const bool skip_approx_experiments) {
  if (skip_approx_experiments) {
    return;
  }
  const size_t max_degree{GetMaxDegree(graph)};
  std::cerr << "\n";
  std::cerr << "*****************************\n";
  std::cerr << "** ApproxJaccardSimilarity **\n";
  std::cerr << "*****************************\n";

  indexed_scan::Index exact_index{graph, scan::JaccardSimilarity{}};
  const QueryInfo best_exact_clustering{
    SearchForClusters(graph, exact_index, max_degree)};
  std::cerr << "Best exact params "
    << ParametersToString(
        best_exact_clustering.mu, best_exact_clustering.epsilon)
    << ", modularity " << best_exact_clustering.modularity << '\n';
  exact_index = indexed_scan::Index{};  // destruct index to save memory

  for (uint32_t num_samples{64}; 2 * num_samples < max_degree; num_samples *= 2) {
    std::cerr << "\n";
    std::cerr << "Samples=" << num_samples << '\n';
    std::cerr << "----------\n";
    PrintClock();
    OutputBuildApproximateIndexTimes<Graph, scan::ApproxJaccardSimilarity>(
          graph, num_samples, num_index_rounds);
    PrintClock();
    OutputApproximateQuality<Graph, scan::ApproxJaccardSimilarity>(
          graph,
          best_exact_clustering,
          max_degree,
          num_samples,
          num_index_rounds);
  }
}

// Executes SCAN on the unweighted input graph and reports stats on the
// execution.
template <class Graph>
void RunScanUnweighted(Graph& graph, const commandLine& params) {
  {
    // get the graph into cache so timings are consistent
    indexed_scan::Index{&graph, scan::ApproxCosineSimilarity{1, 0}};
  }
  const bool is_serial{params.getOptionValue("--serial")};
  constexpr size_t kIndexRounds{5};
  constexpr size_t kClusterRounds{5};
  RunCosineExperiments(&graph, kIndexRounds, kClusterRounds, is_serial);
  RunJaccardExperiments(&graph, kIndexRounds, is_serial);
}

// Executes SCAN on the weighted input graph and reports stats on the execution.
template <class Graph>
void RunScanWeighted(Graph& graph, const commandLine& params) {
  {
    // get the graph into cache so timings are consistent
    indexed_scan::Index{&graph, scan::ApproxCosineSimilarity{1, 0}};
  }
  const bool is_serial{params.getOptionValue("--serial")};
  constexpr size_t kIndexRounds{5};
  constexpr size_t kClusterRounds{5};
  RunCosineExperiments(&graph, kIndexRounds, kClusterRounds, is_serial);
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
  if (is_graph_float_weighted) {
    auto graph{gbbs::gbbs_io::read_weighted_symmetric_graph<float>(
        input_graph_file, should_mmap_graph, is_graph_binary)};
    gbbs::RunScanWeighted(graph, params);
  } else {
    auto graph{gbbs::gbbs_io::read_unweighted_symmetric_graph(
        input_graph_file, should_mmap_graph, is_graph_binary)};
    gbbs::RunScanUnweighted(graph, params);
  }
}
