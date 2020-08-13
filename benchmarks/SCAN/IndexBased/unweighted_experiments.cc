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

#define PY_SSIZE_T_CLEAN
#include <python3.6m/Python.h>

#include "benchmarks/SCAN/IndexBased/scan.h"
#include "benchmarks/SCAN/IndexBased/similarity_measure.h"
#include "benchmarks/SCAN/IndexBased/utils.h"
#include "gbbs/gbbs.h"
#include "pbbslib/assert.h"

namespace gbbs {

namespace {

PyObject* py_module_name;
PyObject* py_module;
PyObject* py_ari_func;

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

double AdjustedRandIndex(const scan::Clustering& a, const scan::Clustering& b) {
  PyObject* py_args = PyTuple_New(2);
  for (const auto& [i, arg] : {std::make_pair(0, a), std::make_pair(1, b)}) {
    PyObject* py_list = PyList_New(arg.size());
    for (size_t j = 0; j < arg.size(); j++) {
      PyObject* py_long = PyLong_FromLong(arg[j]);
      PyList_SetItem(py_list, j, py_long);
    }
    PyTuple_SetItem(py_args, i, py_list);
  }
  PyObject* py_return_value = PyObject_CallObject(py_ari_func, py_args);
  ASSERT(py_return_value != nullptr);
  double adjusted_rand_index = PyFloat_AsDouble(py_return_value);
  Py_DECREF(py_return_value);
  Py_DECREF(py_args);
  return adjusted_rand_index;
}

struct QueryInfo {
  uint64_t mu;
  float epsilon;
  scan::Clustering clusters;
  double modularity;
};

// get a good clustering, where goodness is measured based on modularity
template <class Graph>
QueryInfo SearchForClusters(Graph* graph, const indexed_scan::Index index, const size_t max_degree) {
  double best_modularity{-2.0};
  uint64_t best_mu{0};
  float best_epsilon{-1.0};
  scan::Clustering best_clusters{};
  for (size_t mu{2}; mu < max_degree; mu *= 2) {
    for (double epsilon{0.01}; epsilon < 0.995; epsilon += 0.01) {
      constexpr bool kDeterministic{true};  // for consistency
      scan::Clustering clusters{index.Cluster(mu, epsilon, kDeterministic)};
      const double modularity{scan::Modularity(graph, clusters)};
      if (modularity > best_modularity) {
        best_modularity = modularity;
        best_mu = mu;
        best_epsilon = epsilon;
        best_clusters = std::move(clusters);
      }
    }
  }
  return QueryInfo {
    .mu = best_mu,
    .epsilon = best_epsilon,
    .clusters = std::move(best_clusters),
    .modularity = best_modularity};
}

template <class Graph>
std::vector<double> Modularities_5(Graph* graph, const indexed_scan::Index index) {
  constexpr uint64_t kMu{5};
  constexpr bool kDeterministic{true};  // for consistency
  std::vector<double> modularities;
  for (double epsilon{0.01}; epsilon < 0.995; epsilon += 0.01) {
    const scan::Clustering clusters{index.Cluster(kMu, epsilon, kDeterministic)};
    const double modularity{scan::Modularity(graph, clusters)};
    modularities.emplace_back(modularity);
  }
  return modularities;
}

}  // namespace

// Executes SCAN on the input graph and reports stats on the execution.
template <class Graph>
double RunScan(Graph& graph, commandLine parameters) {
  Py_Initialize();
  py_module_name = PyUnicode_DecodeFSDefault("sklearn.metrics");
  py_module = PyImport_Import(py_module_name);
  py_ari_func = PyObject_GetAttrString(py_module, "adjusted_rand_score");
  ASSERT(py_module_name != nullptr);
  ASSERT(py_module != nullptr);
  ASSERT(py_ari_func != nullptr);
  ASSERT(PyCallable_Check(py_ari_func));

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

    const QueryInfo best_exact_query{SearchForClusters(&graph, exact_index, max_degree)};
    std::cerr << "Best exact params " << ParametersToString(best_exact_query.mu, best_exact_query.epsilon) << ", modularity " << best_exact_query.modularity << '\n';
    const std::vector<double> exact_modularities_5{Modularities_5(&graph, exact_index)};
    std::cerr << "exact mod 5:";
    for (auto mod : exact_modularities_5) {
      std::cerr << " " << mod;
    }
    std::cerr << '\n';

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

      double total_modularity_at_exact = 0.0;
      double total_modularity_at_approx = 0.0;
      double total_ari = 0.0;
      std::vector<double> total_approx_modularities_5(exact_modularities_5.size());
      for (size_t seed{0}; seed < index_rounds; seed++) {
        const indexed_scan::Index approx_index{&graph, scan::ApproxCosineSimilarity{num_samples, seed}};
        {
          constexpr bool kDeterministic{true};  // for consistency
          const scan::Clustering clusters{approx_index.Cluster(best_exact_query.mu, best_exact_query.epsilon, kDeterministic)};
          const double modularity{scan::Modularity(&graph, clusters)};
          const double ari{AdjustedRandIndex(best_exact_query.clusters, clusters)};
          total_modularity_at_exact += modularity;
          total_ari += ari;
          std::cerr << "At exact params: modularity=" << modularity << " ARI=" << ari << '\n';
        }
        {
          const QueryInfo best_approx_query{SearchForClusters(&graph, approx_index, max_degree)};
          total_modularity_at_approx += best_approx_query.modularity;
          std::cerr << "At best approx params " << ParametersToString(best_approx_query.mu, best_approx_query.epsilon) << " modularity=" << best_approx_query.modularity << '\n';
        }
        const std::vector<double> approx_modularities_5{Modularities_5(&graph, approx_index)};
        for (size_t i = 0; i < approx_modularities_5.size(); i++) {
          total_approx_modularities_5[i] += approx_modularities_5[i];
        }
      }
      std::cerr << "Mean modularity @ exact params=" << total_modularity_at_exact / index_rounds << " mean ari=" << total_ari / index_rounds << '\n';
      std::cerr << "Mean modularity @ approx best params=" << total_modularity_at_approx / index_rounds << '\n';
      std::cerr << "approx mod 5:";
      for (auto mod : total_approx_modularities_5) {
        std::cerr << ' ' << mod / index_rounds;
      }
      std::cerr << '\n';
    }
  }

  {
    std::cerr << "*****************************\n";
    std::cerr << "** ApproxJaccardSimilarity **\n";
    std::cerr << "*****************************\n";

    const indexed_scan::Index exact_index{&graph, scan::JaccardSimilarity{}};
    const QueryInfo best_exact_query{SearchForClusters(&graph, exact_index, max_degree)};
    std::cerr << "Best exact params " << ParametersToString(best_exact_query.mu, best_exact_query.epsilon) << ", modularity " << best_exact_query.modularity << '\n';
    const std::vector<double> exact_modularities_5{Modularities_5(&graph, exact_index)};
    std::cerr << "exact mod 5:";
    for (auto mod : exact_modularities_5) {
      std::cerr << " " << mod;
    }
    std::cerr << '\n';

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

      double total_modularity_at_exact = 0.0;
      double total_modularity_at_approx = 0.0;
      double total_ari = 0.0;
      std::vector<double> total_approx_modularities_5(exact_modularities_5.size());
      for (size_t seed{0}; seed < index_rounds; seed++) {
        const indexed_scan::Index approx_index{&graph, scan::ApproxJaccardSimilarity{num_samples, seed}};
        {
          constexpr bool kDeterministic{true};  // for consistency
          const scan::Clustering clusters{approx_index.Cluster(best_exact_query.mu, best_exact_query.epsilon, kDeterministic)};
          const double modularity{scan::Modularity(&graph, clusters)};
          const double ari{AdjustedRandIndex(best_exact_query.clusters, clusters)};
          total_modularity_at_exact += modularity;
          total_ari += ari;
          std::cerr << "At exact params: modularity=" << modularity << " ARI=" << ari << '\n';
        }
        {
          const QueryInfo best_approx_query{SearchForClusters(&graph, approx_index, max_degree)};
          total_modularity_at_approx += best_approx_query.modularity;
          std::cerr << "At best approx params " << ParametersToString(best_approx_query.mu, best_approx_query.epsilon) << " modularity=" << best_approx_query.modularity << '\n';
        }
        const std::vector<double> approx_modularities_5{Modularities_5(&graph, approx_index)};
        for (size_t i = 0; i < approx_modularities_5.size(); i++) {
          total_approx_modularities_5[i] += approx_modularities_5[i];
        }
      }
      std::cerr << "Mean modularity @ exact params=" << total_modularity_at_exact / index_rounds << " mean ari=" << total_ari / index_rounds << '\n';
      std::cerr << "Mean modularity @ approx best params=" << total_modularity_at_approx / index_rounds << '\n';
      std::cerr << "approx mod 5:";
      for (auto mod : total_approx_modularities_5) {
        std::cerr << ' ' << mod / index_rounds;
      }
      std::cerr << '\n';
    }
  }

  Py_DECREF(py_ari_func);
  Py_DECREF(py_module);
  Py_DECREF(py_module_name);
  Py_FinalizeEx();
  return 0.0;
}

static constexpr bool kMutatesGraph{false};

}  // namespace gbbs

generate_symmetric_once_main(gbbs::RunScan, gbbs::kMutatesGraph);
