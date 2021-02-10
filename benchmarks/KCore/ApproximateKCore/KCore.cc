// Usage (static):
// numactl -i all ./KCore -rounds 3 -s -eps 0.4 -delta 3 <static graph>
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -eps: epsilon
//     -delta: delta
//     -use_pow: indicates whether to give less granular approximations
//               by using 2^{peeled_bkt} as the coreness estimate
//     -nb : the number of buckets to use in the bucketing implementation
//     -stats : indicates whether to output comparisons to exact coreness
//              values

// Usage (dynamic):
// numactl -i all ./KCore -rounds 3 -s -eps 0.4 -delta 3 -i <dynamic graph> -start_size 0 -end_size 10000 -b 1000 <static graph>
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -i : file path to the dynamic graph, with a single dynamic edge on each
//          line, in the format "<+/-> u v", where + denotes insertion, -
//          denotes deletion, and (u, v) are the vertex ids denoting the edge
//          (e.g., "+ 0 1")
//     -start_size : indicates the index of the dynamic edge given by -i
//                   in which to begin performing dynamic updates (dynamic edges
//                   prior to this index are processed, but the time
//                   taken to load the dynamic edges prior to this index is
//                   not included in the batch time in the output)
//     -end_size : indicates the index of the dynamic edge given by -i in which
//                 to stop performing dynamic updates (exclusive); default
//                 value is the last edge in the dynamic graph
//     -b : batch size (output includes time per batch, from start_size to
//          end_size)
//     -eps: epsilon
//     -delta: delta
//     -use_pow: indicates whether to give less granular approximations
//               by using 2^{peeled_bkt} as the coreness estimate
//     -nb : the number of buckets to use in the bucketing implementation
//     -stats : indicates whether to output comparisons to exact coreness
//              values
// Note: The static graph is ignored if -i is specified.

#include "KCore.h"
#include "gbbs/dynamic_graph_io.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

namespace gbbs {

// Given approximate cores, output comparisons to exact k-core
template <class Graph>
void print_stats(sequence<uintE>& cores, double eps, Graph& dynamic_graph,
  size_t num_buckets, bool use_stats){
  if (use_stats) {
    uintE max_core = parlay::reduce(cores, parlay::maxm<uintE>());
    std::cout << "### Coreness Estimate: " << max_core << std::endl;

    double mult_appx = (2 + 2*eps);
    auto true_cores = KCore(dynamic_graph, num_buckets);
    auto n = dynamic_graph.n;

    uintE max_true_core = parlay::reduce(true_cores, parlay::maxm<uintE>());
    std::cout << "### Coreness Exact: " << max_true_core << std::endl;

    double total_error = 0.0;
    double max_error = 0.0;
    double min_error = std::numeric_limits<double>::max();
    double denominator = 0;
    for (size_t i=0; i<n; i++) {
      double true_core = true_cores[i];
      double appx_core = cores[i];
      if (true_core != 0 && appx_core != 0) {
        auto this_error = std::max(true_core, appx_core) /
          std::min(true_core, appx_core);
        total_error += this_error;
        max_error = std::max(max_error, this_error);
        min_error = std::min(min_error, this_error);
        denominator++;
      }
    }
    min_error = min_error == std::numeric_limits<double>::max() ? 0 : min_error;
    auto avg_error = (denominator == 0) ? 0 : (total_error / denominator);

    std::cout << "### Per Vertex Average Coreness Error: " << avg_error
      << std::endl;
    std::cout << "### Per Vertex Min Coreness Error: " << min_error << std::endl;
    std::cout << "### Per Vertex Max Coreness Error: " << max_error << std::endl;
  }
}

template <class Graph>
double KCore_runner(Graph& G, commandLine P) {
  size_t num_buckets = P.getOptionLongValue("-nb", 16);
  double eps = P.getOptionDoubleValue("-eps", 0.2);
  double delta = P.getOptionDoubleValue("-delta", 0.1);
  bool use_pow = P.getOptionValue("-use_pow");
  std::cout << "### Application: KCore" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  if (num_buckets !=
    static_cast<size_t>((1 << pbbslib::log2_up(num_buckets)))) {
    std::cout << "Number of buckets must be a power of two."
              << "\n";
    exit(-1);
  }
  assert(P.getOption("-s"));

  // Options to use the dynamic graph instead of the static graph
  const std::string kInputFlag{"-i"};
  const char* const input_file{P.getOptionValue(kInputFlag)};
  bool use_stats = P.getOptionValue("-stats");
  long start_size = P.getOptionLongValue("-start_size", 0);
  long end_size = P.getOptionLongValue("-end_size", 0);

  // Options not exposed to general users
  long num_dynamic_edges = P.getOptionLongValue("-num_dynamic_edges", 0);
  long batch_size = P.getOptionLongValue("-b", 1);
  const char* const init_graph_file(P.getOptionValue("-init_graph_file"));

  bool use_dynamic = (input_file && input_file[0]);

  // Run the static algorithm
  if (!use_dynamic) {
    timer t; t.start();
    auto cores = approximate_kcore::KCore(G, num_buckets, eps, delta, use_pow);
    double tt = t.stop();
    if (use_stats) print_stats(cores, eps, G, num_buckets, use_stats);
    std::cout << "### Running Time: " << tt << std::endl;
    return tt;
  }
  
  // Load dynamic edges
  using W = typename Graph::weight_type;
  BatchDynamicEdges<W> batch_edge_list =
    read_batch_dynamic_edge_list<W>(input_file);
  
  // Load initial graph if it exists
  size_t offset = 0;
  if (init_graph_file) {
    BatchDynamicEdges<W> init_graph_list =
      read_batch_dynamic_edge_list<W>(init_graph_file);
    offset = prepend_dynamic_edge_list(batch_edge_list, init_graph_list);
  }

  // Run the dynamic algorithm given num_dynamic_edges and start_size(not
  // exposed to general users)
  if (num_dynamic_edges != 0) {
    auto batch = batch_edge_list.edges;

    // Compute correct end index of the batch
    num_dynamic_edges += offset;
    if (num_dynamic_edges > batch.size()) num_dynamic_edges = batch.size();

    // Compute the correct start index of the batch
    start_size += offset;

    // Separate insertions and deletions in the batch
    auto compare_seq = parlay::delayed_seq<bool>(
      num_dynamic_edges - start_size, [&] (size_t i) {
        return !batch[i + start_size].insert;
    });
    auto split = parlay::internal::split_two(
      parlay::make_slice(batch.data() + start_size,
      batch.data() + num_dynamic_edges), compare_seq);
    // Store the index of the first deletion in the batch
    size_t j = split.second + start_size;
    parallel_for(start_size, num_dynamic_edges, [&] (size_t i) {
      batch[i] = split.first[i - start_size];
    });

    timer t; t.start();
    // Insert the batch if there exist insertions
    if (j != start_size) {
      symmetric_graph<symmetric_vertex, W> dynamic_graph_insert =
        dynamic_edge_list_to_symmetric_graph(batch_edge_list, j);
      auto cores_insert = approximate_kcore::KCore(
        dynamic_graph_insert, num_buckets, eps, delta, use_pow);
      print_stats(
        cores_insert, eps, dynamic_graph_insert, num_buckets, use_stats);
    }

    // Delete the batch if there exist deletions
    if (j != num_dynamic_edges){
      symmetric_graph<symmetric_vertex, W> dynamic_graph = 
        dynamic_edge_list_to_symmetric_graph(batch_edge_list, num_dynamic_edges);
      auto cores = approximate_kcore::KCore(
        dynamic_graph, num_buckets, eps, delta, use_pow);
      print_stats(cores, eps, dynamic_graph, num_buckets, use_stats);
    }
    double tt = t.stop();

    std::cout << "### Batch Running Time: " << tt << std::endl;
    std::cout << "### Batch Num: " << num_dynamic_edges - offset << std::endl;
    return tt;
  }

  // Run the dynamic algorithm given start_size, end_size, and batch_size
  auto batch = batch_edge_list.edges;

  // Compute the correct end index of the batch
  end_size += offset;
  if (end_size == offset || end_size > batch.size()) end_size = batch.size();

  // Compute the correct start index of the batch
  start_size += offset;

  // Separate the insertions and deletions in each batch of size batch_size
  // between start_size and end_size
  auto end_seq = parlay::sequence<uintE>(3 + 2 * ((end_size - start_size) /
    batch_size), (uintE)0);
  end_seq[0] = start_size;
  for (size_t i = start_size; i < end_size; i += batch_size) {
    // Compute the index immediately following the last edge in the batch
    num_dynamic_edges = std::min((size_t) end_size, i + batch_size);

    auto compare_seq = parlay::delayed_seq<bool>(
      num_dynamic_edges - i, [&] (size_t k) {
        return !batch[k + i].insert;
    });
    auto split = parlay::internal::split_two(parlay::make_slice(
      batch.data() + i, batch.data() + num_dynamic_edges), compare_seq);
  
    // Store the index of the first deletion in the batch
    size_t j = split.second + i;
    parallel_for(i, num_dynamic_edges, [&] (size_t k) {
      batch[k] = split.first[k - i];
    });

    // Store in order the indices of the first deletion and insertion, for each
    // batch
    end_seq[1 + 2 * ((i - start_size) / batch_size)] = (j == i) ?
      UINT_E_MAX : j;
    end_seq[1 + 2 * ((i - start_size) / batch_size) + 1] =
      (j == num_dynamic_edges) ? UINT_E_MAX : num_dynamic_edges;
  }

  timer t1; t1.start();
  for (size_t i = 0; i < end_seq.size() - 1; i++) {
    // Retrieve the indices separating insertions and deletions
    size_t start = end_seq[i];
    size_t end = end_seq[i + 1];
    if (i != 0 && (start == 0 || end == 0)) break;

    timer t; t.start();
    // Process the batch if the batch is of non-zero size
    if (end != UINT_E_MAX) {
      auto dynamic_graph = dynamic_edge_list_to_symmetric_graph(batch_edge_list, end);
      auto cores = approximate_kcore::KCore(dynamic_graph, num_buckets, eps, delta, use_pow);
      print_stats(cores, eps, dynamic_graph, num_buckets, use_stats);
    }
    double tt = t.stop();

    // Output the batch running time and index if the batch is of non-zero
    // size
    if (end != UINT_E_MAX) {
      std::cout << "### Batch Running Time: " << tt << std::endl;
      if ((end - offset) % batch_size == 0 || end == batch.size())
        std::cout << "### Batch Num: " << end - offset << std::endl;
    }
  }
  double tt1 = t1.stop();
  
  return tt1;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::KCore_runner, false);
