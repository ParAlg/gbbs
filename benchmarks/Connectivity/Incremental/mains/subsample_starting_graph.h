#pragma once

#ifndef EMPTY_STARTING_GRAPH

template <class Graph, bool provides_initial_graph>
void run_all_tests(Graph& G, size_t n, pbbs::sequence<incremental_update>& updates, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P);

template <class Graph>
double Run(Graph& G, commandLine P) {
  using W = typename Graph::weight_type;
  int test_num = P.getOptionIntValue("-t", -1);
  int rounds = P.getOptionIntValue("-r", 5);

  double update_pct = P.getOptionDoubleValue("-update_pct", 0.1); /* percentage of edges to sample for updates */
  double deletion_pct = P.getOptionDoubleValue("-deletion_pct", 0.6) + update_pct; /* percentage of edges to delete */

  /* fraction of updates that are insertions */
  double insert_to_query = P.getOptionDoubleValue("-insert_to_query", 0.5);


  /* idea: hash every edge to a double in (0, 1). */

  size_t batch_size = P.getOptionLongValue("-batch_size", 1000000); /* batch size */

  /* 1) sample edges to use as updates */

  auto hash_to_double = [&] (const uintE& u, const uintE& v) -> double {
    auto min_v = std::min(u, v);
    auto max_v = std::max(u, v);
    size_t hashed_v = pbbs::hash64((static_cast<size_t>(min_v) << 32UL) + static_cast<size_t>(max_v));
    return static_cast<double>(hashed_v) / static_cast<double>(std::numeric_limits<size_t>::max());
  };

  auto update_pred =  [&] (const uintE& u, const uintE& v, const W& wgh) {
    return hash_to_double(u,v) < update_pct; /* return in sample */
  };
  auto updates_arr = sample_edges(G, update_pred);
  auto updates = pbbs::sequence((std::tuple<uintE, uintE>*)updates_arr.E, updates_arr.m);
  updates_arr.E = nullptr; /* relinquish memory */

  /* 2) call filter_graph to delete all deletions + updates from G */
  auto delete_pred =  [&] (const uintE& u, const uintE& v, const W& wgh) {
    return hash_to_double(u,v) > deletion_pct; /* keep */
  };
  /* Can also try using filter_edges */
  auto FG = filter_graph(G, delete_pred);
  std::cout << "### Initial graph size in edges = " << FG.m << std::endl;

  /* select a static_alg */

  auto annotated_updates = annotate_updates(updates, insert_to_query);

  size_t n = FG.n;
  run_all_tests<decltype(FG), true>(FG, n, annotated_updates, batch_size, insert_to_query, rounds, P);

  return 1.0;
}

//mmapcopy
generate_symmetric_once_main(Run, true);

#endif
