#pragma once

namespace connectit {

  template <class Graph, bool provides_initial_graph>
  bool run_multiple_shiloach_vishkin(
      Graph& G,
      size_t n,
      pbbs::sequence<std::tuple<uintE, uintE>>& updates,
      size_t batch_size,
      size_t insert_to_query,
      size_t rounds,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine& P) {
      auto alg = shiloachvishkin_cc::SVAlgorithm<Graph>(G);
      return run_abstract_alg<Graph, decltype(alg), provides_initial_graph>(G, n, updates, batch_size, insert_to_query, alg);
    };

    auto name = "shiloach_vishkin";
    return run_multiple(G, rounds, name, P, test);
  }

} // connectit
