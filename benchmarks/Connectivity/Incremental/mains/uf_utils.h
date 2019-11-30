#pragma once

namespace connectit {
  template<
    class Graph,
    UniteOption    unite_option,
    FindOption     find_option,
    bool provides_initial_graph>
  bool run_multiple_uf_alg(
      Graph& G,
      size_t n,
      pbbs::sequence<std::tuple<uintE, uintE>>& updates,
      size_t batch_size,
      size_t insert_to_query,
      size_t rounds,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine P) {
      /* Create initial parents array */

      auto find = get_find_function<find_option>();
      auto unite = get_unite_function<unite_option, decltype(find)>(n, find);
      using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
      auto alg = UF(G, unite, find);

      return run_abstract_alg<Graph, decltype(alg), provides_initial_graph>(G, n, updates, batch_size, insert_to_query, alg);
    };

    auto name = uf_options_to_string<no_sampling, find_option, unite_option>();
    return run_multiple(G, rounds, name, P, test);
  }

  template<
    class Graph,
    UniteOption    unite_option,
    FindOption     find_option,
    SpliceOption   splice_option,
    bool provides_initial_graph>
  bool run_multiple_uf_alg(
      Graph& G,
      size_t n,
      pbbs::sequence<std::tuple<uintE, uintE>>& updates,
      size_t batch_size,
      size_t insert_to_query,
      size_t rounds,
      commandLine& P) {
    static_assert(unite_option == unite_rem_cas || unite_option == unite_rem_lock);

    auto test = [&] (Graph& G, commandLine P) {
      auto find = get_find_function<find_option>();
      auto splice = get_splice_function<splice_option>();
      auto unite = get_unite_function<unite_option, decltype(find), decltype(splice), find_option>(G.n, find, splice);
      using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
      auto alg = UF(G, unite, find);

      return run_abstract_alg<Graph, decltype(alg), provides_initial_graph>(G, n, updates, batch_size, insert_to_query, alg);
    };

    auto name = uf_options_to_string<no_sampling, find_option, unite_option>();
    return run_multiple(G, rounds, name, P, test);
  }

  template <
    class Graph,
    JayantiFindOption find_option,
    bool provides_initial_graph>
  bool run_multiple_jayanti_alg(
      Graph& G,
      size_t n,
      pbbs::sequence<std::tuple<uintE, uintE>>& updates,
      size_t batch_size,
      size_t insert_to_query,
      size_t rounds,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine& P) {
      auto find = get_jayanti_find_function<find_option>();
      using UF = jayanti_rank::JayantiTBUnite<Graph, decltype(find)>;
      auto alg = UF(G, n, find);
      return run_abstract_alg<Graph, decltype(alg), provides_initial_graph>(G, n, updates, batch_size, insert_to_query, alg);
    };

    auto name = jayanti_options_to_string<no_sampling, find_option>();
    return run_multiple(G, rounds, name, P, test);
  }
} // connectit
