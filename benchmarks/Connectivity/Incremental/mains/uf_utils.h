#pragma once

#include "benchmarks/Connectivity/connectit.h"

namespace gbbs {
namespace connectit {
template <class Graph, UniteOption unite_option, FindOption find_option,
          bool provides_initial_graph>
bool run_multiple_uf_alg(Graph& G, size_t n,
                         sequence<incremental_update>& updates,
                         size_t batch_size, size_t insert_to_query,
                         size_t rounds, commandLine& P) {
  auto test = [&](Graph& graph, commandLine params) {
    /* Create initial parents array */

    auto find_fn = get_find_function<find_option>();
    auto unite_fn =
        get_unite_function<unite_option, decltype(find_fn), find_option>(
            n, find_fn);
    using UF =
        union_find::UFAlgorithm<decltype(find_fn), decltype(unite_fn), Graph>;
    auto alg = UF(graph, unite_fn, find_fn);

    /* NOTE: none of the algorithms going through this step require reordering
     * a batch */
    static_assert(unite_option == unite || unite_option == unite_early ||
                  unite_option == unite_nd);

    bool check = params.getOptionValue("-check");
    return run_abstract_alg<Graph, decltype(alg), provides_initial_graph,
                            /* reorder_batch = */ false>(
        graph, n, updates, batch_size, insert_to_query, check, alg);
  };

  auto name = uf_options_to_string(no_sampling, find_option, unite_option);
  return run_multiple(G, rounds, name, P, test);
}

template <class Graph, UniteOption unite_option, FindOption find_option,
          SpliceOption splice_option, bool provides_initial_graph>
bool run_multiple_uf_alg(Graph& G, size_t n,
                         sequence<incremental_update>& updates,
                         size_t batch_size, size_t insert_to_query,
                         size_t rounds, commandLine& P) {
  static_assert(unite_option == unite_rem_cas ||
                unite_option == unite_rem_lock);

  auto test = [&](Graph& graph, commandLine params) {
    auto find = get_find_function<find_option>();
    auto splice = get_splice_function<splice_option>();
    auto unite =
        get_unite_function<unite_option, decltype(find), decltype(splice),
                           find_option>(n, find, splice);
    using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
    auto alg = UF(graph, unite, find);

    bool check = params.getOptionValue("-check");
    if
      constexpr(splice_option == splice_atomic) {
        return run_abstract_alg<Graph, decltype(alg), provides_initial_graph,
                                /* reorder_batch = */ true>(
            graph, n, updates, batch_size, insert_to_query, check, alg);
      }
    else {
      return run_abstract_alg<Graph, decltype(alg), provides_initial_graph,
                              /* reorder_batch = */ false>(
          graph, n, updates, batch_size, insert_to_query, check, alg);
    }
  };

  auto name = uf_options_to_string(no_sampling, find_option, unite_option,
                                   splice_option);
  return run_multiple(G, rounds, name, P, test);
}

template <class Graph, JayantiFindOption find_option,
          bool provides_initial_graph>
bool run_multiple_jayanti_alg(Graph& G, size_t n,
                              sequence<incremental_update>& updates,
                              size_t batch_size, size_t insert_to_query,
                              size_t rounds, commandLine& P) {
  auto test = [&](Graph& graph, commandLine& params) {
    auto find = get_jayanti_find_function<find_option>();
    using UF = jayanti_rank::JayantiTBUnite<Graph, decltype(find)>;
    auto alg = UF(graph, n, find);
    bool check = params.getOptionValue("-check");
    return run_abstract_alg<Graph, decltype(alg), provides_initial_graph,
                            /* reorder_batch = */ false>(
        graph, n, updates, batch_size, insert_to_query, check, alg);
  };

  auto name = jayanti_options_to_string(no_sampling, find_option);
  return run_multiple(G, rounds, name, P, test);
}
}  // namespace connectit
}  // namespace gbbs
