#pragma once

namespace gbbs {
namespace connectit {
template <class Graph, LiuTarjanConnectOption connect_option,
          LiuTarjanUpdateOption update_option,
          LiuTarjanShortcutOption shortcut_option,
          LiuTarjanAlterOption alter_option, bool provides_initial_graph>
bool run_multiple_liu_tarjan_alg(Graph& G, size_t n,
                                 sequence<incremental_update>& updates,
                                 size_t batch_size, size_t insert_to_query,
                                 size_t rounds, commandLine& P) {
  auto test = [&](Graph& graph, commandLine& params) {
    auto alg_connect = lt::get_connect_function<connect_option>();
    auto alg_update = lt::get_update_function<update_option>();
    auto alg_shortcut = lt::get_shortcut_function<shortcut_option>();
    auto alg_alter = lt::get_alter_function<alter_option>();

    using LT = lt::LiuTarjanAlgorithm<decltype(alg_connect), connect_option,
                                      decltype(alg_update), update_option,
                                      decltype(alg_shortcut), shortcut_option,
                                      decltype(alg_alter), alter_option, Graph>;
    auto alg = LT(graph, n, alg_connect, alg_update, alg_shortcut, alg_alter);
    bool check = params.getOptionValue("-check");
    return run_abstract_alg<Graph, decltype(alg), provides_initial_graph,
                            /* reorder_batch = */ true>(
        graph, n, updates, batch_size, insert_to_query, check, alg);
  };

  auto name =
      liu_tarjan_options_to_string(no_sampling, connect_option, update_option,
                                   shortcut_option, alter_option);
  return run_multiple(G, rounds, name, P, test);
}
}  // namespace connectit
}  // namespace gbbs
