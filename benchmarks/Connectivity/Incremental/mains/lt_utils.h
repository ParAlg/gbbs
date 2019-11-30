#pragma once

namespace connectit {
  template <
    class Graph,
    LiuTarjanConnectOption  connect_option,
    LiuTarjanUpdateOption   update_option,
    LiuTarjanShortcutOption shortcut_option,
    LiuTarjanAlterOption    alter_option,
    bool provides_initial_graph>
  bool run_multiple_liu_tarjan_alg(
      Graph& G,
      size_t n,
      pbbs::sequence<std::tuple<uintE, uintE>>& updates,
      size_t batch_size,
      size_t insert_to_query,
      size_t rounds,
      commandLine& P) {

    auto test = [&] (Graph& G, commandLine& P) {
      auto connect = lt::get_connect_function<connect_option>();
      auto update = lt::get_update_function<update_option>();
      auto shortcut = lt::get_shortcut_function<shortcut_option>();

      static_assert(alter_option == no_alter);

      using LT = lt::LiuTarjanAlgorithm<
        decltype(connect),
        connect_option,
        decltype(update),
        update_option,
        decltype(shortcut),
        shortcut_option,
        Graph>;
      auto alg = LT(G, n, connect, update, shortcut);
      return run_abstract_alg<Graph, decltype(alg), provides_initial_graph>(G, n, updates, batch_size, insert_to_query, alg);
    };

    auto name = liu_tarjan_options_to_string<no_sampling,connect_option,update_option,shortcut_option,alter_option>();
    return run_multiple(G, rounds, name, P, test);
  }
} // connectit
