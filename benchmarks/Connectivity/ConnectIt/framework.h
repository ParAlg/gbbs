#pragma once

#include "benchmarks/Connectivity/LabelPropagation/Connectivity.h"
#include "benchmarks/Connectivity/LiuTarjan/Connectivity.h"
#include "benchmarks/Connectivity/ShiloachVishkin/Connectivity.h"
#include "benchmarks/Connectivity/UnionFind/Connectivity.h"
#include "benchmarks/Connectivity/common.h"

#include "sampling.h"

namespace gbbs {
namespace connectit {

/* Union-Find Strategies */
template <FindOption find_option>
auto get_find_function() {
  if
    constexpr(find_option == find_naive) { return find_variants::find_naive; }
  else if
    constexpr(find_option == find_compress) {
      return find_variants::find_compress;
    }
  else if
    constexpr(find_option == find_atomic_split) {
      return find_variants::find_atomic_split;
    }
  else if
    constexpr(find_option == find_atomic_halve) {
      return find_variants::find_atomic_halve;
    }
};

template <SpliceOption splice_option>
auto get_splice_function() {
  if
    constexpr(splice_option == split_atomic_one) {
      return splice_variants::split_atomic_one;
    }
  else if
    constexpr(splice_option == halve_atomic_one) {
      return splice_variants::halve_atomic_one;
    }
  else {
    return splice_variants::splice_atomic;
  }
};

template <UniteOption unite_option, class Find, FindOption find_option>
auto get_unite_function([[maybe_unused]] size_t n, Find& find) {
  if
    constexpr(unite_option == unite) {
      return unite_variants::Unite<Find>(find);
    }
  else if
    constexpr(unite_option == unite_early) {
      return unite_variants::UniteEarly<Find, find_option>(find);
    }
  else if
    constexpr(unite_option == unite_nd) {
      return unite_variants::UniteND<Find>(n, find);
    }
  else {
    std::cout << "Unsupported unite option" << std::endl;
    abort();
  }
};

template <UniteOption unite_option, class Find, class Splice,
          FindOption find_option>
auto get_unite_function([[maybe_unused]] size_t n, Find& find, Splice& splice) {
  if
    constexpr(unite_option == unite_rem_lock) {
      return unite_variants::UniteRemLock<decltype(splice), decltype(find),
                                          find_option>(find, splice, n);
    }
  else if
    constexpr(unite_option == unite_rem_cas) {
      return unite_variants::UniteRemCAS<decltype(splice), decltype(find),
                                         find_option>(find, splice);
    }
  else {
    std::cout << "Unsupported unite option" << std::endl;
    abort();
  }
};

/* Jayanti strategies */
template <JayantiFindOption find_option>
auto get_jayanti_find_function() {
  if
    constexpr(find_option == find_twotrysplit) {
      return jayanti_rank::find_twotry_splitting;
    }
  else {
    return jayanti_rank::find;
  }
}

template <class Graph, class Algorithm, AlgorithmType algorithm_type,
          SamplingOption sampling_option>
sequence<parent> compose_algorithm_and_sampling(Graph& G, commandLine& P,
                                                Algorithm& alg) {
  if
    constexpr(sampling_option == sample_kout) {
      using KOut = KOutSamplingTemplate<Graph>;
      auto sample = KOut(G, P);
      auto connectivity =
          SamplingAlgorithmTemplate<Graph, KOut, Algorithm, algorithm_type,
                                    sampling_option>(G, sample, alg);
      return connectivity.components();
    }
  else if
    constexpr(sampling_option == sample_bfs) {
      using BFS = BFSSamplingTemplate<Graph>;
      auto sample = BFS(G, P);
      auto connectivity =
          SamplingAlgorithmTemplate<Graph, BFS, Algorithm, algorithm_type,
                                    sampling_option>(G, sample, alg);
      return connectivity.components();
    }
  else if
    constexpr(sampling_option == sample_ldd) {
      using LDD = LDDSamplingTemplate<Graph>;
      auto sample = LDD(G, P);
      auto connectivity =
          SamplingAlgorithmTemplate<Graph, LDD, Algorithm, algorithm_type,
                                    sampling_option>(G, sample, alg);
      return connectivity.components();
    }
  else {
    static_assert(sampling_option == no_sampling);
    auto connectivity =
        NoSamplingAlgorithmTemplate<Graph, Algorithm, algorithm_type>(G, alg);
    return connectivity.components();
  }
}

/* Selects the sampling strategy, and calls the appropriate dispatcher */
template <class Graph, SamplingOption sampling_option, FindOption find_option,
          UniteOption unite_option>
sequence<parent> run_uf_alg(Graph& G, commandLine& P) {
  size_t n = G.n;
  auto find = get_find_function<find_option>();
  auto unite =
      get_unite_function<unite_option, decltype(find), find_option>(n, find);
  using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
  auto alg = UF(G, unite, find);
  return compose_algorithm_and_sampling<Graph, decltype(alg), union_find_type,
                                        sampling_option>(G, P, alg);
}

template <class Graph, SamplingOption sampling_option, FindOption find_option,
          UniteOption unite_option, SpliceOption splice_option>
sequence<parent> run_uf_alg(Graph& G, commandLine& P) {
  static_assert(unite_option == unite_rem_cas ||
                unite_option == unite_rem_lock);
  auto find = get_find_function<find_option>();
  auto splice = get_splice_function<splice_option>();
  auto unite =
      get_unite_function<unite_option, decltype(find), decltype(splice),
                         find_option>(G.n, find, splice);
  using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
  auto alg = UF(G, unite, find);
  return compose_algorithm_and_sampling<Graph, decltype(alg), union_find_type,
                                        sampling_option>(G, P, alg);
}

template <class Graph, SamplingOption sampling_option,
          JayantiFindOption find_option>
sequence<parent> run_jayanti_alg(Graph& G, commandLine& P) {
  size_t n = G.n;
  auto find = get_jayanti_find_function<find_option>();
  using UF = jayanti_rank::JayantiTBUnite<Graph, decltype(find)>;
  auto alg = UF(G, n, find);
  return compose_algorithm_and_sampling<Graph, decltype(alg), union_find_type,
                                        sampling_option>(G, P, alg);
}

/* Run a connectivity algorithm with no other internal options, switching only
 * on the sampling parameter (e.g., Shiloach-Vishkin) */
template <class Graph, SamplingOption sampling_option,
          template <class G> class Algorithm, AlgorithmType algorithm_type>
sequence<parent> run_sample_only_alg(Graph& G, commandLine& P) {
  using ALG = Algorithm<Graph>;
  auto alg = ALG(G);
  return compose_algorithm_and_sampling<Graph, decltype(alg), algorithm_type,
                                        sampling_option>(G, P, alg);
}

/* LiuTarjan: Selects the sampling strategy, and calls the appropriate
 * dispatcher */
template <
    class Graph, SamplingOption sampling_option,
    LiuTarjanConnectOption connect_option, LiuTarjanUpdateOption update_option,
    LiuTarjanShortcutOption shortcut_option, LiuTarjanAlterOption alter_option>
sequence<parent> run_liu_tarjan_alg(Graph& G, commandLine& P) {
  auto alg_connect = lt::get_connect_function<connect_option>();
  auto alg_update = lt::get_update_function<update_option>();
  auto alg_shortcut = lt::get_shortcut_function<shortcut_option>();
  auto alg_alter = lt::get_alter_function<alter_option>();

  if
    constexpr(alter_option == no_alter) { /* no alter */
      using LT =
          lt::LiuTarjanAlgorithm<decltype(alg_connect), connect_option,
                                 decltype(alg_update), update_option,
                                 decltype(alg_shortcut), shortcut_option,
                                 decltype(alg_alter), alter_option, Graph>;
      auto alg = LT(G, G.n, alg_connect, alg_update, alg_shortcut, alg_alter);

      return compose_algorithm_and_sampling<Graph, decltype(alg),
                                            liu_tarjan_type, sampling_option>(
          G, P, alg);
    }
  else { /* using alter: unsupported on cpu platform */
    abort();
  }
}

/* LiuTarjan: Selects the sampling strategy, and calls the appropriate
 * dispatcher */
template <
    class Graph, SamplingOption sampling_option,
    LiuTarjanConnectOption connect_option, LiuTarjanUpdateOption update_option,
    LiuTarjanShortcutOption shortcut_option, LiuTarjanAlterOption alter_option>
sequence<parent> run_liu_tarjan_alg(
    Graph& G, sequence<std::tuple<uintE, uintE, typename Graph::weight_type>>&&
                  mutable_graph,
    commandLine& P) {
  auto alg_connect = lt::get_connect_function<connect_option>();
  auto alg_update = lt::get_update_function<update_option>();
  auto alg_shortcut = lt::get_shortcut_function<shortcut_option>();
  auto alg_alter = lt::get_alter_function<alter_option>();

  using LT = lt::LiuTarjanAlgorithmCOO<
      decltype(alg_connect), connect_option, decltype(alg_update),
      update_option, decltype(alg_shortcut), shortcut_option,
      decltype(alg_alter), alter_option, typename Graph::weight_type>;
  auto alg = LT(std::move(mutable_graph), G.n, alg_connect, alg_update,
                alg_shortcut, alg_alter);

  return compose_algorithm_and_sampling<Graph, decltype(alg), liu_tarjan_type,
                                        sampling_option>(G, P, alg);
}

}  // namesapce connectit
}  // namespace gbbs
