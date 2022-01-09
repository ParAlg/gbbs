#pragma once

#include "benchmarks/SpanningForest/LabelPropagation/SpanningForest.h"
#include "benchmarks/SpanningForest/LiuTarjan/SpanningForest.h"
#include "benchmarks/SpanningForest/ShiloachVishkin/SpanningForest.h"
#include "benchmarks/SpanningForest/UnionFind/SpanningForest.h"
#include "benchmarks/SpanningForest/common.h"

/* main connectit header */
#include "benchmarks/Connectivity/connectit.h"

/* sampling include */
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

template <UniteOption unite_option, class Find>
auto get_unite_function(size_t n, Find& find) {
  if
    constexpr(unite_option == unite) {
      return unite_variants::Unite<Find>(find);
    }
  else if
    constexpr(unite_option == unite_early) {
      return unite_variants::UniteEarly();
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
auto get_unite_function(size_t n, Find& find, Splice& splice) {
  if
    constexpr(unite_option == unite_rem_lock) {
      return unite_variants::UniteRemLock(find, splice, n);
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

/* Sampling algorithm returns */
template <class Graph, class Algorithm, AlgorithmType algorithm_type,
          SamplingOption sampling_option>
sequence<edge> compose_algorithm_and_sampling(Graph& G, commandLine& P,
                                              Algorithm& alg) {
  if
    constexpr(sampling_option == sample_kout) {
      using KOut = KOutSamplingTemplate<Graph>;
      auto sample = KOut(G, P);
      auto comp_alg =
          SamplingAlgorithmTemplate<Graph, KOut, Algorithm, algorithm_type,
                                    sampling_option>(G, sample, alg);
      return comp_alg.spanning_forest();
    }
  else if
    constexpr(sampling_option == sample_bfs) {
      using BFS = BFSSamplingTemplate<Graph>;
      auto sample = BFS(G, P);
      auto comp_alg =
          SamplingAlgorithmTemplate<Graph, BFS, Algorithm, algorithm_type,
                                    sampling_option>(G, sample, alg);
      return comp_alg.spanning_forest();
    }
  else if
    constexpr(sampling_option == sample_ldd) {
      using LDD = LDDSamplingTemplate<Graph>;
      auto sample = LDD(G, P);
      auto comp_alg =
          SamplingAlgorithmTemplate<Graph, LDD, Algorithm, algorithm_type,
                                    sampling_option>(G, sample, alg);
      return comp_alg.spanning_forest();
    }
  else {
    static_assert(sampling_option == no_sampling);
    auto comp_alg =
        NoSamplingAlgorithmTemplate<Graph, Algorithm, algorithm_type>(G, alg);
    return comp_alg.spanning_forest();
  }
}

/* Runs union-find algorithms (no-splice) */
template <class Graph, SamplingOption sampling_option, FindOption find_option,
          UniteOption unite_option>
sequence<edge> run_uf_alg(Graph& G, commandLine& P) {
  size_t n = G.n;
  auto find = get_find_function<find_option>();
  auto unite = get_unite_function<unite_option, decltype(find)>(n, find);
  using UF = union_find_sf::UFAlgorithm<decltype(find), decltype(unite), Graph>;
  auto alg = UF(G, unite, find);
  return compose_algorithm_and_sampling<Graph, decltype(alg), union_find_type,
                                        sampling_option>(G, P, alg);
}

/* Runs union-find algorithms (splice-based) */
template <class Graph, SamplingOption sampling_option, FindOption find_option,
          UniteOption unite_option, SpliceOption splice_option>
sequence<edge> run_uf_alg(Graph& G, commandLine& P) {
  static_assert(unite_option == unite_rem_cas ||
                unite_option == unite_rem_lock);
  auto find = get_find_function<find_option>();
  auto splice = get_splice_function<splice_option>();
  auto unite =
      get_unite_function<unite_option, decltype(find), decltype(splice),
                         find_option>(G.n, find, splice);
  using UF = union_find_sf::UFAlgorithm<decltype(find), decltype(unite), Graph>;
  auto alg = UF(G, unite, find);
  return compose_algorithm_and_sampling<Graph, decltype(alg), union_find_type,
                                        sampling_option>(G, P, alg);
}

/* Runs union-find algorithms (jayanti-tarjan) */
template <class Graph, SamplingOption sampling_option,
          JayantiFindOption find_option>
sequence<edge> run_jayanti_alg(Graph& G, commandLine& P) {
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
sequence<edge> run_sample_only_alg(Graph& G, commandLine& P) {
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
sequence<edge> run_liu_tarjan_alg(Graph& G, commandLine& P) {
  auto connect = lt::get_connect_function<connect_option>();
  auto update = lt::get_update_function<update_option>();
  auto shortcut = lt::get_shortcut_function<shortcut_option>();

  if
    constexpr(alter_option == no_alter) { /* no alter */
      using LT =
          lt::LiuTarjanAlgorithm<decltype(connect), connect_option,
                                 decltype(update), update_option,
                                 decltype(shortcut), shortcut_option, Graph>;
      auto alg = LT(G, G.n, connect, update, shortcut);

      return compose_algorithm_and_sampling<Graph, decltype(alg),
                                            liu_tarjan_type, sampling_option>(
          G, P, alg);
    }
  else { /* using alter: unsupported on cpu platform */
    abort();
  }
}

}  // namesapce connectit
}  // namespace gbbs
