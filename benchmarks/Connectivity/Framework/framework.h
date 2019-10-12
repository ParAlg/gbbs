#pragma once

#include "benchmarks/Connectivity/UnionFind/Connectivity.h"
#include "benchmarks/Connectivity/ShiloachVishkin/Connectivity.h"
#include "benchmarks/Connectivity/LabelPropagation/Connectivity.h"

#include "sampling.h"

namespace connectit {

  enum SamplingOption {
    kout, bfs, ldd, no_sampling
  };

  enum FindOption {
    find_compress, find_naive, find_split, find_halve, find_atomic_split, find_atomic_halve
  };

  enum UniteOption {
    unite, unite_early, unite_nd, unite_rem_lock, unite_rem_cas
  };

  template <FindOption find_option>
  std::string find_to_string() {
    if constexpr (find_option == find_compress) {
      return "find_compress";
    } else if constexpr (find_option == find_naive) {
      return "find_naive";
    } else if constexpr (find_option == find_split) {
      return "find_split";
    } else if constexpr (find_option == find_halve) {
      return "find_halve";
    } else if constexpr (find_option == find_atomic_split) {
      return "find_atomic_split";
    } else if constexpr (find_option == find_atomic_halve) {
      return "find_atomic_halve";
    } else {
      abort();
      return "";
    }
  }

  template <UniteOption unite_option>
  std::string unite_to_string() {
    if constexpr (unite_option == unite) {
      return "unite";
    } else if constexpr (unite_option == unite_early) {
      return "unite_early";
    } else if constexpr (unite_option == unite_nd) {
      return "unite_nd";
    } else if constexpr (unite_option == unite_rem_lock) {
      return "unite_rem_lock";
    } else {
      abort();
    }
    return "";
  }

  template <SamplingOption sampling_option>
  std::string sampling_to_string() {
    if constexpr (sampling_option == kout) {
      return "kout";
    } else if constexpr (sampling_option == bfs) {
      return "bfs";
    } else if constexpr (sampling_option == ldd) {
      return "ldd";
    } else {
      return "no_sampling";
    }
  }

  template <SamplingOption sampling_option, FindOption find_option, UniteOption unite_option>
  std::string options_to_string() {
    return "uf; sample="
      + sampling_to_string<sampling_option>()
      + "; unite=" + unite_to_string<unite_option>()
      + "; find=" + find_to_string<find_option>();
  };

  template <FindOption find_option>
  auto get_find_function() {
    if constexpr (find_option == find_naive) {
      return find_variants::find_naive;
    } else if constexpr (find_option == find_compress) {
      return find_variants::find_compress;
    } else if constexpr (find_option == find_atomic_split) {
      return find_variants::find_atomic_split;
    } else if constexpr (find_option == find_atomic_halve) {
      return find_variants::find_atomic_halve;
    } else if constexpr (find_option == find_split) {
      return find_variants::find_split;
    } else if constexpr (find_option == find_halve) {
      return find_variants::find_halve;
    }
  };

  template <UniteOption unite_option, class Find>
  auto get_unite_function(size_t n, Find& find) {
    if constexpr (unite_option == unite) {
      return unite_variants::Unite<Find>(find);
    } else if constexpr (unite_option == unite_rem_lock) {
      return unite_variants::UniteRemLock(n);
    } else if constexpr (unite_option == unite_early) {
      return unite_variants::UniteEarly();
    } else if constexpr (unite_option == unite_nd) {
      return unite_variants::UniteND<Find>(n, find);
    } else {
      std::cout << "Unsupported unite option" << std::endl;
      abort();
    }
  };

  /* Selects the sampling strategy, and calls the appropriate dispatcher */
  template <
    class Graph,
    SamplingOption sampling_option,
    FindOption find_option,
    UniteOption unite_option>
  pbbs::sequence<uintE> run_uf_algorithm(
      Graph& G,
      commandLine& P) {
    size_t n = G.n;
    auto find = get_find_function<find_option>();
    auto unite = get_unite_function<unite_option, decltype(find)>(n, find);
    using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
    auto alg = UF(G, unite, find);

    if constexpr (sampling_option == kout) {
      using Afforest = AfforestSamplingTemplate<decltype(find), decltype(unite), Graph>;
      auto sample = Afforest(G, unite, find, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, Afforest, UF>(G, sample, alg);
      return connectivity.components();
    } else if constexpr (sampling_option == bfs) {
      using BFS = BFSSamplingTemplate<Graph>;
      auto sample = BFS(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, BFS, UF>(G, sample, alg);
      return connectivity.components();
    } else if constexpr (sampling_option == ldd) {
      using LDD = LDDSamplingTemplate<Graph>;
      auto sample = LDD(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, LDD, UF>(G, sample, alg);
      return connectivity.components();
    } else {
      static_assert(sampling_option == no_sampling);
      auto connectivity = NoSamplingAlgorithmTemplate<Graph, UF>(G, alg);
      return connectivity.components();
    }
  }

} // namesapce connectit
