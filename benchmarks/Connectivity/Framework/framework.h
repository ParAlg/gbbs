#pragma once

#include "benchmarks/Connectivity/UnionFind/Connectivity.h"
#include "benchmarks/Connectivity/ShiloachVishkin/Connectivity.h"
#include "benchmarks/Connectivity/LabelPropagation/Connectivity.h"
#include "benchmarks/Connectivity/LiuTarjan/Connectivity.h"

#include "sampling.h"

namespace connectit {

  enum SamplingOption {
    kout, bfs, ldd, no_sampling
  };

  /* Union-Find options */
  enum FindOption {
    find_compress, find_naive, find_split, find_halve, find_atomic_split, find_atomic_halve
  };
  enum UniteOption {
    unite, unite_early, unite_nd, unite_rem_lock, unite_rem_cas
  };

  /* RemCAS-specific options */
  enum SpliceOption {
    split_atomic_one, halve_atomic_one, splice, splice_atomic
  };

  /* Jayanti-specific options */
  enum JayantiFindOption {
    find_twotrysplit, find_simple
  };

  /* LiuTarjan-specific options */

  using LiuTarjanConnectOption = liu_tarjan::LiuTarjanConnectOption;
  using LiuTarjanUpdateOption = liu_tarjan::LiuTarjanUpdateOption;
  using LiuTarjanShortcutOption = liu_tarjan::LiuTarjanShortcutOption;
  using LiuTarjanAlterOption = liu_tarjan::LiuTarjanAlterOption;

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

  template <SpliceOption splice_option>
  auto splice_to_string() {
    if constexpr (splice_option == split_atomic_one) {
      return "split_atomic_one";
    } else if constexpr (splice_option == halve_atomic_one) {
      return "halve_atomic_one";
    } else if constexpr (splice_option == splice) {
      return "splice";
    } else {
      return "splice_atomic";
    }
  };

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
    } else if constexpr (unite_option == unite_rem_cas) {
      return "unite_rem_cas";
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
  std::string uf_options_to_string() {
    return "uf; sample="
      + sampling_to_string<sampling_option>()
      + "; unite=" + unite_to_string<unite_option>()
      + "; find=" + find_to_string<find_option>();
  };


  template <SamplingOption sampling_option, FindOption find_option, UniteOption unite_option, SpliceOption splice_option>
  std::string uf_options_to_string() {
    return "uf; sample="
      + sampling_to_string<sampling_option>()
      + "; unite=" + unite_to_string<unite_option>()
      + "; find=" + find_to_string<find_option>() +
      + "; splice=" + splice_to_string<splice_option>();
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


  template <SpliceOption splice_option>
  auto get_splice_function() {
    if constexpr (splice_option == split_atomic_one) {
      return splice_variants::split_atomic_one;
    } else if constexpr (splice_option == halve_atomic_one) {
      return splice_variants::halve_atomic_one;
    } else if constexpr (splice_option == splice) {
      return splice_variants::splice;
    } else {
      return splice_variants::splice_atomic;
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


  template <
    class Graph,
    class Algorithm,
    SamplingOption sampling_option,
    FindOption find_option, /* for afforest */
    UniteOption unite_option /* for afforest */,
    SpliceOption splice_option /* for afforest */>
  pbbs::sequence<uintE> compose_algorithm_and_sampling(Graph& G, commandLine& P, Algorithm& alg) {
    size_t n = G.n;
    if constexpr (sampling_option == kout) {
      auto find = get_find_function<find_option>();
      if constexpr (unite_option == unite_rem_cas) {
        auto splice = get_splice_function<splice_option>();
        auto unite = unite_variants::UniteRemCAS<decltype(splice), decltype(find)>(splice, find);
        using Afforest = AfforestSamplingTemplate<decltype(find), decltype(unite), Graph>;
        auto sample = Afforest(G, find, unite, P);
        auto connectivity = SamplingAlgorithmTemplate<Graph, Afforest, Algorithm>(G, sample, alg);
        return connectivity.components();
      } else {
        auto unite = get_unite_function<unite_option, decltype(find)>(n, find);
        using Afforest = AfforestSamplingTemplate<decltype(find), decltype(unite), Graph>;
        auto sample = Afforest(G, find, unite, P);
        auto connectivity = SamplingAlgorithmTemplate<Graph, Afforest, Algorithm>(G, sample, alg);
        return connectivity.components();
      }
    } else if constexpr (sampling_option == bfs) {
      using BFS = BFSSamplingTemplate<Graph>;
      auto sample = BFS(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, BFS, Algorithm>(G, sample, alg);
      return connectivity.components();
    } else if constexpr (sampling_option == ldd) {
      using LDD = LDDSamplingTemplate<Graph>;
      auto sample = LDD(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, LDD, Algorithm>(G, sample, alg);
      return connectivity.components();
    } else {
      static_assert(sampling_option == no_sampling);
      auto connectivity = NoSamplingAlgorithmTemplate<Graph, Algorithm>(G, alg);
      return connectivity.components();
    }
  }


  /* Selects the sampling strategy, and calls the appropriate dispatcher */
  template <
    class Graph,
    SamplingOption sampling_option,
    FindOption find_option,
    UniteOption unite_option>
  pbbs::sequence<uintE> run_uf_alg(
      Graph& G,
      commandLine& P) {
    size_t n = G.n;
    auto find = get_find_function<find_option>();
    auto unite = get_unite_function<unite_option, decltype(find)>(n, find);
    using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
    auto alg = UF(G, unite, find);
    compose_algorithm_and_sampling<
      Graph,
      decltype(alg),
      sampling_option,
      find_option,
      unite_option,
      splice>(G, P, alg);
  }

  template <
    class Graph,
    SamplingOption sampling_option,
    FindOption find_option,
    UniteOption unite_option,
    SpliceOption splice_option>
  pbbs::sequence<uintE> run_uf_alg(
      Graph& G,
      commandLine& P) {
    static_assert(unite_option == unite_rem_cas);
    size_t n = G.n;
    auto find = get_find_function<find_option>();
    auto splice = get_splice_function<splice_option>();
    auto unite = unite_variants::UniteRemCAS<decltype(splice), decltype(find)>(splice, find);
    using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
    auto alg = UF(G, unite, find);
    compose_algorithm_and_sampling<
      Graph,
      decltype(alg),
      sampling_option,
      find_option,
      unite_option,
      splice_option>(G, P, alg);
  }


  /* Jayanti strategies */

  template <JayantiFindOption find_option>
  auto jayanti_find_to_string() {
    if constexpr (find_option == find_twotrysplit) {
      return "find_twotrysplitting";
    } else {
      return "find_simple";
    }
  }

  template <SamplingOption sampling_option, JayantiFindOption find_option>
  std::string jayanti_options_to_string() {
    return "jayanti; sample="
      + sampling_to_string<sampling_option>()
      + "; find=" + jayanti_find_to_string<find_option>();
  }

  template <JayantiFindOption find_option>
  auto get_jayanti_find_function() {
    if constexpr (find_option == find_twotrysplit) {
      return jayanti_rank::find_twotry_splitting;
    } else {
      return jayanti_rank::find;
    }
  }

  template <
    class Graph,
    SamplingOption sampling_option,
    JayantiFindOption find_option>
  pbbs::sequence<uintE> run_jayanti_alg(
      Graph& G,
      commandLine& P) {
    size_t n = G.n;
    auto find = get_jayanti_find_function<find_option>();
    using UF = jayanti_rank::JayantiTBUnite<Graph, decltype(find)>;
    auto alg = UF(G, find);

    if constexpr (sampling_option == kout) {
      auto fc = find_variants::find_compress;
      auto unite = unite_variants::UniteND<decltype(fc)>(n, fc);
      using Afforest = AfforestSamplingTemplate<decltype(fc), decltype(unite), Graph>;
      auto sample = Afforest(G, fc, unite, P);
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


  template <
    class Graph,
    SamplingOption sampling_option,
    template <class G> class Algorithm>
  pbbs::sequence<uintE> run_sample_only_alg(
      Graph& G,
      commandLine& P) {
    size_t n = G.n;
    using ALG = Algorithm<Graph>;
    auto alg = ALG(G);

    if constexpr (sampling_option == kout) {
      auto fc = find_variants::find_compress;
      auto unite = unite_variants::UniteND<decltype(fc)>(n, fc);
      using Afforest = AfforestSamplingTemplate<decltype(fc), decltype(unite), Graph>;
      auto sample = Afforest(G, fc, unite, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, Afforest, ALG>(G, sample, alg);
      return connectivity.components();
    } else if constexpr (sampling_option == bfs) {
      using BFS = BFSSamplingTemplate<Graph>;
      auto sample = BFS(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, BFS, ALG>(G, sample, alg);
      return connectivity.components();
    } else if constexpr (sampling_option == ldd) {
      using LDD = LDDSamplingTemplate<Graph>;
      auto sample = LDD(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, LDD, ALG>(G, sample, alg);
      return connectivity.components();
    } else {
      static_assert(sampling_option == no_sampling);
      auto connectivity = NoSamplingAlgorithmTemplate<Graph, ALG>(G, alg);
      return connectivity.components();
    }
  }


//  /* Selects the sampling strategy, and calls the appropriate dispatcher */
//  template <
//    class Graph,
//    LiuTarjanConnectOption connect_option,
//    LiuTarjanUpdateOption update_option,
//    LiuTarjanShortcutOption shortcut_option>
//  pbbs::sequence<uintE> run_liu_tarjan_alg(
//      Graph& G,
//      commandLine& P) {
//    size_t n = G.n;
//    auto connect = liu_tarjan::get_connect_function<connect_option>();
//    auto update = liu_tarjan::get_update_function<update_option>();
//    auto shortcut = liu_tarjan::get_shortcut_function<shortcut_option>();
//
//    using LT = union_find::LiuTarjanAlgorithm<
//      decltype(connect),
//      connect_option,
//      decltype(update),
//      update_option,
//      decltype(shortcut),
//      shortcut_option
//      Graph>;
//    auto alg = LT(G, connect, update, shortcut);
//
//    if constexpr (sampling_option == kout) {
//      using Afforest = AfforestSamplingTemplate<decltype(find), decltype(unite), Graph>;
//      auto sample = Afforest(G, find, unite, P);
//      auto connectivity = SamplingAlgorithmTemplate<Graph, Afforest, LT>(G, sample, alg);
//      return connectivity.components();
//    } else if constexpr (sampling_option == bfs) {
//      using BFS = BFSSamplingTemplate<Graph>;
//      auto sample = BFS(G, P);
//      auto connectivity = SamplingAlgorithmTemplate<Graph, BFS, LT>(G, sample, alg);
//      return connectivity.components();
//    } else if constexpr (sampling_option == ldd) {
//      using LDD = LDDSamplingTemplate<Graph>;
//      auto sample = LDD(G, P);
//      auto connectivity = SamplingAlgorithmTemplate<Graph, LDD, LT>(G, sample, alg);
//      return connectivity.components();
//    } else {
//      static_assert(sampling_option == no_sampling);
//      auto connectivity = NoSamplingAlgorithmTemplate<Graph, LT>(G, alg);
//      return connectivity.components();
//    }
//  }

} // namesapce connectit
