#pragma once

#include "benchmarks/Connectivity/UnionFind/Connectivity.h"
#include "benchmarks/Connectivity/ShiloachVishkin/Connectivity.h"
#include "benchmarks/Connectivity/LabelPropagation/Connectivity.h"
#include "benchmarks/Connectivity/LiuTarjan/Connectivity.h"
#include "benchmarks/Connectivity/common.h"

#include "sampling.h"

namespace connectit {

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
    if constexpr (sampling_option == sample_kout) {
      return "kout";
    } else if constexpr (sampling_option == sample_bfs) {
      return "bfs";
    } else if constexpr (sampling_option == sample_ldd) {
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
    }
  };


  template <SpliceOption splice_option>
  auto get_splice_function() {
    if constexpr (splice_option == split_atomic_one) {
      return splice_variants::split_atomic_one;
    } else if constexpr (splice_option == halve_atomic_one) {
      return splice_variants::halve_atomic_one;
    } else {
      return splice_variants::splice_atomic;
    }
  };

  template <UniteOption unite_option, class Find>
  auto get_unite_function(size_t n, Find& find) {
    if constexpr (unite_option == unite) {
      return unite_variants::Unite<Find>(find);
    } else if constexpr (unite_option == unite_early) {
      return unite_variants::UniteEarly();
    } else if constexpr (unite_option == unite_nd) {
      return unite_variants::UniteND<Find>(n, find);
    } else {
      std::cout << "Unsupported unite option" << std::endl;
      abort();
    }
  };

  template <UniteOption unite_option, class Find, class Splice, FindOption find_option>
  auto get_unite_function(size_t n, Find& find, Splice& splice) {
    if constexpr (unite_option == unite_rem_lock) {
      return unite_variants::UniteRemLock(find, splice, n);
    } else if constexpr (unite_option == unite_rem_cas) {
      return unite_variants::UniteRemCAS<decltype(splice), decltype(find), find_option>(find, splice);
    } else {
      std::cout << "Unsupported unite option" << std::endl;
      abort();
    }
  };


  template <
    class Graph,
    class Algorithm,
    AlgorithmType algorithm_type,
    SamplingOption sampling_option,
    FindOption find_option, /* for afforest */
    UniteOption unite_option /* for afforest */,
    SpliceOption splice_option /* for afforest */>
  pbbs::sequence<parent> compose_algorithm_and_sampling(Graph& G, commandLine& P, Algorithm& alg) {
    if constexpr (sampling_option == sample_kout) {
      auto find = get_find_function<find_option>();
      if constexpr (unite_option == unite_rem_cas || unite_option == unite_rem_lock) {
        auto splice = get_splice_function<splice_option>();
        auto unite = get_unite_function<unite_option, decltype(find), decltype(splice), find_option>(G.n, find, splice);
        using Afforest = AfforestSamplingTemplate<decltype(find), decltype(unite), Graph>;
        auto sample = Afforest(G, find, unite, P);
        auto connectivity = SamplingAlgorithmTemplate<Graph, Afforest, Algorithm, algorithm_type, sampling_option>(G, sample, alg);
        return connectivity.components();
      } else {
        auto unite = get_unite_function<unite_option, decltype(find)>(G.n, find);
        using Afforest = AfforestSamplingTemplate<decltype(find), decltype(unite), Graph>;
        auto sample = Afforest(G, find, unite, P);
        auto connectivity = SamplingAlgorithmTemplate<Graph, Afforest, Algorithm, algorithm_type, sampling_option>(G, sample, alg);
        return connectivity.components();
      }
    } else if constexpr (sampling_option == sample_bfs) {
      using BFS = BFSSamplingTemplate<Graph>;
      auto sample = BFS(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, BFS, Algorithm, algorithm_type, sampling_option>(G, sample, alg);
      return connectivity.components();
    } else if constexpr (sampling_option == sample_ldd) {
      using LDD = LDDSamplingTemplate<Graph>;
      auto sample = LDD(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, LDD, Algorithm, algorithm_type, sampling_option>(G, sample, alg);
      return connectivity.components();
    } else {
      static_assert(sampling_option == no_sampling);
      auto connectivity = NoSamplingAlgorithmTemplate<Graph, Algorithm, algorithm_type>(G, alg);
      return connectivity.components();
    }
  }


  /* Selects the sampling strategy, and calls the appropriate dispatcher */
  template <
    class Graph,
    SamplingOption sampling_option,
    FindOption find_option,
    UniteOption unite_option>
  pbbs::sequence<parent> run_uf_alg(
      Graph& G,
      commandLine& P) {
    size_t n = G.n;
    auto find = get_find_function<find_option>();
    auto unite = get_unite_function<unite_option, decltype(find)>(n, find);
    using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
    auto alg = UF(G, unite, find);
    return compose_algorithm_and_sampling<
      Graph,
      decltype(alg),
      union_find_type,
      sampling_option,
      find_option,
      unite_option,
      splice_atomic>(G, P, alg);
  }

  template <
    class Graph,
    SamplingOption sampling_option,
    FindOption find_option,
    UniteOption unite_option,
    SpliceOption splice_option>
  pbbs::sequence<parent> run_uf_alg(
      Graph& G,
      commandLine& P) {
    static_assert(unite_option == unite_rem_cas || unite_option == unite_rem_lock);
    auto find = get_find_function<find_option>();
    auto splice = get_splice_function<splice_option>();
    auto unite = get_unite_function<unite_option, decltype(find), decltype(splice), find_option>(G.n, find, splice);
    using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
    auto alg = UF(G, unite, find);
    return compose_algorithm_and_sampling<
      Graph,
      decltype(alg),
      union_find_type,
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
  pbbs::sequence<parent> run_jayanti_alg(
      Graph& G,
      commandLine& P) {
    size_t n = G.n;
    auto find = get_jayanti_find_function<find_option>();
    using UF = jayanti_rank::JayantiTBUnite<Graph, decltype(find)>;
    auto alg = UF(G, n, find);

    if constexpr (sampling_option == sample_kout) {
      auto fc = find_variants::find_compress;
      auto unite = unite_variants::UniteND<decltype(fc)>(n, fc);
      using Afforest = AfforestSamplingTemplate<decltype(fc), decltype(unite), Graph>;
      auto sample = Afforest(G, fc, unite, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, Afforest, UF, union_find_type, sampling_option>(G, sample, alg);
      return connectivity.components();
    } else if constexpr (sampling_option == sample_bfs) {
      using BFS = BFSSamplingTemplate<Graph>;
      auto sample = BFS(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, BFS, UF, union_find_type, sampling_option>(G, sample, alg);
      return connectivity.components();
    } else if constexpr (sampling_option == sample_ldd) {
      using LDD = LDDSamplingTemplate<Graph>;
      auto sample = LDD(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, LDD, UF, union_find_type, sampling_option>(G, sample, alg);
      return connectivity.components();
    } else {
      static_assert(sampling_option == no_sampling);
      auto connectivity = NoSamplingAlgorithmTemplate<Graph, UF, union_find_type>(G, alg);
      return connectivity.components();
    }
  }


  template <
    class Graph,
    SamplingOption sampling_option,
    template <class G> class Algorithm,
    AlgorithmType algorithm_type>
  pbbs::sequence<parent> run_sample_only_alg(
      Graph& G,
      commandLine& P) {
    using ALG = Algorithm<Graph>;
    auto alg = ALG(G);

    if constexpr (sampling_option == sample_kout) {
      auto fc = find_variants::find_compress;
      auto unite = unite_variants::UniteND<decltype(fc)>(G.n, fc);
      using Afforest = AfforestSamplingTemplate<decltype(fc), decltype(unite), Graph>;
      auto sample = Afforest(G, fc, unite, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, Afforest, ALG, algorithm_type, sampling_option>(G, sample, alg);
      return connectivity.components();
    } else if constexpr (sampling_option == sample_bfs) {
      using BFS = BFSSamplingTemplate<Graph>;
      auto sample = BFS(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, BFS, ALG, algorithm_type, sampling_option>(G, sample, alg);
      return connectivity.components();
    } else if constexpr (sampling_option == sample_ldd) {
      using LDD = LDDSamplingTemplate<Graph>;
      auto sample = LDD(G, P);
      auto connectivity = SamplingAlgorithmTemplate<Graph, LDD, ALG, algorithm_type, sampling_option>(G, sample, alg);
      return connectivity.components();
    } else {
      static_assert(sampling_option == no_sampling);
      auto connectivity = NoSamplingAlgorithmTemplate<Graph, ALG, algorithm_type>(G, alg);
      return connectivity.components();
    }
  }



  /* Selects the sampling strategy, and calls the appropriate dispatcher */

  template <LiuTarjanConnectOption connect_option>
  auto connect_to_string() {
    if constexpr (connect_option == simple_connect) {
      return "connect";
    } else if constexpr (connect_option == parent_connect) {
      return "parent_connect";
    } else if constexpr (connect_option == extended_connect) {
      return "extended_connect";
    } else {
      abort();
    }
  }

  template <LiuTarjanUpdateOption update_option>
  auto update_to_string() {
    if constexpr (update_option == simple_update) {
      return "simple_update";
    } else if constexpr (update_option == root_update) {
      return "root_update";
    } else {
      abort();
    }
  }

  template <LiuTarjanShortcutOption shortcut_option>
  auto shortcut_to_string() {
    if constexpr (shortcut_option == shortcut) {
      return "shortcut";
    } else if constexpr (shortcut_option == full_shortcut) {
      return "full_shortcut";
    } else {
      abort();
    }
  }

  template <LiuTarjanAlterOption alter_option>
  auto alter_to_string() {
    if constexpr (alter_option == alter) {
      return "alter";
    } else {
      return "no_alter";
    }
  }


  template <
    SamplingOption sampling_option,
    LiuTarjanConnectOption connect_option,
    LiuTarjanUpdateOption update_option,
    LiuTarjanShortcutOption shortcut_option,
    LiuTarjanAlterOption alter_option>
  std::string liu_tarjan_options_to_string() {
    return "liu_tarjan; sample=" + sampling_to_string<sampling_option>()
      + "; connect=" + connect_to_string<connect_option>()
      + "; update=" + update_to_string<update_option>()
      + "; shortcut=" + shortcut_to_string<shortcut_option>()
      + "; alter=" + alter_to_string<alter_option>();
  }

  template <
    class Graph,
    SamplingOption          sampling_option,
    LiuTarjanConnectOption  connect_option,
    LiuTarjanUpdateOption   update_option,
    LiuTarjanShortcutOption shortcut_option,
    LiuTarjanAlterOption    alter_option>
  pbbs::sequence<parent> run_liu_tarjan_alg(
      Graph& G,
      commandLine& P) {
    auto connect = lt::get_connect_function<connect_option>();
    auto update = lt::get_update_function<update_option>();
    auto shortcut = lt::get_shortcut_function<shortcut_option>();

    if constexpr (alter_option == no_alter) { /* no alter */
      using LT = lt::LiuTarjanAlgorithm<
        decltype(connect),
        connect_option,
        decltype(update),
        update_option,
        decltype(shortcut),
        shortcut_option,
        Graph>;
      auto alg = LT(G, G.n, connect, update, shortcut);

      return compose_algorithm_and_sampling<
        Graph,
        decltype(alg),
        liu_tarjan_type,
        sampling_option,
        find_compress,
        unite,
        splice_atomic>(G, P, alg);
    } else { /* using alter */
      abort();
    }
  }

} // namesapce connectit
