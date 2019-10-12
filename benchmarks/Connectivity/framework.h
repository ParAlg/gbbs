#pragma once

#include "UnionFind/Connectivity.h"
#include "ShiloachVishkin/Connectivity.h"
#include "LabelPropagation/Connectivity.h"

namespace connectit {


  template <
    class Graph,
    template <class F, class U, class G> class Sample,
    template <class Find> class Unite>
  pbbs::sequence<uintE> run_uf_algorithm(
      Graph& G,
      std::string& find,
      commandLine& P,
      bool use_hooks = false) {
  if (find_arg == "find_compress") {
    auto find = find_variants::find_compress;
    auto unite = Unite<decltype(find)>(G.n, find);
    auto sample = Sample<decltype(find), decltype(unite), Graph>(G, unite, find, P, use_hooks);
    auto algorithm = UFAlgorithm<decltype(find), decltype(unite)>(G, use_hooks);
    auto connectivity_algorithm = sampling_algorithm_template<Graph, decltype(sampler), decltype(algorithm)>();
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
    return q.components();
  } else if (find_arg == "find_naive") {
    auto find = find_variants::find_naive;
    auto unite = Unite<decltype(find)>(G.n, find);
    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
    return q.components();
  } else if (find_arg == "find_split") {
    auto find = find_variants::find_split;
    auto unite = Unite<decltype(find)>(G.n, find);
    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
    return q.components();
  } else if (find_arg == "find_halve") {
    auto find = find_variants::find_halve;
    auto unite = Unite<decltype(find)>(G.n, find);
    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
    return q.components();
  } else if (find_arg == "find_atomic_split") {
    auto find = find_variants::find_atomic_split;
    auto unite = Unite<decltype(find)>(G.n, find);
    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
    return q.components();
  } else if (find_arg == "find_atomic_halve") {
    auto find = find_variants::find_atomic_halve;
    auto unite = Unite<decltype(find)>(G.n, find);
    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
    return q.components();
  }
  return pbbs::sequence<uintE>();
}

  }

  template <class Graph,
            template <class F, class U, class G> class Sample>
  pbbs::sequence<uintE> run_uf_algorithm(
      Graph& G,
      std::string& unite,
      std::string& find,
      commandLine& P) {
    pbbs::sequence<uintE> components;
    if (unite == "unite") {
      components =
        union_find::select_framework_algorithm<
          unite_variants::Unite,
          Sample,
          Graph>(G, find, sampling_rounds, /* use_hooks = */false);
    } else if (unite == "unite_early") {
      components =
        union_find::select_framework_algorithm<
          unite_variants::UniteEarly,
          Sample,
          Graph>(G, find, sampling_rounds, /* use_hooks = */false);
    } else if (unite == "unite_nd") {
      components =
        union_find::select_framework_algorithm<
          unite_variants::UniteND,
          Sample,
          Graph>(G, find, sampling_rounds, /* use_hooks = */true);
    } else if (unite == "unite_rem_lock") {
      components =
        union_find::select_framework_algorithm<
          unite_variants::UniteRemLock,
          Sample,
          Graph>(G, find, sampling_rounds, /* use_hooks = */false);
    } else if (unite == "unite_rem_cas") {

    }
    return components;
  }

  template <class Graph>
  pbbs::sequence<uintE> run_algorithm(
      Graph& G,
      std::string& algorithm_type,
      std::string& sample,
      std::string& unite,
      std::string& find,
      commandLine& P) {
    if (algorithm_type == "uf") {
      if (sample == "kout") {
        return run_uf_algorithm<Graph, AfforestSamplingTemplate>(G, unite, find, P);
      } else if (sample == "bfs") {
        return run_uf_algorithm<Graph, union_find::UnionFindSampledBFSTemplate>(G, unite, find, P);
      } else if (sample == "ldd") {
        return run_uf_algorithm<Graph, union_find::UnionFindLDDTemplate>(G, unite, find, P);
      } else { // no sample
        return run_uf_algorithm<Graph, union_find::UnionFindTemplate>(G, unite, find, P);
      }
    } else if (algorithm_type == "liu_tarjan") {


    } else {
      std::cout << "Unknown algorithm type: " << algorithm_type << std::endl;
      abort();
    }
  }

} // namesapce connectit
