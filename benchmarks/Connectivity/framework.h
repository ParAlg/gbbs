#pragma once

#include "UnionFind/Connectivity.h"
#include "ShiloachVishkin/Connectivity.h"
#include "LabelPropagation/Connectivity.h"

namespace connectit {
  template <class Graph,
            template <class F, class U, class G> class Sample>
  pbbs::sequence<uintE> run_algorithm_sampler(
      Graph& G,
      std::string& unite,
      std::string& find,
      commandLine& P) {
    size_t sampling_rounds = P.getOptionLongValue("-sample_rounds", 2L);
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
    if (algorithm_type == "union_find") {
      if (sample == "kout") {
        return run_algorithm_sampler<Graph, union_find::UnionFindSampleTemplate>(G, unite, find, P);
      } else if (sample == "bfs") {
        return run_algorithm_sampler<Graph, union_find::UnionFindSampledBFSTemplate>(G, unite, find, P);
      } else if (sample == "ldd") {
        return run_algorithm_sampler<Graph, union_find::UnionFindLDDTemplate>(G, unite, find, P);
      } else { // no sample
        return run_algorithm_sampler<Graph, union_find::UnionFindTemplate>(G, unite, find, P);
      }
    } else if (algorithm_type == "liu_tarjan") {

    } else {
      std::cout << "Unknown algorithm type: " << algorithm_type << std::endl;
      abort();
    }
  }

} // namesapce connectit
