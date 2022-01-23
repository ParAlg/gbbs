#pragma once

// connectit header
#include "benchmarks/Connectivity/connectit.h"

// sampling strategies
#include "bfs_sampling.h"
#include "kout_sampling.h"
#include "ldd_sampling.h"

#include "utils.h"

/* ************************ Sampling Template ***************************/

namespace gbbs {
template <class Graph, class Sampler, class Algorithm,
          AlgorithmType algorithm_type, SamplingOption sampling_option>
struct SamplingAlgorithmTemplate {
  Graph& G;
  Sampler& sampler;
  Algorithm& algorithm;
  SamplingAlgorithmTemplate(Graph& G, Sampler& sampler, Algorithm& algorithm)
      : G(G), sampler(sampler), algorithm(algorithm) {}

  auto spanning_forest() {
    timer sample_t;
    sample_t.start();

    auto[Parents, Edges] = sampler.initial_spanning_forest();
    sample_t.stop();
    sample_t.next("sample time");

    parent frequent_comp;
    double pct;
    std::tie(frequent_comp, pct) = connectit::sample_frequent_element(Parents);

    algorithm.initialize(Parents, Edges);

    algorithm.template compute_spanning_forest<sampling_option>(Parents, Edges,
                                                                frequent_comp);

    return Edges;
    /* filter empty_edge pairs from Edge */
    return parlay::filter(Edges,
                          [&](const edge& e) { return e != empty_edge; });
  }
};

template <class Graph, class Algorithm, AlgorithmType algorithm_type>
struct NoSamplingAlgorithmTemplate {
  Graph& G;
  Algorithm& algorithm;
  NoSamplingAlgorithmTemplate(Graph& G, Algorithm& algorithm)
      : G(G), algorithm(algorithm) {}

  auto spanning_forest() {
    size_t n = G.n;
    auto Parents = sequence<parent>(n, [&](size_t i) { return i; });
    auto Edges = sequence<edge>(n, empty_edge);
    algorithm.initialize(Parents, Edges);
    algorithm.template compute_spanning_forest<no_sampling>(Parents, Edges);
    return Edges;
    return parlay::filter(Edges,
                          [&](const edge& e) { return e != empty_edge; });
  }
};

}  // namespace gbbs
