#pragma once

#include "benchmarks/SpanningForest/check.h"

namespace gbbs {
namespace connectit {
template <class Graph, SamplingOption sampling_option, UniteOption unite_option,
          FindOption find_option>
bool run_multiple_uf_alg(Graph& G, size_t rounds, sequence<edge>& correct,
                         commandLine& P) {
  auto test = [&](Graph& G, commandLine P, sequence<edge>& correct) {
    timer tt;
    tt.start();
    auto edges =
        run_uf_alg<Graph, sampling_option, find_option, unite_option>(G, P);
    double t = tt.stop();
    if (P.getOptionValue("-check")) {
      spanning_forest::check_spanning_forest(G.n, correct, edges);
    }
    return t;
  };
  auto name =
      uf_options_to_string<sampling_option, find_option, unite_option>();
  return run_multiple(G, rounds, correct, name, P, test);
}

template <class Graph, SamplingOption sampling_option, UniteOption unite_option,
          FindOption find_option, SpliceOption splice_option>
bool run_multiple_uf_alg(Graph& G, size_t rounds, sequence<edge>& correct,
                         commandLine& P) {
  auto test = [&](Graph& G, commandLine P, sequence<edge>& correct) {
    timer tt;
    tt.start();
    auto edges = run_uf_alg<Graph, sampling_option, find_option, unite_option,
                            splice_option>(G, P);
    double t = tt.stop();
    if (P.getOptionValue("-check")) {
      spanning_forest::check_spanning_forest(G.n, correct, edges);
    }
    return t;
  };
  auto name = uf_options_to_string<sampling_option, find_option, unite_option,
                                   splice_option>();
  return run_multiple(G, rounds, correct, name, P, test);
}

template <class Graph, SamplingOption sampling_option,
          JayantiFindOption find_option>
bool run_multiple_jayanti_alg(Graph& G, size_t rounds, sequence<edge>& correct,
                              commandLine& P) {
  auto test = [&](Graph& G, commandLine P, sequence<edge>& correct) {
    timer tt;
    tt.start();
    auto edges = run_jayanti_alg<Graph, sampling_option, find_option>(G, P);
    double t = tt.stop();
    if (P.getOptionValue("-check")) {
      spanning_forest::check_spanning_forest(G.n, correct, edges);
    }
    return t;
  };
  auto name = jayanti_options_to_string<sampling_option, find_option>();
  return run_multiple(G, rounds, correct, name, P, test);
}

}  // namespace connectit
}  // namespace gbbs
