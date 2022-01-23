#pragma once

#include "benchmarks/Connectivity/common.h"

namespace gbbs {
namespace connectit {
template <class Graph, SamplingOption sampling_option, UniteOption unite_option,
          FindOption find_option>
bool run_multiple_uf_alg(Graph& G, size_t rounds, sequence<parent>& correct,
                         commandLine& P) {
  auto test = [&](Graph& graph, commandLine params,
                  sequence<parent>& correct_cc) {
    timer tt;
    tt.start();
    auto CC = run_uf_alg<Graph, sampling_option, find_option, unite_option>(
        graph, params);
    double t = tt.stop();
    if (params.getOptionValue("-check")) {
      cc_check(correct_cc, CC);
    }
    return t;
  };
  auto name = uf_options_to_string(sampling_option, find_option, unite_option);
  return run_multiple(G, rounds, correct, name, P, test);
}

template <class Graph, SamplingOption sampling_option, UniteOption unite_option,
          FindOption find_option, SpliceOption splice_option>
bool run_multiple_uf_alg(Graph& G, size_t rounds, sequence<parent>& correct,
                         commandLine& P) {
  auto test = [&](Graph& graph, commandLine params,
                  sequence<parent>& correct_cc) {
    timer tt;
    tt.start();
    auto CC = run_uf_alg<Graph, sampling_option, find_option, unite_option,
                         splice_option>(graph, params);
    double t = tt.stop();
    if (params.getOptionValue("-check")) {
      cc_check(correct_cc, CC);
    }
    return t;
  };
  auto name = uf_options_to_string(sampling_option, find_option, unite_option,
                                   splice_option);
  return run_multiple(G, rounds, correct, name, P, test);
}

template <class Graph, SamplingOption sampling_option,
          JayantiFindOption find_option>
bool run_multiple_jayanti_alg(Graph& G, size_t rounds,
                              sequence<parent>& correct, commandLine& P) {
  auto test = [&](Graph& graph, commandLine params,
                  sequence<parent>& correct_cc) {
    timer tt;
    tt.start();
    auto CC =
        run_jayanti_alg<Graph, sampling_option, find_option>(graph, params);
    double t = tt.stop();
    if (params.getOptionValue("-check")) {
      cc_check(correct_cc, CC);
    }
    return t;
  };
  auto name = jayanti_options_to_string(sampling_option, find_option);
  return run_multiple(G, rounds, correct, name, P, test);
}

}  // namespace connectit
}  // namespace gbbs
