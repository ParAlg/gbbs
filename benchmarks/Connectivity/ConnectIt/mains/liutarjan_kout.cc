// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "benchmarks/Connectivity/BFSCC/Connectivity.h"
#include "benchmarks/Connectivity/ConnectIt/framework.h"
#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"
#include "benchmarks/Connectivity/common.h"

#include "bench_utils.h"
#include "uf_utils.h"

namespace gbbs {
namespace connectit {

template <class Graph, SamplingOption sampling_option,
          LiuTarjanConnectOption connect_option,
          LiuTarjanUpdateOption simple_update_option,
          LiuTarjanShortcutOption shortcut_option,
          LiuTarjanAlterOption alter_option>
bool run_multiple_liu_tarjan_alg(Graph& G, size_t rounds,
                                 sequence<parent>& correct, commandLine& P) {
  if
    constexpr(alter_option == no_alter) {
      auto test = [&](Graph& graph, commandLine params,
                      sequence<parent>& correct_cc) {
        timer tt;
        tt.start();
        auto CC = run_liu_tarjan_alg<Graph, sampling_option, connect_option,
                                     simple_update_option, shortcut_option,
                                     alter_option>(graph, params);
        double t = tt.stop();
        if (params.getOptionValue("-check")) {
          cc_check(correct_cc, CC);
        }
        return t;
      };
      auto name = liu_tarjan_options_to_string(sampling_option, connect_option,
                                               simple_update_option,
                                               shortcut_option, alter_option);
      return run_multiple(G, rounds, correct, name, P, test);
    }
  else {
    auto test = [&](Graph& graph, commandLine params,
                    sequence<parent>& correct_cc) {
      if (graph.m > 10000000000UL) {
        return 0.0; /* too large to run in COO */
      }
      // Convert graph to COO format:
      auto mutable_edges = graph.edges();
      // using W = typename Graph::weight_type;
      // if (sizeof(W) > 0) {
      //  std::cout << "Do not use with a weighted graph!" << std::endl; //
      //  (trivial to in fact use it with a weighted graph, but avoid for
      //  benchmarking as it adds orthogonal overheads)
      //  abort();
      //}
      timer tt;
      tt.start();
      auto CC = run_liu_tarjan_alg<Graph, sampling_option, connect_option,
                                   simple_update_option, shortcut_option,
                                   alter_option>(
          graph, std::move(mutable_edges), params);
      double t = tt.stop();
      if (params.getOptionValue("-check")) {
        cc_check(correct_cc, CC);
      }
      return t;
    };
    auto name = liu_tarjan_options_to_string(sampling_option, connect_option,
                                             simple_update_option,
                                             shortcut_option, alter_option);
    return run_multiple(G, rounds, correct, name, P, test);
  }
}

template <class Graph>
void liutarjan_CUSA(Graph& G, int rounds, commandLine& P,
                    sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, simple_connect, simple_update,
                              shortcut, alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_CRSA(Graph& G, int rounds, commandLine& P,
                    sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, simple_connect, root_update,
                              shortcut, alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_PUSA(Graph& G, int rounds, commandLine& P,
                    sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, simple_update,
                              shortcut, alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_PRSA(Graph& G, int rounds, commandLine& P,
                    sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, root_update,
                              shortcut, alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_PUS(Graph& G, int rounds, commandLine& P,
                   sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, simple_update,
                              shortcut, no_alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_PRS(Graph& G, int rounds, commandLine& P,
                   sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, root_update,
                              shortcut, no_alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_EUSA(Graph& G, int rounds, commandLine& P,
                    sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, extended_connect,
                              simple_update, shortcut, alter>(G, rounds,
                                                              correct, P);
}

template <class Graph>
void liutarjan_EUS(Graph& G, int rounds, commandLine& P,
                   sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, extended_connect,
                              simple_update, shortcut, no_alter>(G, rounds,
                                                                 correct, P);
}

template <class Graph>
void liutarjan_CUFA(Graph& G, int rounds, commandLine& P,
                    sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, simple_connect, simple_update,
                              full_shortcut, alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_CRFA(Graph& G, int rounds, commandLine& P,
                    sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, simple_connect, root_update,
                              full_shortcut, alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_PUFA(Graph& G, int rounds, commandLine& P,
                    sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, simple_update,
                              full_shortcut, alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_PRFA(Graph& G, int rounds, commandLine& P,
                    sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, root_update,
                              full_shortcut, alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_PUF(Graph& G, int rounds, commandLine& P,
                   sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, simple_update,
                              full_shortcut, no_alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_PRF(Graph& G, int rounds, commandLine& P,
                   sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, root_update,
                              full_shortcut, no_alter>(G, rounds, correct, P);
}

template <class Graph>
void liutarjan_EUFA(Graph& G, int rounds, commandLine& P,
                    sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, extended_connect,
                              simple_update, full_shortcut, alter>(G, rounds,
                                                                   correct, P);
}

template <class Graph>
void liutarjan_EUF(Graph& G, int rounds, commandLine& P,
                   sequence<parent>& correct) {
  run_multiple_liu_tarjan_alg<Graph, sample_kout, extended_connect,
                              simple_update, full_shortcut, no_alter>(
      G, rounds, correct, P);
}
}

template <class Graph>
double Benchmark_runner(Graph& G, commandLine P) {
  int rounds = P.getOptionIntValue("-r", 5);

  auto correct = sequence<parent>();
  if (P.getOptionValue("-check")) {
    correct = workefficient_cc::CC(G, 0.2, false, true);
    RelabelDet(correct);
  }
  run_tests(
      G, rounds, P, correct, connectit::liutarjan_PUSA<Graph>,
      {
          connectit::liutarjan_CUSA<Graph>, connectit::liutarjan_CRSA<Graph>,
          connectit::liutarjan_PUSA<Graph>, connectit::liutarjan_PRSA<Graph>,
          connectit::liutarjan_PUS<Graph>, connectit::liutarjan_PRS<Graph>,
          connectit::liutarjan_EUSA<Graph>, connectit::liutarjan_EUS<Graph>,

          connectit::liutarjan_CUFA<Graph>, connectit::liutarjan_CRFA<Graph>,
          connectit::liutarjan_PUFA<Graph>, connectit::liutarjan_PRFA<Graph>,
          connectit::liutarjan_PUF<Graph>, connectit::liutarjan_PRF<Graph>,
          connectit::liutarjan_EUFA<Graph>, connectit::liutarjan_EUF<Graph>,
      });
  return 1.0;
}
}  // namespace gbbs

generate_symmetric_once_main(gbbs::Benchmark_runner, false);
