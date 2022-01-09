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

namespace gbbs {
namespace connectit {
template <class Graph, SamplingOption sampling_option,
          template <class G> class Algorithm, AlgorithmType algorithm_type>
bool run_multiple_sample_only_alg(Graph& G, size_t rounds,
                                  sequence<parent>& correct, commandLine& P,
                                  std::string name) {
  auto test = [&](Graph& graph, commandLine params,
                  sequence<parent>& correct_cc) {
    timer tt;
    tt.start();
    auto CC =
        run_sample_only_alg<Graph, sampling_option, Algorithm, algorithm_type>(
            graph, params);
    double t = tt.stop();
    if (params.getOptionValue("-check")) {
      cc_check(correct_cc, CC);
    }
    return t;
  };
  auto test_name = name + "; " + sampling_to_string(sampling_option);
  return run_multiple(G, rounds, correct, test_name, P, test);
}

template <class Graph>
void shiloachvishkin_nosample(Graph& G, int rounds, commandLine& P,
                              sequence<parent>& correct) {
  run_multiple_sample_only_alg<Graph, no_sampling,
                               shiloachvishkin_cc::SVAlgorithm,
                               shiloach_vishkin_type>(G, rounds, correct, P,
                                                      "shiloach_vishkin");
}

template <class Graph>
void shiloachvishkin_kout(Graph& G, int rounds, commandLine& P,
                          sequence<parent>& correct) {
  run_multiple_sample_only_alg<Graph, sample_kout,
                               shiloachvishkin_cc::SVAlgorithm,
                               shiloach_vishkin_type>(G, rounds, correct, P,
                                                      "shiloach_vishkin");
}

template <class Graph>
void shiloachvishkin_bfs(Graph& G, int rounds, commandLine& P,
                         sequence<parent>& correct) {
  run_multiple_sample_only_alg<Graph, sample_bfs,
                               shiloachvishkin_cc::SVAlgorithm,
                               shiloach_vishkin_type>(G, rounds, correct, P,
                                                      "shiloach_vishkin");
}

template <class Graph>
void shiloachvishkin_ldd(Graph& G, int rounds, commandLine& P,
                         sequence<parent>& correct) {
  run_multiple_sample_only_alg<Graph, sample_ldd,
                               shiloachvishkin_cc::SVAlgorithm,
                               shiloach_vishkin_type>(G, rounds, correct, P,
                                                      "shiloach_vishkin");
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
  run_tests(G, rounds, P, correct, connectit::shiloachvishkin_nosample<Graph>,
            {connectit::shiloachvishkin_nosample<Graph>,
             connectit::shiloachvishkin_kout<Graph>,
             connectit::shiloachvishkin_bfs<Graph>,
             connectit::shiloachvishkin_ldd<Graph>});
  return 1.0;
}
}  // namespace gbbs

generate_symmetric_once_main(gbbs::Benchmark_runner, false);
