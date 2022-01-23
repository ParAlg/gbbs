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
#include "benchmarks/Connectivity/LabelPropagation/Connectivity.h"
#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"
#include "benchmarks/Connectivity/common.h"

#include "bench_utils.h"

namespace gbbs {
namespace connectit {

template <class Graph>
double t_bfs_cc(Graph& G, commandLine P, sequence<parent>& correct) {
  time(t, auto CC = bfs_cc::CC(G));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

template <class Graph>
void bfscc_nosample(Graph& G, int rounds, commandLine& P,
                    sequence<parent>& correct) {
  run_multiple(G, rounds, correct, "bfs_cc", P, t_bfs_cc<Graph>);
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
  run_tests(G, rounds, P, correct, connectit::bfscc_nosample<Graph>,
            {connectit::bfscc_nosample<Graph>});
  return 1.0;
}
}  // namespace gbbs

generate_symmetric_once_main(gbbs::Benchmark_runner, false);
