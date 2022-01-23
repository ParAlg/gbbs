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

#include "benchmarks/SpanningForest/BFSSF/SpanningForest.h"
#include "benchmarks/SpanningForest/Framework/framework.h"
#include "benchmarks/SpanningForest/common.h"

#include "bench_utils.h"
#include "uf_utils.h"

namespace gbbs {
namespace connectit {
template <class Graph>
void unite_find_compress(Graph& G, int rounds, commandLine& P,
                         sequence<edge>& correct) {
  run_multiple_uf_alg<Graph, sample_bfs, unite, find_compress>(G, rounds,
                                                               correct, P);
}

template <class Graph>
void unite_find_naive(Graph& G, int rounds, commandLine& P,
                      sequence<edge>& correct) {
  run_multiple_uf_alg<Graph, sample_bfs, unite, find_naive>(G, rounds, correct,
                                                            P);
}

template <class Graph>
void unite_find_atomic_split(Graph& G, int rounds, commandLine& P,
                             sequence<edge>& correct) {
  run_multiple_uf_alg<Graph, sample_bfs, unite, find_atomic_split>(G, rounds,
                                                                   correct, P);
}

template <class Graph>
void unite_find_atomic_halve(Graph& G, int rounds, commandLine& P,
                             sequence<edge>& correct) {
  run_multiple_uf_alg<Graph, sample_bfs, unite, find_atomic_halve>(G, rounds,
                                                                   correct, P);
}
}

template <class Graph>
double Benchmark_runner(Graph& G, commandLine P) {
  int test_num = P.getOptionIntValue("-t", -1);
  int rounds = P.getOptionIntValue("-r", 5);

  auto correct = sequence<edge>();
  if (P.getOptionValue("-check")) {
    correct = bfs_sf::SpanningForestDet(G);
  }
  run_tests(G, rounds, P, correct, connectit::unite_find_naive<Graph>,
            {connectit::unite_find_compress<Graph>,
             connectit::unite_find_naive<Graph>,
             connectit::unite_find_atomic_split<Graph>,
             connectit::unite_find_atomic_halve<Graph>});
  return 1.0;
}
}  // namespace gbbs

generate_symmetric_once_main(gbbs::Benchmark_runner, false);
