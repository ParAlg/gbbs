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

// Usage:
// numactl -i all ./StronglyConnectedComponents -beta 1.5 -rounds 2 -s -m twitter_J
// flags:
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -beta <value> : the base of the exponent to use (controls how quickly
//                     vertices are added)
//     -rounds : the number of times to run the algorithm
//     -stats : print the #sccs, and the #vertices in the largest scc

#include "StronglyConnectedComponents.h"

namespace gbbs {
namespace filter_based_scc {

template <class Seq>
inline size_t num_done(Seq& labels) {
  auto im_f = [&](size_t i) {
    return (size_t)(labels[i] != kUnfinished);
  };
  auto im = parlay::delayed_seq<size_t>(labels.size(), im_f);

  return parlay::reduce(im);
}

template <class Seq>
inline size_t num_scc(Seq& labels) {
  size_t n = labels.size();
  auto flags = sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  parallel_for(0, n, kDefaultGranularity, [&] (size_t i) {
    size_t label = labels[i];
    if ((label != kUnfinished) && !flags[label]) {
      flags[label] = 1;
    }
  });
  parlay::scan_inplace(flags);
  size_t n_scc = flags[n];
  std::cout << "n_scc = " << flags[n] << "\n";
  return n_scc;
}

template <class Seq>
inline void scc_stats(Seq& labels) {
  size_t n = labels.size();
  auto flags = sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    size_t label = labels[i];
    if (label != kUnfinished)
      flags[label]++;
  }
  size_t maxv = parlay::reduce_max(flags);
  std::cout << "Largest StronglyConnectedComponents has " << maxv << " vertices"
            << "\n";
  for (size_t i = 0; i < n; i++) {
    if (flags[i] == maxv) {
      std::cout << "max flag = " << i << "\n";
    }
  }
}


template <class Graph>
double StronglyConnectedComponents_runner(Graph& G, commandLine P) {
  double beta = P.getOptionDoubleValue("-beta", 1.1);
  std::cout << "### Application: StronglyConnectedComponents (Strongly Connected Components)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -beta = " << beta << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  assert(!P.getOption("-s"));
  timer scc_t;
  scc_t.start();
  auto labels = StronglyConnectedComponents(G, beta);
  double tt = scc_t.stop();
  if (P.getOption("-stats")) {
    num_scc(labels);
    scc_stats(labels);
  }

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}
}  // namespace filter_based_scc
}  // namespace gbbs

generate_asymmetric_main(gbbs::filter_based_scc::StronglyConnectedComponents_runner, false);
