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
// numactl -i all ./MaximalIndependentSet -stats -rounds 4 -s
// com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc
//     -specfor : run the speculative_for based algorithm from pbbs

#include "MaximalIndependentSet.h"

namespace gbbs {

template <class Graph>
double MaximalIndependentSet_runner(Graph& G, commandLine P) {
  bool spec_for = P.getOption("-specfor");
  std::cout << "### Application: MaximalIndependentSet" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -specfor (deterministic reservations) = "
            << spec_for << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  assert(P.getOption("-s"));
  double tt = 0.0;

  // Code below looks duplicated; this is because the return types of specfor
  // and rootset are different
  if (spec_for) {
    timer t;
    t.start();
    auto MaximalIndependentSet =
        MaximalIndependentSet_spec_for::MaximalIndependentSet(G);
    // in spec_for, MaximalIndependentSet[i] == 1 indicates that i was chosen
    tt = t.stop();
    auto size_f = [&](size_t i) { return (MaximalIndependentSet[i] == 1); };
    auto size_imap = parlay::delayed_seq<size_t>(G.n, size_f);
    if (P.getOptionValue("-stats")) {
      std::cout << "MaximalIndependentSet size: " << parlay::reduce(size_imap)
                << "\n";
    }
    if (P.getOptionValue("-verify")) {
      verify_MaximalIndependentSet(G, size_imap);
    }
  } else {
    timer t;
    t.start();
    auto MaximalIndependentSet =
        MaximalIndependentSet_rootset::MaximalIndependentSet(G);
    tt = t.stop();
    auto size_f = [&](size_t i) { return MaximalIndependentSet[i]; };
    auto size_imap = parlay::delayed_seq<size_t>(G.n, size_f);
    if (P.getOptionValue("-stats")) {
      std::cout << "MaximalIndependentSet size: " << parlay::reduce(size_imap)
                << "\n";
    }
    if (P.getOptionValue("-verify")) {
      verify_MaximalIndependentSet(G, size_imap);
    }
  }

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

generate_main(gbbs::MaximalIndependentSet_runner, false);
