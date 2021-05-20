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

#include "ParallelHAC.h"

namespace gbbs {


struct ActualWeight {
  struct data {};
  template <class Graph, class WeightType=float>
  struct GetWeight {
    using weight_type = typename Graph::weight_type;
    using underlying_weight_type = typename Graph::weight_type;
    Graph& G;

    GetWeight(Graph& G) : G(G) {}

    weight_type id() { return 0; }

    static weight_type linkage (const weight_type& lhs, const weight_type& rhs) {
      return lhs + rhs / static_cast<weight_type>(2);
    }

    // Convert an underlying weight to an initial edge weight for this edge.
    weight_type get_weight(const uintE& u, const uintE& v, const underlying_weight_type& wgh) const {
      return wgh;
    }

  };
};


template <class Graph>
double HAC_runner(Graph& G, commandLine P) {

  bool heap_based = P.getOptionValue("-heapbased");
  string linkage_opt = P.getOptionValue("-linkage", "complete");

  std::cout << "### Application: HAC" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: heap-based = " << heap_based << " linkage = " << linkage_opt << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t; t.start();
  double tt;

  auto Weights = ActualWeight::template GetWeight<Graph>(G);
  clustering::ParallelUPGMA(G, Weights);

  return tt;
}

}  // namespace gbbs

generate_symmetric_weighted_main(gbbs::HAC_runner, false);
