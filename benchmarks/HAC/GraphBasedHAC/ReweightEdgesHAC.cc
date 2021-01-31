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
// numactl -i all ./HAC -s -m -rounds 3 twitter_wgh_SJ
// flags:
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#include "HAC.h"
#include "HAC_configuration.h"

namespace gbbs {

template <class W>
W min_linkage(W w1, W w2) {
  return std::min(w1, w2);
}

template <class W>
W max_linkage(W w1, W w2) {
  return std::max(w1, w2);
}

template <class W>
W weighted_avg_linkage(W w1, W w2) {
  return (w1 + w2) / static_cast<W>(2);
}

//template <class Graph, class GetWeight = EmptyToLogW>
//struct DissimilarityMinLinkage : GetWeight::template GetWeight<Graph> {
//
//  using base = typename GetWeight::template GetWeight<Graph>;
//  using weight_type = typename base::weight_type;
//
//  using base::base;
//
//  static constexpr bool similarity_clustering = false;
//
//  // Used to specify whether we are doing similarity of dissimilarity
//  // clustering. Similarity means taking max (heavier weights are more similar)
//  // and dissimilarity means taking min (smaller edges are "closer")
//  static weight_type augmented_combine(const weight_type& lhs, const weight_type& rhs) {
//    return std::min(lhs, rhs);   // similarity
//  }
//
//  static weight_type id() {
//    return std::numeric_limits<weight_type>::max();
//  }
//
//  // The linkage function.
//  static weight_type linkage(const weight_type& lhs, const weight_type& rhs) {
//    return std::min(lhs, rhs);
//  }
//};


template <class Graph>
double HAC_runner(Graph& G, commandLine P) {
  std::cout << "### Application: HAC" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: " << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  // using W = typename Graph::weight_type;
  // auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
  //   std::cout << u << " " << v << std::endl;
  // };
  // G.get_vertex(0).out_neighbors().map(map_f, false);

  // std::cout << "1 edges" << std::endl;
  // G.get_vertex(1).out_neighbors().map(map_f, false);

  timer t; t.start();

  // auto W = MinLinkage<Graph>(G);
  auto W = MaxLinkage<Graph>(G);
  //auto W = WeightedAverageLinkage<Graph>(G);
  // auto W = AverageLinkage<Graph>(G);
  auto ml = [] (float w1, float w2) { return std::min(w1, w2); };
  greedy_exact::HAC(G, W, ml);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

generate_symmetric_main(gbbs::HAC_runner, false);
