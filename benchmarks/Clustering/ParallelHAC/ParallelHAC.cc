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
#include "HACConfiguration.h"

namespace gbbs {

template <class Weights, class Dendrogram>
void WriteAvgDendrogramToDisk(Weights& wgh, Dendrogram& dendrogram,
                           const std::string& of) {
  ofstream out;
  out.open(of);
  size_t wrote = 0;
  out.precision(std::numeric_limits<double>::max_digits10);
  for (size_t i = 0; i < dendrogram.size(); i++) {
    if (dendrogram[i].first != i) {
      if (dendrogram[i].first != UINT_E_MAX) {
        out << i << " " << dendrogram[i].first << " "
            << dendrogram[i].second << std::endl;
      }
      wrote++;
    }
  }
  std::cout << "Wrote " << wrote << " parent-pointers." << std::endl;
}

template <class Graph>
double HAC_runner(Graph& G, commandLine P) {

  string linkage_opt = P.getOptionValue("-linkage", "complete");
  double epsilon = P.getOptionDoubleValue("-epsilon", 0.1);
  bool no_fine_grained_merge = P.getOptionValue("-no_fine_grained_merge");

  std::cout << "### Application: HAC" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: linkage = " << linkage_opt << " epsilon = " << epsilon << " fine_grained_merge = " << (!no_fine_grained_merge) << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t; t.start();
  double tt;

  auto Weights = ApproxAverageLinkage<Graph, SimilarityClustering, ActualWeight>(G);
  auto dendrogram = clustering::ParallelUPGMA(G, Weights, epsilon, !no_fine_grained_merge);
  auto of = P.getOptionValue("-of", "");
  if (of != "") {
    // write merges
    WriteAvgDendrogramToDisk(Weights, dendrogram, of);
    exit(0);
  }

  return tt;
}

}  // namespace gbbs

//generate_symmetric_weighted_main(gbbs::HAC_runner, false);
generate_symmetric_float_weighted_main(gbbs::HAC_runner);
