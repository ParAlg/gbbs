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

#include "HeapBased.h"
#include "NNChainBased.h"
#include "HAC_configuration.h"

#include "AvgLinkageUtils/HeapBased.h"

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
W weightedavg_linkage(W w1, W w2) {
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


template <class Weights, class Dendrogram>
void WriteDendrogramToDisk(Weights& wgh, Dendrogram& dendrogram, const std::string& of) {
  ofstream out;
  out.open(of);
  size_t wrote = 0;
  for (size_t i=0; i<dendrogram.size(); i++) {
    if (dendrogram[i].first != i) {
      if (dendrogram[i].first != UINT_E_MAX) {
        out << i << " " << dendrogram[i].first << " " << Weights::AsString(dendrogram[i].second) << std::endl;
      }
//      else {
//        out << i << " " << i << " " << 0 << std::endl;  // cluster root
//      }
      wrote++;
    }
  }
  std::cout << "Wrote " << wrote << " parent-pointers. " << std::endl;
}

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


//  using W = typename Graph::weight_type;
//  auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
//    std::cout << u << " " << v << " " << wgh << std::endl;
//  };
//  G.get_vertex(4).out_neighbors().map(map_f, false);
//  std::cout << "0 edges" << std::endl;
//  G.get_vertex(0).out_neighbors().map(map_f, false);
//  exit(0);


  timer t; t.start();
  double tt;

if (heap_based) {
    if (linkage_opt == "weightedavg") {
      auto Wghs = WeightedAverageLinkage<Graph, DissimilarityClustering, ActualWeight>(G);
      auto dendrogram = heap_based::HAC(G, Wghs);
      tt = t.stop();
      std::cout << "### Running Time: " << tt << std::endl;
      auto of = P.getOptionValue("-of", "");
      if (of != "") {
        // write merges
        WriteDendrogramToDisk(Wghs, dendrogram, of);
        exit(0);
      }
    } else if (linkage_opt == "complete") {
      auto Wghs = MinLinkage<Graph, DissimilarityClustering, ActualWeight>(G);
      auto dendrogram = heap_based::HAC(G, Wghs);
      tt = t.stop();
      std::cout << "### Running Time: " << tt << std::endl;
      auto of = P.getOptionValue("-of", "");
      if (of != "") {
        // write merges
        WriteDendrogramToDisk(Wghs, dendrogram, of);
        exit(0);
      }
    } else if (linkage_opt == "single") {
      auto Wghs = MaxLinkage<Graph, DissimilarityClustering, ActualWeight>(G);
      auto dendrogram = heap_based::HAC(G, Wghs);
      tt = t.stop();
      std::cout << "### Running Time: " << tt << std::endl;
      auto of = P.getOptionValue("-of", "");
      if (of != "") {
        // write merges
        WriteDendrogramToDisk(Wghs, dendrogram, of);
        exit(0);
      }
    } else if (linkage_opt == "normalizedavg") {
      auto Wghs = NormAverageLinkage<Graph, SimilarityClustering, ActualWeight>(G);
      auto dendrogram = heap_based::HAC(G, Wghs);
      tt = t.stop();
      std::cout << "### Running Time: " << tt << std::endl;
      auto of = P.getOptionValue("-of", "");
      if (of != "") {
        // write merges
        WriteDendrogramToDisk(Wghs, dendrogram, of);
        exit(0);
      }
    } else if (linkage_opt == "avg") {
      std::cerr << "Average-linkage only supported for the similarity setting" << std::endl;
      exit(-1);
    } else {
      std::cerr << "Unknown linkage option: " << linkage_opt << std::endl;
      exit(-1);
    }
  } else {
    if (linkage_opt == "weightedavg") {
      auto Wghs = WeightedAverageLinkage<Graph, DissimilarityClustering, ActualWeight>(G);
      auto dendrogram = nn_chain::HAC(G, Wghs);
      tt = t.stop();
      std::cout << "### Running Time: " << tt << std::endl;
      auto of = P.getOptionValue("-of", "");
      if (of != "") {
        // write merges
        WriteDendrogramToDisk(Wghs, dendrogram, of);
        exit(0);
      }
    } else if (linkage_opt == "complete") {
      auto Wghs = MinLinkage<Graph, DissimilarityClustering, ActualWeight>(G);
      auto dendrogram = nn_chain::HAC(G, Wghs);
      tt = t.stop();
      std::cout << "### Running Time: " << tt << std::endl;
      auto of = P.getOptionValue("-of", "");
      if (of != "") {
        // write merges
        WriteDendrogramToDisk(Wghs, dendrogram, of);
        exit(0);
      }
    } else if (linkage_opt == "single") {
      auto Wghs = MaxLinkage<Graph, DissimilarityClustering, ActualWeight>(G);
      auto dendrogram = nn_chain::HAC(G, Wghs);
      tt = t.stop();
      std::cout << "### Running Time: " << tt << std::endl;
      auto of = P.getOptionValue("-of", "");
      if (of != "") {
        // write merges
        WriteDendrogramToDisk(Wghs, dendrogram, of);
        exit(0);
      }
    } else if (linkage_opt == "normalizedavg") {
      auto Wghs = NormAverageLinkage<Graph, DissimilarityClustering, ActualWeight>(G);
      auto dendrogram = nn_chain::HAC(G, Wghs);
      tt = t.stop();
      std::cout << "### Running Time: " << tt << std::endl;
      auto of = P.getOptionValue("-of", "");
      if (of != "") {
        // write merges
        WriteDendrogramToDisk(Wghs, dendrogram, of);
        exit(0);
      }
    } else if (linkage_opt == "appx-avg") {
      std::cout << "The approximate average linkage algorithm only supports the -heapbased option." << std::endl;
      exit(-1);
    } else {
      std::cerr << "Unknown linkage option: " << linkage_opt << std::endl;
      exit(-1);
    }
  }




  return tt;
}

}  // namespace gbbs

//  using W = typename Graph::weight_type;
//  auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
//    std::cout << u << " " << v << " " << wgh << std::endl;
//  };
//  G.get_vertex(19).out_neighbors().map(map_f, false);
//
//  std::cout << "32 edges" << std::endl;
//  G.get_vertex(32).out_neighbors().map(map_f, false);
//  exit(0);


generate_symmetric_weighted_main(gbbs::HAC_runner, false);
