#include "HAC_configuration.h"
#include "HeapBased.h"
#include "NNChainBased.h"
#include "AvgLinkageUtils/HeapBased.h"


template <class Sim, class Graph>
sequence<std::pair<uintE, typename Weights::weight_type>> HAC(Graph& G, std::string linkage) {
  if (linkage == "complete") {
    if constexpr (std::is_same<Sim, Sim>()) {
      auto Wghs = MinLinkage<Graph, Sim, ActualWeight>(G);
      auto dendrogram = heap_based::HAC(G, Wghs);
    } else {
      auto Wghs = MinLinkage<Graph, DissimilarityClustering, ActualWeight>(G);
      auto dendrogram = heap_based::HAC(G, Wghs);
    }
  } else if (linkage == "single") {
    if constexpr (std::is_same<Sim, Sim>()) {
      auto Wghs = MaxLinkage<Graph, Sim, ActualWeight>(G);
      auto dendrogram = heap_based::HAC(G, Wghs);
    } else {
      auto Wghs = MaxLinkage<Graph, DissimilarityClustering, ActualWeight>(G);
      auto dendrogram = heap_based::HAC(G, Wghs);
    }
  } else if (linkage == "weightedavg") {
      auto Wghs =
          WeightedAverageLinkage<Graph, Sim, ActualWeight>(G);
      auto dendrogram = heap_based::HAC(G, Wghs);
  } else if (linkage == "avg") {
      auto Wghs =
          ApproxAverageLinkage<Graph, Sim, ActualWeight>(G);
      double epsilon = 0.1;
      auto dendrogram = approx_average_linkage::HAC(G, Wghs, epsilon);
  } else {
    std::cout << "Unknown linkage function" << std::endl;
    exit(0);
  }
}

template <class Graph>
void HAC(Graph& G, std::string linkage, bool similarity=true) {
  if (similarity) return RunHAC<SimilarityClustering>(G, linkage);
  return RunHAC<DissimilarityClustering>(G, linkage);
}
