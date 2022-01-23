#pragma once

#include "gbbs/gbbs.h"

#include "AvgLinkageUtils/HeapBased.h"
#include "HAC_configuration.h"
#include "HeapBased.h"
#include "NNChainBased.h"

namespace gbbs {

template <class Sim, class Graph>
sequence<std::pair<uintE, typename Graph::weight_type>> RunHAC(
    Graph& G, std::string linkage) {
  sequence<std::pair<uintE, typename Graph::weight_type>> dendrogram;
  if (linkage == "complete") {
    if
      constexpr(std::is_same<Sim, SimilarityClustering>()) {
        auto Wghs = MinLinkage<Graph, Sim, ActualWeight>(G);
        dendrogram = nn_chain::HAC(G, Wghs);
      }
    else {
      auto Wghs = MaxLinkage<Graph, DissimilarityClustering, ActualWeight>(G);
      dendrogram = nn_chain::HAC(G, Wghs);
    }
  } else if (linkage == "single") {
    if
      constexpr(std::is_same<Sim, SimilarityClustering>()) {
        auto Wghs = MaxLinkage<Graph, Sim, ActualWeight>(G);
        dendrogram = nn_chain::HAC(G, Wghs);
      }
    else {
      auto Wghs = MinLinkage<Graph, DissimilarityClustering, ActualWeight>(G);
      dendrogram = nn_chain::HAC(G, Wghs);
    }
  } else {
    std::cout << "Unknown linkage function" << std::endl;
    exit(0);
  }
  return dendrogram;
}

template <class Graph>
sequence<std::pair<uintE, typename Graph::weight_type>> HAC(
    Graph& G, std::string linkage, bool similarity = true) {
  std::cout << "linkage = " << linkage << " similarity = " << similarity
            << std::endl;
  if (similarity) return RunHAC<SimilarityClustering>(G, linkage);
  return RunHAC<DissimilarityClustering>(G, linkage);
}

}  // namespace gbbs
