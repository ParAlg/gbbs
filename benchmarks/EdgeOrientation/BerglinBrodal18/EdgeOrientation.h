#pragma once

#include <queue>
// #include <unordered_set>
// #include <vector>
#include <stdlib.h>
#include <math.h>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

namespace gbbs {

struct EdgeQueues {
  size_t n;  // number of vertices;
  size_t k;
  bool count_flips = false;

  //using level = std::unordered_set<uintE>;
  using Q = std::queue<uintE>;
  using edge_type = std::pair<uintE, uintE>;

  parlay::sequence<Q> L;

  EdgeQueues(size_t _n, size_t _k, bool _count_flips) : n(_n), k(_k), count_flips(_count_flips) {

    L = parlay::sequence<Q>(n);
  }

  bool insert_edge(edge_type e) {
    auto[u, v] = e;

    return true;
  }

  bool delete_edge(edge_type e) {
    auto[u, v] = e;


    return true;
  }


};

template <class Graph>
inline void RunEdgeOrientation(Graph& G, EdgeQueues& q) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  for (size_t i = 0; i < n; i++) {
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (u < v) {
        q.insert_edge({u, v});
      }
    };
    G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
  }

  std::cout << "Finished all insertions!" << std::endl;
  // layers.check_invariants();
  // std::cout << "Coreness estimate = " << layers.max_coreness() << std::endl;

  for (size_t i = 0; i < n; i++) {
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (u < v) {
        q.delete_edge({u, v});
      }
    };
    G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
  }

  std::cout << "Finished all deletions!" << std::endl;
  // layers.check_invariants();
  // std::cout << "Coreness estimate = " << layers.max_coreness() << std::endl;

  // std::cout << "Total level increases and decreases: " << layers.total_work
  //           << std::endl;
}

template <class W>
inline void RunEdgeOrientation(BatchDynamicEdges<W>& batch_edge_list, long batch_size,
  bool count_flips, EdgeQueues& q, size_t offset) {
  auto batch = batch_edge_list.edges;
  if (offset != 0) {
    for (size_t i = 0; i < offset; i++) {
      if (batch[i].insert) q.insert_edge({batch[i].from, batch[i].to});
      else q.delete_edge({batch[i].from, batch[i].to});
    }
  }
  
  for (size_t i = offset; i < batch.size(); i += batch_size) {
    timer t; t.start();
    for (size_t j = i; j < std::min(batch.size(), i + batch_size); j++) {
      if (batch[j].insert) q.insert_edge({batch[j].from, batch[j].to});
      else q.delete_edge({batch[j].from, batch[j].to});
    }

    double tt = t.stop();
    std::cout << "### Batch Running Time: " << tt << std::endl;
    // std::cout << "### Batch Num: " <<
    //   std::min(batch.size(), i + batch_size) - offset << std::endl;
    // std::cout << "### Coreness Estimate: " << layers.max_coreness() << std::endl;
  }
}

template <class Graph, class W>
inline void RunEdgeOrientation(Graph& G, BatchDynamicEdges<W> batch_edge_list,
  long batch_size, bool count_flips, size_t k) {
  uintE max_vertex = std::max( static_cast<uintE>(G.n), batch_edge_list.max_vertex);
  auto q = EdgeQueues(max_vertex, k, count_flips);
  // if (G.n > 0) RunLDS(G, layers);
  if (batch_edge_list.max_vertex > 0)
    RunEdgeOrientation(batch_edge_list, batch_size, count_flips, q, 0);
}

}  // namespace gbbs
