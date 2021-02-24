#pragma once

#include <queue>
#include <set>
// #include <unordered_set>
// #include <vector>
#include <stdlib.h>
#include <math.h>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
// #include "benchmarks/KCore/JulienneDBS17/KCore.h"

namespace gbbs {

struct EdgeQueues {
  size_t n;  // number of vertices;
  size_t m = 0;
  size_t k;
  size_t num_flips = 0;

  using Q = std::set<uintE>; //can't use queue?
  using edge_type = std::pair<uintE, uintE>;

  parlay::sequence<Q> L;
  parlay::sequence<std::set<uintE>> A;
  parlay::sequence<size_t> out_degrees;
  size_t max_out_degree = 0;


  EdgeQueues(size_t _n, size_t _k, size_t max_degree) : n(_n), k(_k){
    L = parlay::sequence<Q>(n);
    A = parlay::sequence<std::stack<uintE>>(max_degree);
    out_degrees = parlay::sequence<size_t>(n, (size_t) 0);
  }

  inline void flip_k(){
    for(size_t i = 0; i < min(k, m); ++i){
      uintE v = pop(A[max_out_degree]);
      uintE u = remove_neighbor(v);
      insert_neighbor(u, v);
    }
    num_flips += k;
  }

  // qq can't be empty
  inline uintE pop(std::set<uintE> qq){
    if(qq.empty()){
      std::cout << "popped from empty queue!" << std::endl;
      return -1;
    }
    return *(qq.begin());
  }

  inline uintE pop(std::queue<uintE> qq){
     if(qq.empty()){
      std::cout << "popped from empty queue!" << std::endl;
      return -1;
    }
    return qq.pop();
  }

  inline uintE insert(std::set<uintE> qq, uintE u){
    qq.insert(u);
  }

  inline uintE insert(std::queue<uintE> qq, uintE u){
    qq.push(u);
  }

  inline uintE remove_neighbor(uintE v){
    uintE u = pop(L[v]);
    A[out_degrees[v]].erase(v);
    out_degrees[v] -= 1;
    A[out_degrees[v]].insert(v);
    while(A[max_out_degree].empty()){ max_out_degree--;}
    return u;
  }

  inline uintE remove_neighbor(uintE v, uintE u){
    L[v].erase(u);
    A[out_degrees[v]].erase(v);
    out_degrees[v] -= 1;
    A[out_degrees[v]].insert(v);
    while(A[max_out_degree].empty()){ max_out_degree--;}
    return u;
  }

  inline void insert_neighbor(uintE v, uintE u){
    insert(L[v], u);
    A[out_degrees[v]].erase(v);
    out_degrees[v] += 1;
    A[out_degrees[v]].insert(v);
    if(out_degrees[v] > max_out_degree) max_out_degree = out_degrees[v];
  }

  bool insert_edge(edge_type e) {
    auto[v, u] = e;
    m++;
    insert_neighbor(v, u);
    flip_k();
    return true;
  }

  bool delete_edge(edge_type e) {
    auto[v, u] = e;
    m--;
    remove_neighbor(v, u);
    flip_k();
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

  for (size_t i = 0; i < n; i++) {
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (u < v) {
        q.delete_edge({u, v});
      }
    };
    G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
  }

  std::cout << "Finished all deletions!" << std::endl;
  std::cout << "### Num Flips: " << q.num_flips << std::endl;
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

    if(count_flips){
      std::cout << "### Num Flips: " << q.num_flips << std::endl;
      std::cout << "### Max Out Degree: " << q.max_out_degree << std::endl;
    }
  }
}

template <class Graph, class W>
inline void RunEdgeOrientation(Graph& G, BatchDynamicEdges<W> batch_edge_list,
  long batch_size, bool count_flips, size_t k, size_t max_degree) {
  uintE max_vertex = std::max( static_cast<uintE>(G.n), batch_edge_list.max_vertex);
  auto q = EdgeQueues(max_vertex, k, max_degree);
  // if (G.n > 0) RunLDS(G, layers);
  if (batch_edge_list.max_vertex > 0)
    RunEdgeOrientation(batch_edge_list, batch_size, count_flips, q, 0);
}

}  // namespace gbbs
