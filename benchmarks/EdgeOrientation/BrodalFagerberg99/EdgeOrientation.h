#pragma once

#include <queue>
#include <set>
#include <stack>
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
  size_t Delta;
  size_t num_flips = 0;

  using Q = std::set<uintE>; //can't use queue?
  using edge_type = std::pair<uintE, uintE>;

  parlay::sequence<Q> L;
  parlay::sequence<size_t> out_degrees;


  EdgeQueues(size_t _n, size_t _Delta, size_t max_degree) : n(_n), Delta(_Delta){
    L = parlay::sequence<Q>(n);
    out_degrees = parlay::sequence<size_t>(n, (size_t) 0);
  }
  // qq can't be empty
  inline uintE pop(std::set<uintE> &qq){
    if(qq.empty()){
      std::cout << "popped from empty queue!" << std::endl;
      return -1;
    }
    uintE u = *(qq.begin());
    qq.erase(u);
    return u;
  }

  inline uintE pop(std::queue<uintE> &qq){
     if(qq.empty()){
      std::cout << "popped from empty queue!" << std::endl;
      return -1;
    }
    uintE u = qq.front();
    qq.pop();
    return u;
  }

  inline uintE pop(std::stack<uintE> &qq){
     if(qq.empty()){
      std::cout << "popped from empty stack!" << std::endl;
      return -1;
    }
    uintE u = qq.top();
    qq.pop();
    return u;
  }

  inline void insert(std::set<uintE> &qq, uintE u){
    qq.insert(u);
  }

  inline void insert(std::queue<uintE> &qq, uintE u){
    qq.push(u);
  }

  inline void insert(std::stack<uintE> &qq, uintE u){
    qq.push(u);
  }


  inline bool remove_neighbor(uintE v, uintE u){
    if(L[v].find(u)!=L[v].end()){
      L[v].erase(u);
      out_degrees[v] -= 1;
      return true;
    }
    return false;
  }

  // add u to be v's neighbor
  inline void insert_neighbor(uintE v, uintE u){
    insert(L[v], u);
    // A[out_degrees[v]].erase(v);
    out_degrees[v] += 1;
    // A[out_degrees[v]].insert(v);
    // if(out_degrees[v] > max_out_degree) max_out_degree = out_degrees[v];
    if(out_degrees[v] == Delta+1){
      std::stack<uintE> S = std::stack<uintE>();
      insert(S, v);
      size_t count = 0;
      while(!S.empty() && count < n){
        uintE w = pop(S);
        if(out_degrees[w] <= Delta) continue;
        count++;
        num_flips += L[w].size();
        for (auto it=L[w].begin(); it!=L[w].end(); ++it){
          uintE x = *it;
          insert(L[x], w);
           out_degrees[x] += 1;
          if(out_degrees[x] == Delta+1) insert(S, x);
        }
        L[w].clear();
        out_degrees[w] = 0;
      }
      if(count == n){
        std::cout << "while loop reached n" << std::endl;
        exit(1);
      }
    }
  }

  bool insert_edge(edge_type e) {
    auto[v, u] = e;
    m++;
    insert_neighbor(v, u);
    return true;
  }

  bool delete_edge(edge_type e) {
    auto[v, u] = e;
    m--;
    remove_neighbor(v, u);
    remove_neighbor(u, v);
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

  //TODO: add erase
  // for (size_t i = 0; i < n; i++) {
  //   auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
  //     if (u < v) {
  //       q.delete_edge({u, v});
  //     }
  //   };
  //   G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
  // }

  // std::cout << "Finished all deletions!" << std::endl;
  std::cout << "### Num Flips: " << q.num_flips << std::endl;
}

template <class W>
inline void RunEdgeOrientation(BatchDynamicEdges<W>& batch_edge_list, long batch_size,
  bool count_flips, EdgeQueues& q, size_t offset) {
  auto batch = batch_edge_list.edges;
  if (offset != 0) {
    for (size_t i = 0; i < offset; i++) {
      if (batch[i].insert) q.insert_edge({batch[i].from, batch[i].to});
      // else q.delete_edge({batch[i].fro .m, batch[i].to}); //TODO: add erase
    }
  }
  int batch_num = 0;
  for (size_t i = offset; i < batch.size(); i += batch_size) {
    std::cout << "### Batch Num: " << batch_num++ << std::endl;
    timer t; t.start();
    for (size_t j = i; j < std::min(batch.size(), i + batch_size); j++) {
      if (batch[j].insert) q.insert_edge({batch[j].from, batch[j].to});
      // else q.delete_edge({batch[j].from, batch[j].to});  //TODO: add erase
    }

    double tt = t.stop();
    std::cout << "### Batch Running Time: " << tt << std::endl;

    if(count_flips){
      std::cout << "### Num Flips: " << q.num_flips << std::endl;
      std::cout << "### Max Out Degree: " << *parlay::max_element(q.out_degrees, std::less<size_t>()) << std::endl;
    }
  }
}

template <class Graph, class W>
inline void RunEdgeOrientation(Graph& G, BatchDynamicEdges<W> batch_edge_list,
  long batch_size, bool count_flips, size_t k, size_t max_degree) {
  uintE max_vertex = std::max( static_cast<uintE>(G.n), batch_edge_list.max_vertex);
  std::cout << "Delta = " << k << std::endl;
  auto q = EdgeQueues(max_vertex, k, max_degree);
  // if (G.n > 0) RunLDS(G, layers);
  if (batch_edge_list.max_vertex > 0)
    RunEdgeOrientation(batch_edge_list, batch_size, count_flips, q, 0);
}

}  // namespace gbbs
