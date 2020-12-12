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

#pragma once

#include "gbbs/gbbs.h"
#include "benchmarks/SpanningForest/common.h"

namespace gbbs {
namespace bfs_sf {

template <class W>
struct BFS_SpanningForest_F {
  parent* Parents;
  BFS_SpanningForest_F(parent* _Parents) : Parents(_Parents) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    Parents[d] = s;
    return 1;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) {
    return (pbbs::atomic_compare_and_swap(&Parents[d], static_cast<parent>(UINT_E_MAX), static_cast<parent>(s)));
  }
  inline bool cond(const uintE& d) { return (Parents[d] == UINT_E_MAX); }
};

/* nondeterministic version */
template <class Graph>
void BFS_SpanningForest(Graph& G, uintE src, pbbs::sequence<parent>& parents) {
  using W = typename Graph::weight_type;
  vertexSubset Frontier(G.n, src);
  size_t reachable = 0; size_t rounds = 0;
  while (!Frontier.isEmpty()) {
    reachable += Frontier.size();
    vertexSubset output =
        edgeMap(G, Frontier, BFS_SpanningForest_F<W>(parents.begin()), -1, sparse_blocked | dense_parallel);
    Frontier.del();
    Frontier = output;
    rounds++;
  }
  Frontier.del();
}


template <class Graph>
inline pbbs::sequence<edge> SpanningForest(Graph& G) {
  size_t n = G.n;
  auto parents = pbbs::sequence<parent>(n, UINT_E_MAX);
  for (size_t i=0; i<n; i++) {
    if (parents[i] == UINT_E_MAX) {
      parents[i] = i;
      BFS_SpanningForest(G, i, parents);
    }
  }
  return spanning_forest::parents_to_edges(parents);
}

/* deterministic version */
template <class W>
struct BFS_SpanningForest_Det_F {
  pbbs::sequence<parent>& Parents;
  pbbs::sequence<bool>& visited;
  BFS_SpanningForest_Det_F(pbbs::sequence<parent>& _Parents, pbbs::sequence<bool>& visited) : Parents(_Parents), visited(visited) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (s < Parents[d]) {
      Parents[d] = s;
    }
    return false;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) {
    pbbs::write_min<parent>(&Parents[d], static_cast<parent>(s), std::less<parent>());
    return false;
  }
  inline bool cond(const uintE& d) { return (!visited[d]); }
};

template <class W>
struct BFS_SpanningForest_Det_F_2 {
  pbbs::sequence<parent>& Parents;
  pbbs::sequence<bool>& visited;
  BFS_SpanningForest_Det_F_2(pbbs::sequence<parent>& _Parents, pbbs::sequence<bool>& visited) : Parents(_Parents), visited(visited) {}
  inline bool update(const uintE& s, const uintE& d, const W& wgh) {  // Update
    if (Parents[d] == s) {
      visited[d] = true;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {  // Atomic version of Update
    if (Parents[d] == s) {
      visited[d] = true;
      return true;
    }
    return false;
  }
  // Cond function checks if vertex has been visited yet
  inline bool cond(uintE d) { return !visited[d]; }
};

template <class Graph>
void BFS_SpanningForest_Det(Graph& G, uintE src, pbbs::sequence<parent>& parents, pbbs::sequence<bool>& visited) {
  using W = typename Graph::weight_type;
  vertexSubset Frontier(G.n, src);
  size_t reachable = 0; size_t rounds = 0;
  while (!Frontier.isEmpty()) {
    reachable += Frontier.size();
    edgeMap(G, Frontier, BFS_SpanningForest_Det_F<W>(parents, visited), -1, sparse_blocked | dense_parallel);
    auto output = edgeMap(G, Frontier, BFS_SpanningForest_Det_F_2<W>(parents, visited), -1, sparse_blocked | dense_parallel);
    Frontier.del();
    Frontier = output;
    rounds++;
  }
  Frontier.del();
}

template <class Graph>
inline pbbs::sequence<edge> SpanningForestDet(Graph& G) {
  size_t n = G.n;
  auto parents = pbbs::sequence<parent>(n, UINT_E_MAX);
  auto visited = pbbs::sequence<bool>(G.n, false);
  for (size_t i=0; i<n; i++) {
    if (parents[i] == UINT_E_MAX) {
      parents[i] = i;
      visited[i] = true;
      BFS_SpanningForest_Det(G, i, parents, visited);
    }
  }
  return spanning_forest::parents_to_edges(parents);
}

}  // namespace bfs_sf
}  // namespace gbbs
