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

#include "benchmarks/Connectivity/common.h"
#include "gbbs/gbbs.h"

namespace gbbs {
namespace bfs_cc {

template <class W>
struct BFS_ComponentLabel_F {
  parent* Parents;
  uintE src;
  BFS_ComponentLabel_F(parent* _Parents, uintE src)
      : Parents(_Parents), src(src) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (Parents[d] != src) {
      Parents[d] = src;
      return 1;
    } else {
      return 0;
    }
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) {
    return (gbbs::atomic_compare_and_swap(&Parents[d],
                                          static_cast<parent>(UINT_E_MAX),
                                          static_cast<parent>(src)));
  }
  inline bool cond(const uintE& d) { return (Parents[d] == UINT_E_MAX); }
};

/* Returns a mapping from either i --> i, if i is not reached by the BFS, or
 * i --> src, if i is reachable from src in the BFS */
template <class Graph>
void BFS_ComponentLabel(Graph& G, uintE src, sequence<parent>& parents) {
  using W = typename Graph::weight_type;
  if (G.get_vertex(src).out_degree() > 0) {
    vertexSubset Frontier(G.n, src);
    size_t reachable = 0;
    size_t rounds = 0;
    parents[src] = src;
    while (!Frontier.isEmpty()) {
      reachable += Frontier.size();
      vertexSubset output =
          edgeMap(G, Frontier, BFS_ComponentLabel_F<W>(parents.begin(), src),
                  -1, sparse_blocked | dense_parallel);
      if (output.size() > 0) {
        std::cout << "output.size = " << output.size() << std::endl;
      }
      Frontier = std::move(output);
      rounds++;
    }
  }
}

template <class Graph>
inline sequence<parent> CC(Graph& G) {
  size_t n = G.n;
  auto parents = sequence<parent>(n, UINT_E_MAX);
  for (size_t i = 0; i < n; i++) {
    if (parents[i] == UINT_E_MAX) {
      BFS_ComponentLabel(G, i, parents);
    }
  }
  return parents;
}

}  // namespace bfs_cc
}  // namespace gbbs
