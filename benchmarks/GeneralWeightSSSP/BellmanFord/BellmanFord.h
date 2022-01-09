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

namespace gbbs {

template <class W, class Distance>
struct BF_F {
  Distance* SP;
  intE* Visited;
  BF_F(Distance* _SP, intE* _Visited) : SP(_SP), Visited(_Visited) {}
  inline bool update(const uintE& s, const uintE& d, const W& edgeLen) {
    Distance newDist;
    if
      constexpr(std::is_same<W, gbbs::empty>()) { newDist = SP[s] + 1; }
    else {
      newDist = SP[s] + edgeLen;
    }
    if (SP[d] > newDist) {
      SP[d] = newDist;
      if (Visited[d] == 0) {
        Visited[d] = 1;
        return 1;
      }
    }
    return 0;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& edgeLen) {
    Distance newDist;
    if
      constexpr(std::is_same<W, gbbs::empty>()) { newDist = SP[s] + 1; }
    else {
      newDist = SP[s] + edgeLen;
    }
    return (gbbs::write_min(&SP[d], newDist) &&
            gbbs::atomic_compare_and_swap(&Visited[d], 0, 1));
  }
  inline bool cond(uintE d) { return cond_true(d); }
};

// reset visited vertices
struct BF_Vertex_F {
  intE* Visited;
  BF_Vertex_F(intE* _Visited) : Visited(_Visited) {}
  inline bool operator()(uintE i) {
    Visited[i] = 0;
    return 1;
  }
};

template <class Graph>
auto BellmanFord(Graph& G, uintE start) {
  using W = typename Graph::weight_type;
  using Distance =
      typename std::conditional<std::is_same<W, gbbs::empty>::value, uintE,
                                W>::type;

  size_t n = G.n;
  auto Visited = sequence<int>(n, 0);
  auto SP = sequence<Distance>(n, std::numeric_limits<Distance>::max());
  SP[start] = 0;

  vertexSubset Frontier(n, start);
  size_t round = 0;
  while (!Frontier.isEmpty()) {
    // Check for a negative weight cycle
    if (round == n) {
      std::cout << " Found negative weight cycle." << std::endl;
      break;
    }
    auto em_f = BF_F<W, Distance>(SP.begin(), Visited.begin());
    auto output =
        edgeMap(G, Frontier, em_f, G.m / 10, sparse_blocked | dense_forward);
    vertexMap(output, BF_Vertex_F(Visited.begin()));
    std::cout << output.size() << "\n";
    Frontier = std::move(output);
    round++;
  }
  auto dist_im_f = [&](size_t i) {
    return (SP[i] == (std::numeric_limits<Distance>::max())) ? 0 : SP[i];
  };
  auto dist_im = parlay::delayed_seq<Distance>(n, dist_im_f);
  std::cout << "max dist = " << parlay::reduce_max(dist_im) << "\n";
  std::cout << "n rounds = " << round << "\n";
  return SP;
}

}  // namespace gbbs
