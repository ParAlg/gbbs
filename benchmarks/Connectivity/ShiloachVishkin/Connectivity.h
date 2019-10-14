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

#include "ligra/ligra.h"

namespace shiloachvishkin_cc {

template <class Graph>
inline sequence<uintE> CC(Graph& G) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  auto parents = pbbs::sequence<uintE>(n, [&] (uintE i) { return i; });
  bool changed = true;
  size_t rounds = 0;
  while (changed) {
    changed = false;
    rounds++;
    parallel_for(0, n, [&] (uintE u) {
      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        uintE p_u = parents[u];
        uintE p_v = parents[v];
        if (p_u != p_v) {
          uintE larger = std::max(u,v);
          uintE smaller = std::min(u,v); // tricks require sign extension
          if (larger == parents[larger]) {
            if (!changed) {
              changed = true;
            }
            pbbs::write_min(&parents[larger], smaller, std::less<uintE>());
          }
        }
      };
      G.get_vertex(u).mapOutNgh(u, map_f);
    }, 1);
  }

  // compress
  parallel_for(0, n, [&] (uintE u) {
    while (parents[u] != parents[parents[u]]) {
      parents[u] = parents[parents[u]];
    }
  });
  std::cout << "# Ran: " << rounds << " many rounds" << std::endl;
  return parents;
}


template <class Graph>
struct SVAlgorithm {
  Graph& GA;
  SVAlgorithm(Graph& GA) : GA(GA) {}

  template <bool provides_frequent_comp>
  void compute_components(pbbs::sequence<uintE>& parents, uintE frequent_comp = UINT_E_MAX) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    bool changed = true;
    size_t rounds = 0;

    /* Can this algorithm make use of checking against. frequent_comp? */
    while (changed) {
      changed = false;
      rounds++;
      parallel_for(0, n, [&] (uintE u) {
        auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
          uintE p_u = parents[u];
          uintE p_v = parents[v];
          if (p_u != p_v) {
            uintE larger = std::max(u,v);
            uintE smaller = std::min(u,v); // tricks require sign extension
            if (larger == parents[larger]) {
              if (!changed) {
                changed = true;
              }
              pbbs::write_min(&parents[larger], smaller, std::less<uintE>());
            }
          }
        };
      }, 1);
    }
  }
};


}  // namespace shiloachvishkin_cc
