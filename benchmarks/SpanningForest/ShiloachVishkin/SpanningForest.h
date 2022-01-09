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

#include "benchmarks/Connectivity/connectit.h"
#include "benchmarks/SpanningForest/common.h"
#include "gbbs/gbbs.h"

namespace gbbs {
namespace shiloachvishkin_sf {

template <class Graph>
struct SVAlgorithm {
  Graph& GA;
  SVAlgorithm(Graph& GA) : GA(GA) {}

  void initialize(sequence<parent>& P, sequence<edge>& E) {}

  template <SamplingOption sampling_option>
  void compute_spanning_forest(sequence<parent>& Parents, sequence<edge>& Edges,
                               parent frequent_comp = UINT_E_MAX) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    bool changed = true;
    size_t rounds = 0;

    /* generate candidates based on frequent_comp (if using sampling) */
    size_t candidates_size = n;
    sequence<uintE> unhooked;
    if
      constexpr(sampling_option != no_sampling) {
        auto all_vertices =
            parlay::delayed_seq<uintE>(n, [&](size_t i) { return i; });
        unhooked = parlay::filter(
            all_vertices, [&](uintE v) { return Parents[v] != frequent_comp; });
        candidates_size = unhooked.size();
      }

    auto candidates =
        parlay::delayed_seq<uintE>(candidates_size, [&](size_t i) {
          if
            constexpr(sampling_option == no_sampling) { return i; }
          else {
            return unhooked[i];
          }
        });

    std::cout << "## Frequent comp = " << frequent_comp << std::endl;
    std::cout << "## Starting loop: candidates.size = " << candidates_size
              << std::endl;
    auto PrevParents = Parents;
    while (changed) {
      changed = false;
      rounds++;
      parallel_for(0, candidates.size(), 1, [&](uintE i) {
        uintE u = candidates[i];
        auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
          parent p_u = PrevParents[u];
          parent p_v = PrevParents[v];
          parent l = std::min(p_u, p_v);
          parent h = std::max(p_u, p_v);
          if (l != h && h == PrevParents[h]) {
            gbbs::write_min<parent>(&Parents[h], l, std::less<parent>());
            if (!changed) {
              changed = true;
            }
          }
        };
        GA.get_vertex(u).out_neighbors().map(map_f);

      });

      parallel_for(0, candidates.size(), 1, [&](uintE i) {
        uintE u = candidates[i];
        auto map_f_2 = [&](const uintE& u, const uintE& v, const W& wgh) {
          parent p_u = PrevParents[u];
          parent p_v = PrevParents[v];
          parent l = std::min(p_u, p_v);
          parent h = std::max(p_u, p_v);
          if (l != h && h == PrevParents[h]) {
            if (Parents[h] == l) {
              gbbs::atomic_compare_and_swap(&Edges[h], empty_edge,
                                            std::make_pair(u, v));
            }
          }
        };
        GA.get_vertex(u).out_neighbors().map(map_f_2);
      });

      // compress
      parallel_for(0, n, [&](uintE u) {
        uintE pathlen = 1;
        while (Parents[u] != Parents[Parents[u]]) {
          Parents[u] = Parents[Parents[u]];
          pathlen++;
        }
        PrevParents[u] = Parents[u];
      });
    }
    std::cout << "#rounds = " << rounds << std::endl;
  }
};

template <class Graph>
inline sequence<edge> SpanningForest(Graph& G) {
  size_t n = G.n;
  auto Parents =
      sequence<parent>::from_function(n, [&](size_t i) { return i; });
  auto Edges = sequence<edge>(n, empty_edge);

  auto alg = SVAlgorithm<Graph>(G);
  alg.template compute_spanning_forest<no_sampling>(Parents, Edges);

  return parlay::filter(Edges, [&](const edge& e) { return e != empty_edge; });
}

}  // namespace shiloachvishkin_sf
}  // namespace gbbs
