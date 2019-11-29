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
#include "benchmarks/Connectivity/common.h"

namespace shiloachvishkin_cc {

template <class Graph>
struct SVAlgorithm {
  Graph& GA;
  SVAlgorithm(Graph& GA) : GA(GA) {}

  void initialize(pbbs::sequence<parent>& P) {}

  template <SamplingOption sampling_option>
  void compute_components(pbbs::sequence<parent>& parents, parent frequent_comp = UINT_E_MAX) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    bool changed = true;
    size_t rounds = 0;

    pbbs::sequence<parent> clusters;
    if constexpr (sampling_option != no_sampling) {
      clusters = parents;
    }

    while (changed) {
      changed = false;
      rounds++;
      parallel_for(0, n, [&] (uintE u) {
        auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
          parent p_u = parents[u];
          parent p_v = parents[v];
          parent l = std::min(p_u, p_v);
          parent h = std::max(p_u, p_v);
          if (l != h  &&    h == parents[h]) {
            pbbs::write_min<parent>(&parents[h], l, std::less<parent>());
            if (!changed) { changed = true; }
          }
        };
        if constexpr (sampling_option != no_sampling) {
          if (clusters[u] != frequent_comp) {
            GA.get_vertex(u).mapOutNgh(u, map_f);
          }
        } else {
          GA.get_vertex(u).mapOutNgh(u, map_f);
        }
      }, 1);

      // compress
      parallel_for(0, n, [&] (uintE u) {
        while (parents[u] != parents[parents[u]]) {
          parents[u] = parents[parents[u]];
        }
      });
    }
    std::cout << "#rounds = " << rounds << std::endl;
  }

  template <class Seq>
  void process_batch(pbbs::sequence<parent>& parents, Seq& batch, size_t insert_to_query) {

    bool changed = true;
    while (changed) {
      changed = false;
      parallel_for(0, batch.size(), [&] (size_t i) {
        parent u, v;
        std::tie(u,v) = batch[i];
        if (i % insert_to_query != 0) { /* update */
          parent p_u = parents[u];
          parent p_v = parents[v];
          if (p_u < p_v && p_u == parents[p_u]) {
            pbbs::write_min<parent>(&parents[p_v], p_u, std::less<parent>());
            if (!changed) { changed = true; }
          }
        } else { /* query, ignore until the find in the next step*/
        }
      });

      // compress the edges in this batch
      parallel_for(0, batch.size(), [&] (uintE i) {
        uintE u, v;
        std::tie(u,v) = batch[i];
        while (parents[u] != parents[parents[u]]) {
          parents[u] = parents[parents[u]];
        }
        while (parents[v] != parents[parents[v]]) {
          parents[v] = parents[parents[v]];
        }
      });
    }
  }

};


}  // namespace shiloachvishkin_cc
