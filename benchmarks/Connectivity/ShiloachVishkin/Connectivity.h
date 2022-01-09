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
namespace shiloachvishkin_cc {

template <class Graph>
struct SVAlgorithm {
  Graph& GA;
  SVAlgorithm(Graph& GA) : GA(GA) {}
  sequence<parent> prev_parents;
  sequence<bool> flags;

  void initialize(sequence<parent>& P) {
    prev_parents = P;
    flags = sequence<bool>(P.size(), false);
  }

  template <SamplingOption sampling_option>
  void compute_components(sequence<parent>& parents,
                          parent frequent_comp = UINT_E_MAX) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    bool changed = true;
    size_t rounds = 0;
    std::cout << "# frequent_comp = " << frequent_comp << std::endl;

    /* generate candidates based on frequent_comp (if using sampling) */
    size_t candidates_size = n;
    sequence<uintE> unhooked;
    if
      constexpr(sampling_option != no_sampling) {
        auto all_vertices =
            parlay::delayed_seq<uintE>(n, [&](size_t i) { return i; });
        unhooked = parlay::filter(
            all_vertices, [&](uintE v) { return parents[v] != frequent_comp; });
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

    while (changed) {
      changed = false;
      rounds++;
      std::cout << "# round = " << rounds << std::endl;
      parallel_for(
          0, candidates.size(),
          [&](uintE i) {
            uintE u = candidates[i];
            auto map_f = [&](const uintE& _u, const uintE& _v, const W& wgh) {
              parent p_u = prev_parents[_u];
              parent p_v = prev_parents[_v];
              parent l = std::min(p_u, p_v);
              parent h = std::max(p_u, p_v);
              if (l != h && h == prev_parents[h]) {
                gbbs::write_min<parent>(&parents[h], l, std::less<parent>());
                if (!changed) {
                  changed = true;
                }
              }
            };
            GA.get_vertex(u).out_neighbors().map(map_f);
          },
          1);

      // compress
      parallel_for(0, n, [&](uintE u) {
        uintE pathlen = 1;
        while (parents[u] != parents[parents[u]]) {
          parents[u] = parents[parents[u]];
          pathlen++;
        }
        prev_parents[u] = parents[u];
        report_pathlen(pathlen);
      });
    }
    std::cout << "#rounds = " << rounds << std::endl;
  }

  template <bool reorder_updates, class Seq>
  void process_batch(sequence<parent>& parents, Seq& updates) {
    static_assert(reorder_updates == false);
    bool changed = true;

    size_t rounds = 0;
    while (changed) {
      rounds++;
      std::cout << "# running round = " << rounds << std::endl;
      changed = false;
      parallel_for(0, updates.size(), [&](size_t i) {
        parent u, v;
        UpdateType utype;
        std::tie(u, v, utype) = updates[i];
        if (utype == insertion_type) { /* update */
          parent p_u = prev_parents[u];
          parent p_v = prev_parents[v];
          parent l = std::min(p_u, p_v);
          parent h = std::max(p_u, p_v);
          if (l != h && h == prev_parents[h]) {
            gbbs::write_min<parent>(&parents[h], l, std::less<parent>());
            if (!changed) {
              changed = true;
            }
          }
        } /* ignore queries for now */
      });

      //      auto diff_map = parlay::delayed_seq<size_t>(parents.size(), [&]
      //      (size_t i) {
      //        return parents[i] != prev_parents[i];
      //      });

      // compress
      parallel_for(0, updates.size(), [&](size_t i) {
        uintE pathlen = 1;
        auto[u, v, utype] = updates[i];
        (void)utype;
        if (flags[u] == false &&
            gbbs::atomic_compare_and_swap(&flags[u], false, true)) {
          while (parents[u] != parents[parents[u]]) {
            parents[u] = parents[parents[u]];
            pathlen++;
          }
          prev_parents[u] = parents[u];
          report_pathlen(pathlen);
        }

        if (flags[v] == false &&
            gbbs::atomic_compare_and_swap(&flags[v], false, true)) {
          while (parents[v] != parents[parents[v]]) {
            parents[v] = parents[parents[v]];
            pathlen++;
          }
          prev_parents[v] = parents[v];
          report_pathlen(pathlen);
        }
      });

      // reset flags
      parallel_for(0, updates.size(), [&](size_t i) {
        auto[u, v, utype] = updates[i];
        (void)utype;
        if (flags[u]) {
          flags[u] = false;
        }
        if (flags[v]) {
          flags[v] = false;
        }
      });

      // compress (also performs queries implicitly on last round)
      //      parallel_for(0, updates.size(), [&] (uintE i) {
      //        uintE pathlen = 1;
      //        auto [u, v, utype] = updates[i];
      //        while (parents[u] != parents[parents[u]]) {
      //          parents[u] = parents[parents[u]];
      //          pathlen++;
      //        }
      //        report_pathlen(pathlen);
      //
      //        pathlen = 1;
      //        while (parents[v] != parents[parents[v]]) {
      //          parents[v] = parents[parents[v]];
      //          pathlen++;
      //        }
      //        prev_parents[u] = parents[u];
      //        prev_parents[v] = parents[v];
      //        report_pathlen(pathlen);
      //      });
    }
  }
};

}  // namespace shiloachvishkin_cc
}  // namespace gbbs
