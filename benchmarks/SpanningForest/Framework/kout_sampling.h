#pragma once

#include "benchmarks/Connectivity/connectit.h"

namespace gbbs {
// returns a component labeling
// Based on the implementation in gapbs/cc.c, Thanks to S. Beamer + M. Sutton
// for the well documented reference implementation of afforest.
template <class G>
struct KOutSamplingTemplate {
  G& GA;
  uint32_t neighbor_rounds;

  KOutSamplingTemplate(G& GA, commandLine& P) : GA(GA) {
    neighbor_rounds = P.getOptionLongValue("-sample_rounds", 2L);
  }

  void link(uintE u, uintE v, sequence<parent>& Parents,
            sequence<edge>& Edges) {
    parent p1 = Parents[u];
    parent p2 = Parents[v];
    while (p1 != p2) {
      parent high = p1 > p2 ? p1 : p2;
      parent low = p1 + (p2 - high);
      parent p_high = Parents[high];
      // Was already 'low' or succeeded in writing 'low'
      if (p_high == low) break;
      if (p_high == high &&
          gbbs::atomic_compare_and_swap(&Parents[high], high, low)) {
        Edges[high] =
            std::make_pair(u, v);  // One write per root, store (u,v) in high.
        break;
      }
      p1 = Parents[Parents[high]];
      p2 = Parents[low];
    }
  }

  auto initial_spanning_forest() {
    using W = typename G::weight_type;
    size_t n = GA.n;
    std::cout << "# neighbor_rounds = " << neighbor_rounds << std::endl;

    auto Parents = sequence<parent>(n, [&](size_t i) { return i; });
    auto Edges = sequence<edge>(n, empty_edge);

    parlay::random rnd;
    uintE granularity = 1024;
    for (uint32_t r = 0; r < neighbor_rounds; r++) {
      if (r == 0) {
        /* First round: sample a directed forest and compress instead of using
         * unite */
        parallel_for(0, n,
                     [&](size_t u) {
                       auto u_vtx = GA.get_vertex(u);
                       if (u_vtx.out_degree() > r) {
                         uintE ngh;
                         W wgh;
                         std::tie(ngh, wgh) =
                             u_vtx.out_neighbors().get_ith_neighbor(r);
                         link(u, ngh, Parents, Edges);
                       }
                     },
                     granularity);
      } else {
        parallel_for(0, n,
                     [&](size_t u) {
                       auto u_rnd = rnd.fork(u);
                       auto u_vtx = GA.get_vertex(u);
                       if (u_vtx.out_degree() > 1) {
                         uintE deg = u_vtx.out_degree() - 1;
                         uintE ngh_idx = 1 + (u_rnd.rand() % deg);
                         auto[ngh, wgh] =
                             u_vtx.out_neighbors().get_ith_neighbor(ngh_idx);
                         link(u, ngh, Parents, Edges);
                       }
                     },
                     granularity);
      }
      // compress nodes fully (turns out this is faster)
      parallel_for(0, n,
                   [&](size_t u) {
                     while (Parents[u] != Parents[Parents[u]]) {
                       Parents[u] = Parents[Parents[u]];
                     }
                   },
                   granularity);
      rnd = rnd.next();
    }
    return std::make_pair(Parents, Edges);
  }
};
}  // namespace gbbs
