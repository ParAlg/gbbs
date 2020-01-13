#pragma once

#include "benchmarks/Connectivity/connectit.h"

// returns a component labeling
// Based on the implementation in gapbs/cc.c, Thanks to S. Beamer + M. Sutton
// for the well documented reference implementation of afforest.
template <class G>
struct KOutSamplingTemplate {
  G& GA;
  uint32_t neighbor_rounds;

  KOutSamplingTemplate(
      G& GA,
      commandLine& P) :
   GA(GA) {
    neighbor_rounds = P.getOptionLongValue("-sample_rounds", 2L);
   }

  void link(uintE u, uintE v, pbbs::sequence<parent>& Parents, pbbs::sequence<edge>& Edges) {
    parent p1 = Parents[u];
    parent p2 = Parents[v];
    while (p1 != p2) {
      parent high = p1 > p2 ? p1 : p2;
      parent low = p1 + (p2 - high);
      parent p_high = Parents[high];
      // Was already 'low' or succeeded in writing 'low'
      if (p_high == low) break;
      if (p_high == high && pbbs::atomic_compare_and_swap(&Parents[high], high, low)) {
        Edges[high] = std::make_pair(u,v); // One write per root, store (u,v) in high.
        break;
      }
      p1 = Parents[Parents[high]];
      p2 = Parents[low];
    }
  }

  auto initial_spanning_forest() {
    using W = typename G::weight_type;
    size_t n = GA.n;
    cout << "# neighbor_rounds = " << neighbor_rounds << endl;

    auto Parents = pbbs::sequence<parent>(n, [&] (size_t i) { return i; });
    auto Edges = pbbs::sequence<edge>(n, empty_edge);

    pbbs::random rnd;
    uintE granularity = 1024;
    for (uint32_t r=0; r<neighbor_rounds; r++) {
      if (r == 0) {
        /* First round: sample a directed forest and compress instead of using
         * unite */
        parallel_for(0, n, [&] (size_t u) {
          auto u_vtx = GA.get_vertex(u);
          if (u_vtx.getOutDegree() > r) {
            uintE ngh; W wgh;
            std::tie(ngh, wgh) = u_vtx.get_ith_out_neighbor(u, r);
            link(u, ngh, Parents, Edges);
          }
        }, granularity);
      } else {
        parallel_for(0, n, [&] (size_t u) {
          auto u_rnd = rnd.fork(u);
          auto u_vtx = GA.get_vertex(u);
          if (u_vtx.getOutDegree() > 0) {
            uintE deglb = (1 << std::max((int)pbbs::log2_up(u_vtx.getOutDegree()), (int)1) - 1) - 1;
            uintE ngh_idx;
            if (deglb == 0) {
              ngh_idx = 0;
            } else {
              ngh_idx = u_rnd.rand() & deglb;
            }
            uintE ngh; W wgh;
            std::tie(ngh, wgh) = u_vtx.get_ith_out_neighbor(u, ngh_idx);
            link(u, ngh, Parents, Edges);
          }
        }, granularity);
      }
      // compress nodes fully (turns out this is faster)
      parallel_for(0, n, [&] (size_t u) {
        while (Parents[u] != Parents[Parents[u]]) {
          Parents[u] = Parents[Parents[u]];
        }
      }, granularity);
      rnd = rnd.next();
    }
    return std::make_pair(Parents, Edges);
   }

};
