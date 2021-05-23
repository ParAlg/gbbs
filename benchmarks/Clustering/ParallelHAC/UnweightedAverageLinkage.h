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

#include <queue>
#include <unordered_set>
#include <vector>

#include "gbbs/gbbs.h"
#include "ClusteredGraph.h"

namespace gbbs {
namespace clustering {

struct vtx_status {
  vtx_status() {}

  bool is_red() { return red; }
  bool is_blue() { return !red; }

  bool red;
};

// At the end of this call, we will have performed merges and ensured that no
// edges exist with weights between [lower_threshold, ...).
template <class ClusteredGraph, class Sim>
void ProcessGraphUnweightedAverage(ClusteredGraph& CG, Sim lower_threshold, Sim max_weight, parlay::random& rnd,
    double eps = 0.05) {
    std::cout << "Thresholds: " << lower_threshold << " and " << max_weight << std::endl;
  using W = typename ClusteredGraph::W;
  // Identify vertices with edges between [lower_threshold, max_weight)
  size_t n = CG.n;
  auto active = sequence<uintE>(n, 0);
  auto colors = sequence<uint8_t>(n, 0);  // 1 is blue, 2 is red
  uint8_t kBlue = 1;
  uint8_t kRed = 2;
  parallel_for(0, n, [&] (size_t i) {
    auto pred_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      assert(wgh.get_weight() <= max_weight);
      return (wgh.get_weight() >= lower_threshold) && (wgh.get_weight() <= max_weight);
    };
    size_t ct = CG.clusters[i].countNeighbors(i, pred_f);
    if (ct > 0) active[i] = ct;
  });

  // use dense iterations for now.
  size_t n_active = pbbslib::reduce(parlay::delayed_seq<uintE>(n, [&] (size_t i) { return active[i] != 0; }));
  std::cout << "nactive = " << n_active << std::endl;

  double one_plus_eps = 1 + eps;

  while (n_active > 0) {

    parallel_for(0, n, [&] (size_t i) {
    if (active[i] > 0) {
      // Decide whether red or blue
      if (rnd.ith_rand(i) & 1) {
        colors[i] = kRed;
      } else {
        colors[i] = kBlue;
      }
    }});

    std::cout << "num_blue = " <<
      pbbslib::reduce(parlay::delayed_seq<uintE>(n, [&] (size_t i) { return colors[i] == kBlue; })) << std::endl;
    std::cout << "num_red = " <<
      pbbslib::reduce(parlay::delayed_seq<uintE>(n, [&] (size_t i) { return colors[i] == kRed; })) << std::endl;

    auto merge_target = sequence<uintE>(n, UINT_E_MAX);

    parallel_for(0, n, [&] (size_t i) {
    if (active[i] > 0 && colors[i] == kBlue) {
      // Select a random heavy-edge incident to the vertex.
      uintE active_i = active[i];
      size_t edge_idx = rnd.ith_rand(n + i) % active_i;
      // seq selection for now (todo: could run in par for high-degrees)
      size_t k = 0;
      uintE ngh_id = std::numeric_limits<uintE>::max();
      W weight;
      auto iter_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        if (wgh.get_weight() >= lower_threshold) {
          if (k == edge_idx) {
            ngh_id = v; weight = wgh;
          } else {
            k++;
          }
        }
      };
      CG.clusters[i].iterate(i, iter_f);
      assert(ngh_id != std::numeric_limits<uintE>::max());

      // Try to join the neighbor's cluster if neighbor is red.
      if (ngh_id == kRed) {
        assert(clusters[ngh_id].active);
        assert(clusters[i].active);
        uintE staleness = CG.clusters[ngh_id].staleness;
        uintE upper_bound = one_plus_eps * staleness;
        uintE ngh_cur_size = CG.clusters[ngh_id].cas_size;
        uintE our_size = CG.clusters[i].num_in_cluster;
        auto old_opt = pbbslib::fetch_and_add_threshold(
            &(CG.clusters[ngh_id].cas_size),
            ngh_cur_size,
            ngh_cur_size + our_size);
        if (old_opt.has_value()) {  // Success
          merge_target[i] = ngh_id;
        }
      }
    }});

    std::cout << "Num merges = " << pbbslib::reduce(parlay::delayed_seq<uintE>(n, [&] (size_t i) { return merge_target[i] != UINT_E_MAX; }))
      << std::endl;
    exit(0);



    rnd = rnd.next();
    // Reset colors.
    // Recompute active and update n_active.
  }
}





template <class ClusteredGraph, class W>
void ProcessEdgesUnweightedAverage(ClusteredGraph& CG, sequence<std::tuple<uintE, uintE, W>>&& edges,
                                   sequence<size_t>& Colors, parlay::random& rnd) {
  using edge = std::tuple<uintE, uintE, W>;
  auto V = sequence<vtx_status>::uninitialized(CG.n);
  auto edges_2 = sequence<edge>::uninitialized(edges.size());

  // Initialize Colors for vertices active in this round.
  auto UpdateColors = [&] () {
    auto latches = sequence<bool>(CG.n, false);
    auto gen_rand = [&] (const uintE& u) {
      // test and set
      if (!latches[u] && pbbslib::atomic_compare_and_swap(&latches[u], false, true)) {
        Colors[u] = rnd.ith_rand(u);  // TODO: make sure to update rnd outside.
      }
    };
    parallel_for(0, edges.size(), [&] (size_t i) {
      auto [u, v, wgh] = edges[i];
      gen_rand(u); gen_rand(v);
    });
  };
  UpdateColors();

  size_t phase = 0;
  constexpr size_t kPhaseMask = 63;

  auto GetColor = [&] (const uintE& u) -> bool {
    auto bit = phase & kPhaseMask;
    auto color = Colors[u];
    return color & (1 << bit);
  };

  std::cout << "Edges.size = " << edges.size() << std::endl;
}

}  // namespace clustering
}  // namespace gbbs
