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
void ProcessGraphUnweightedAverage(ClusteredGraph& CG, Sim lower_threshold, Sim max_weight, parlay::random& rnd) {
    std::cout << "Thresholds: " << lower_threshold << " and " << max_weight << std::endl;
  using W = typename ClusteredGraph::W;
  // Identify vertices with edges between [lower_threshold, max_weight)
  size_t n = CG.n;
  auto active = sequence<uintE>(n, 0);
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
