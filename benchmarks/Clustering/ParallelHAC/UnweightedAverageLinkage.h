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
template <class Weights, class ClusteredGraph, class Sim>
void ProcessGraphUnweightedAverage(ClusteredGraph& CG, Sim lower_threshold, Sim max_weight, parlay::random& rnd,
    double eps = 0.5) {
  std::cout << "Thresholds: " << lower_threshold << " and " << max_weight << std::endl;
  using W = typename ClusteredGraph::W;

  // Identify vertices with edges between [lower_threshold, max_weight)
  size_t n = CG.n;
  std::cout << "n = " << n << std::endl;
  auto active = sequence<uintE>(n, 0);
  auto colors = sequence<uint8_t>(n, 0);  // 1 is blue, 2 is red
  uint8_t kBlue = 1;
  uint8_t kRed = 2;

  auto get_num_active_edges = [&] (uintE u) -> uintE {
    auto pred_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      Sim actual_weight = Weights::get_weight(wgh, u, v, CG);
      if (actual_weight > max_weight) {
        std::cout << "Actual weight = " << actual_weight << " max_weight = " << max_weight << std::endl;
        std::cout << "wgh = " << wgh << std::endl;
        std::cout << "u = " << u << " v = " << v << std::endl;
        std::cout << "u size = " << CG.clusters[u].cluster_size() << std::endl;
        std::cout << "v size = " << CG.clusters[v].cluster_size() << std::endl;
      }
      assert(actual_weight <= max_weight);
      return (actual_weight >= lower_threshold);
    };
    return CG.clusters[u].countNeighbors(u, pred_f);
  };

  parallel_for(0, n, [&] (size_t i) {
    if (CG.clusters[i].active) {
      uintE ct = get_num_active_edges(i);
      if (ct > 0) active[i] = ct;
    }
  });

  // use dense iterations for now.
  size_t n_active = pbbslib::reduce(parlay::delayed_seq<uintE>(n, [&] (size_t i) { return active[i] != 0; }));
  std::cout << "nactive = " << n_active << std::endl;

  double one_plus_eps = 1 + eps;
  size_t rounds = 0;

  while (n_active > 0) {
    std::cout << "Starting round: " << rounds << std::endl;

    parallel_for(0, n, [&] (size_t i) {
    if (active[i] > 0) {
      // Decide whether red or blue
      if (rnd.ith_rand(i) & 1) {
        colors[i] = kRed;
      } else {
        colors[i] = kBlue;
      }
    }});

    std::cout << "num_blue = " << pbbslib::reduce(parlay::delayed_seq<uintE>(n, [&] (size_t i) { return colors[i] == kBlue; })) << std::endl;
    std::cout << "num_red = " << pbbslib::reduce(parlay::delayed_seq<uintE>(n, [&] (size_t i) { return colors[i] == kRed; })) << std::endl;

    auto merge_target = sequence<uintE>(n, UINT_E_MAX);

    parallel_for(0, n, [&] (size_t i) {
    if (active[i] > 0 && colors[i] == kBlue) {
      // Select a random heavy-edge incident to the vertex.
      uintE active_i = active[i];
      size_t edge_idx = rnd.ith_rand(n + i) % active_i;
      // seq selection for now (todo: could run in par for high-degrees)
      size_t k = 0;
      uintE ngh_id = std::numeric_limits<uintE>::max();
      W weight;  // TODO: don't really need to save?
      auto iter_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        if (Weights::get_weight(wgh, u, v, CG) >= lower_threshold) {
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
      if (colors[ngh_id] == kRed) {
        assert(CG.clusters[ngh_id].active);
        assert(CG.clusters[i].active);
        uintE ngh_cur_size = CG.clusters[ngh_id].cluster_size();
        uintE upper_bound = one_plus_eps * ngh_cur_size;
        uintE our_size = CG.clusters[i].cluster_size();

        // Enable to stress-test the merge implementation.
        //merge_target[i] = ngh_id;
        auto old_opt = pbbslib::fetch_and_add_threshold(
            &(CG.clusters[ngh_id].cas_size),
            our_size,
            upper_bound);
        if (old_opt.has_value()) {  // Success
          merge_target[i] = ngh_id;
        }
      }
    }});

    auto pairs = parlay::delayed_seq<std::pair<uintE, uintE>>(n, [&] (size_t i) { return std::make_pair(merge_target[i], i); });
    auto merges = parlay::filter(pairs, [&] (const auto& pair) { return pair.first != UINT_E_MAX; });

    std::cout << "Num merges = " << merges.size() << std::endl;

//    for (size_t i=0; i<CG.clusters.size(); i++) {
//      // check graph consistency (no edges to inactive vertices)
//      if (CG.clusters[i].active) {
//        auto f = [&] (const auto& u, const auto& v, const auto& wgh) {
//          assert(CG.clusters[v].active);
//        };
//        CG.clusters[i].iterate(i, f);
//      }
//    }
//    merges = sequence<std::pair<uintE, uintE>>(4);
//    merges[0] = {4, 0};
//    merges[1] = {4, 1};
//    merges[2] = {4, 2};
//    merges[3] = {4, 3};
    CG.template unite_merge<Sim>(std::move(merges));

    rnd = rnd.next();
    // Reset colors.
    // Recompute active and update n_active.

    timer rt; rt.start();
    parallel_for(0, n, [&] (size_t i) {
      colors[i] = 0;
//      if (active[i] > 0) {  // must have been active before to stay active.
      if (CG.clusters[i].active) {
        active[i] = get_num_active_edges(i);
      } else {
        active[i] = 0;
      }
    });
    n_active = pbbslib::reduce(parlay::delayed_seq<uintE>(n, [&] (size_t i) { return active[i] != 0; }));
    std::cout << "nactive is now = " << n_active << std::endl;
    rt.stop(); rt.reportTotal("reactivate time");

    rounds++;
  }
  std::cout << "Finished bucket." << std::endl;
  size_t rem_active = parlay::reduce(parlay::delayed_seq<uintE>(n, [&] (size_t i) { return CG.clusters[i].active; }));
  std::cout << rem_active << " vertices remain." << std::endl;
}



}  // namespace clustering
}  // namespace gbbs
