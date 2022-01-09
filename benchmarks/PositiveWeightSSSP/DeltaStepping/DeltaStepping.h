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

#include <cmath>
#include "gbbs/bucket.h"
#include "gbbs/gbbs.h"

namespace gbbs {

template <class W, class Distance>
struct Visit_F {
  sequence<std::pair<Distance, bool>>& dists;
  Visit_F(sequence<std::pair<Distance, bool>>& _dists) : dists(_dists) {}

  inline std::optional<Distance> update(const uintE& s, const uintE& d,
                                        const W& w) {
    Distance dist = dists[d].first;
    Distance n_dist;
    if
      constexpr(std::is_same<W, gbbs::empty>()) { n_dist = dists[s].first + 1; }
    else {
      n_dist = dists[s].first + w;
    }
    if (n_dist < dist) {
      if (!dists[d].second) {  // First visitor
        dists[d] = {n_dist, true};
        return std::optional<Distance>(dist);
      }
      dists[d].first = n_dist;
    }
    return std::nullopt;
  }

  inline std::optional<Distance> updateAtomic(const uintE& s, const uintE& d,
                                              const W& w) {
    Distance dist = dists[d].first;
    Distance n_dist;
    if
      constexpr(std::is_same<W, gbbs::empty>()) { n_dist = dists[s].first + 1; }
    else {
      n_dist = dists[s].first + w;
    }
    if (n_dist < dist) {
      gbbs::write_min(&(dists[d].first), n_dist);
      if (!dists[d].second &&
          gbbs::atomic_compare_and_swap(&dists[d].second, false,
                                        true)) {  // First visitor
        return std::optional<Distance>(dist);
      }
    }
    return std::nullopt;
  }

  inline bool cond(const uintE& d) const { return true; }
};

template <class Graph>
auto DeltaStepping(Graph& G, uintE src, double delta,
                   size_t num_buckets = 128) {
  using W = typename Graph::weight_type;
  using Distance =
      typename std::conditional<std::is_same<W, gbbs::empty>::value, uintE,
                                W>::type;
  constexpr Distance kMaxWeight = std::numeric_limits<Distance>::max();
  size_t n = G.n;
  std::cout << "Using delta = " << delta << std::endl;
  auto dists = sequence<std::pair<Distance, bool>>::from_function(
      n, [&](size_t i) { return std::make_pair(kMaxWeight, false); });
  dists[src] = {(Distance)0, false};
  auto bkts = sequence<uintE>(n, UINT_E_MAX);

  auto get_bkt = [&](const Distance& dist) -> uintE {
    return (dist == kMaxWeight) ? UINT_E_MAX : (uintE)(dist / delta);
  };
  auto get_ring = parlay::delayed_seq<uintE>(n, [&](const size_t& v) -> uintE {
    auto d = dists[v].first;
    return (d == kMaxWeight) ? UINT_E_MAX : (uintE)(d / delta);
  });
  auto b = make_vertex_buckets(n, get_ring, increasing, num_buckets);

  auto apply_f = [&](const uintE v, const Distance& oldDist) -> void {
    Distance newDist = dists[v].first;
    dists[v] = {newDist, false};  // reset flag
    // Compute the previous bucket and new bucket for the vertex.
    uintE prev_bkt = get_bkt(oldDist), new_bkt = get_bkt(newDist);
    bkts[v] = b.get_bucket(prev_bkt, new_bkt);
  };

  timer bktt;
  bktt.start();
  auto bkt = b.next_bucket();
  bktt.stop();
  flags fl = no_dense;
  size_t round = 0;
  while (bkt.id != b.null_bkt) {
    round++;
    auto active = vertexSubset(n, std::move(bkt.identifiers));
    // The output of the edgeMap is a vertexSubsetData<Distance> where the value
    // stored with each vertex is its original distance in this round
    auto res = edgeMapData<Distance>(G, active, Visit_F<W, Distance>(dists),
                                     G.m / 20, fl);
    vertexMap(res, apply_f);
    bktt.start();
    if (res.dense()) {
      b.update_buckets(
          [&](size_t i) -> std::optional<std::tuple<uintE, uintE>> {
            if (std::get<0>(res.d[i])) {
              return std::optional<std::tuple<uintE, uintE>>(
                  std::make_tuple(i, bkts[i]));
            } else {
              return std::nullopt;
            }
          },
          G.n);
    } else {
      b.update_buckets(
          [&](size_t i) {
            auto v = std::get<0>(res.s[i]);
            return std::optional<std::tuple<uintE, uintE>>(
                std::make_tuple(v, bkts[v]));
          },
          res.size());
    }
    bkt = b.next_bucket();
    bktt.stop();
  }
  auto get_dist = [&](size_t i) {
    return (dists[i].first == kMaxWeight) ? 0 : dists[i].first;
  };
  auto dist_im = parlay::delayed_seq<Distance>(n, get_dist);
  std::cout << "max_dist = " << parlay::reduce_max(dist_im) << std::endl;
  bktt.next("bucket time");
  auto ret = sequence<Distance>::from_function(
      n, [&](size_t i) { return dists[i].first; });
  return ret;
}

}  // namespace gbbs
