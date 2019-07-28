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
#include "bucket.h"
#include "ligra.h"

namespace wbfs {
constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;

template <class W, class GW>
struct Visit_F {
  sequence<uintE>& dists;
  GW& get_weight;
  Visit_F(sequence<uintE>& dists, GW& get_weight) : dists(dists), get_weight(get_weight) {}

  inline Maybe<uintE> update(const uintE& s, const uintE& d, const W& wgh) {
    uintE edgeLen = get_weight(s, d, wgh);
    uintE oval = dists[d];
    uintE dist = oval | TOP_BIT, n_dist = (dists[s] | TOP_BIT) + edgeLen;
    if (n_dist < dist) {
      if (!(oval & TOP_BIT)) {  // First visitor
        dists[d] = n_dist;
        return Maybe<uintE>(oval);
      }
      dists[d] = n_dist;
    }
    return Maybe<uintE>();
  }

  inline Maybe<uintE> updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    uintE edgeLen = get_weight(s, d, wgh);
    uintE oval = dists[d];
    uintE dist = oval | TOP_BIT;
    uintE n_dist = (dists[s] | TOP_BIT) + edgeLen;
    if (n_dist < dist) {
      if (!(oval & TOP_BIT) &&
          pbbslib::atomic_compare_and_swap(&(dists[d]), oval, n_dist)) {  // First visitor
        return Maybe<uintE>(oval);
      }
      pbbslib::write_min(&(dists[d]), n_dist);
    }
    return Maybe<uintE>();
  }

  inline bool cond(const uintE& d) const { return true; }
};

template <class W, class GW>
Visit_F<W, GW> make_visit_f(sequence<uintE>& dists, GW& get_weight) {
  return Visit_F<W, GW>(dists, get_weight);
}

}  // namespace wbfs

template <class G>
inline sequence<uintE> wBFS(G& GA, uintE src,
                            size_t num_buckets = 128, bool largemem = false,
                            bool no_blocked = false) {
  using W = typename G::weight_type;
//  auto before_state = get_pcm_state();
  timer t;
  t.start();

  timer init;
  init.start();
  size_t n = GA.n;

  auto dists = sequence<uintE>(n, [&](size_t i) { return INT_E_MAX; });
  dists[src] = 0;

  auto get_bkt = [&](const uintE& dist) -> const uintE {
    return (dist == INT_E_MAX) ? UINT_E_MAX : dist;
  };
  auto get_ring = pbbslib::make_sequence<uintE>(n, [&](const size_t& v) -> const uintE {
    auto d = dists[v];
    return (d == INT_E_MAX) ? UINT_E_MAX : d;
  });
  auto b = make_vertex_buckets(n, get_ring, increasing, num_buckets);

  auto apply_f = [&](const uintE v, uintE& oldDist) -> void {
    uintE newDist = dists[v] & wbfs::VAL_MASK;
    dists[v] = newDist;  // Remove the TOP_BIT in the distance.
    // Compute the previous bucket and new bucket for the vertex.
    uintE prev_bkt = get_bkt(oldDist), new_bkt = get_bkt(newDist);
    uintE dest_bkt = b.get_bucket(prev_bkt, new_bkt);
    oldDist = dest_bkt;  // write back
  };

  auto degree_f = [&](size_t i) {
    auto vtx = GA.get_vertex(i);
    return std::max(vtx.getInDegree(), vtx.getOutDegree());
  };
  auto degree_im = pbbslib::make_sequence<size_t>(n, degree_f);
  size_t max_degree = pbbslib::reduce_max(degree_im);
  size_t normalize = 2*max_degree + 1;

  // edge -> custom_weight
  auto get_weight = [&] (const uintE& u, const uintE& v, const W& wgh) -> intE {
    uintE deg_u = GA.get_vertex(u).getOutDegree();
    uintE deg_v = GA.get_vertex(v).getOutDegree();
    return pbbs::log2_up((size_t)((1/static_cast<double>(deg_u + deg_v + 1))*normalize));
  };


  init.stop();
  init.reportTotal("init time");
  timer bt, emt;
  auto bkt = b.next_bucket();
  size_t rd = 0;
  flags fl = dense_forward;
  if (!largemem) fl |= no_dense;
  if (!no_blocked) fl |= sparse_blocked;
  while (bkt.id != b.null_bkt) {
    auto active = vertexSubset(n, bkt.identifiers);
    emt.start();
    // The output of the edgeMap is a vertexSubsetData<uintE> where the value
    // stored with each vertex is its original distance in this round
    auto res =
        edgeMapData<uintE>(GA, active, wbfs::make_visit_f<W>(dists, get_weight), GA.m / 20, fl);
    vertexMap(res, apply_f);
    // update buckets with vertices that just moved
    emt.stop();
    bt.start();
    if (res.dense()) {
      b.update_buckets(res.get_fn_repr(), n);
    } else {
      b.update_buckets(res.get_fn_repr(), res.size());
    }
    res.del();
    active.del();
    bkt = b.next_bucket();
    bt.stop();
    rd++;
  }
  bt.reportTotal("bucket time");
  emt.reportTotal("edge map time");
  auto dist_f = [&](size_t i) { return (dists[i] == INT_E_MAX) ? 0 : dists[i]; };
  auto dist_im = pbbslib::make_sequence<size_t>(n, dist_f);
  std::cout << "max dist = " << pbbslib::reduce_max(dist_im) << "\n";
  std::cout << "n rounds = " << rd << "\n";

  double time_per_iter = t.stop();

//  auto after_state = get_pcm_state();
//  print_pcm_stats(before_state, after_state, 1, time_per_iter);
  return dists;
}
