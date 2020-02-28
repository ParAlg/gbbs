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
#include "ligra/bucket.h"
#include "ligra/ligra.h"

constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;

struct Visit_F {
  sequence<uintE>& dists;
  Visit_F(sequence<uintE>& _dists) : dists(_dists) {}

  inline Maybe<uintE> update(const uintE& s, const uintE& d, const intE& w) {
    uintE oval = dists[d];
    uintE dist = oval | TOP_BIT, n_dist = (dists[s] | TOP_BIT) + w;
    if (n_dist < dist) {
      if (!(oval & TOP_BIT)) {  // First visitor
        dists[d] = n_dist;
        return Maybe<uintE>(oval);
      }
      dists[d] = n_dist;
    }
    return Maybe<uintE>();
  }

  inline Maybe<uintE> updateAtomic(const uintE& s, const uintE& d,
                                   const intE& w) {
    uintE oval = dists[d];
    uintE dist = oval | TOP_BIT;
    uintE n_dist = (dists[s] | TOP_BIT) + w;
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


template <class Graph>
void DeltaStepping(Graph& G, uintE src, uintE delta, size_t num_buckets=128) {
  size_t n = G.n;
  auto dists = pbbs::sequence<uintE>(n, [&] (size_t i) { return INT_E_MAX; });
  dists[src] = 0;

  auto get_bkt = [&] (const uintE& dist) -> uintE {
    return (dist == INT_E_MAX) ? UINT_E_MAX : (dist / delta); };
  auto get_ring = pbbslib::make_sequence<uintE>(n, [&] (const size_t& v) -> uintE {
    auto d = dists[v];
    return (d == INT_E_MAX) ? UINT_E_MAX : (d / delta); });
  auto b = make_vertex_buckets(n, get_ring, increasing, num_buckets);

  auto apply_f = [&] (const uintE v, uintE& oldDist) -> void {
    uintE newDist = dists[v] & VAL_MASK;
    dists[v] = newDist; // Remove the TOP_BIT in the distance.
    // Compute the previous bucket and new bucket for the vertex.
    uintE prev_bkt = get_bkt(oldDist), new_bkt = get_bkt(newDist);
    auto dest = b.get_bucket(prev_bkt, new_bkt);
    oldDist = dest; // write back
  };

  timer bktt;
  bktt.start();
  auto bkt = b.next_bucket();
  bktt.stop();
  flags fl = no_dense;
  while (bkt.id != b.null_bkt) {
    auto active = vertexSubset(n, bkt.identifiers);
    // The output of the edgeMap is a vertexSubsetData<uintE> where the value
    // stored with each vertex is its original distance in this round
    auto res = edgeMapData<uintE>(G, active, Visit_F(dists), G.m/20, fl);
    vertexMap(res, apply_f);
    bktt.start();
    if (res.dense()) {
      b.update_buckets(res.get_fn_repr(), n);
    } else {
      b.update_buckets(res.get_fn_repr(), res.size());
    }
    res.del();
    bkt = b.next_bucket();
    bktt.stop();
  }
  auto get_dist = [&] (size_t i) { return (dists[i] == INT_E_MAX) ? 0 : dists[i]; };
  auto dist_im = pbbs::delayed_seq<uintE>(n, get_dist);
  cout << "max_dist = " << pbbslib::reduce_max(dist_im) << endl;
  bktt.reportTotal("bucket time");
}

template <class Graph>
void Compute(Graph& G, commandLine P) {
  uintE src = P.getOptionLongValue("-src",0);
  uintE delta = P.getOptionLongValue("-delta",1);
  size_t num_buckets = P.getOptionLongValue("-nb", 128);
  if (num_buckets != (1 << pbbs::log2_up(num_buckets))) {
    cout << "Please specify a number of buckets that is a power of two" << endl;
    exit(-1);
  }
  cout << "### Application: Delta-Stepping" << endl;
  cout << "### Graph: " << P.getArgument(0) << endl;
  cout << "### Buckets: " << num_buckets << endl;
  cout << "### n: " << G.n << endl;
  cout << "### m: " << G.m << endl;
  cout << "### delta = " << delta << endl;
  DeltaStepping(G, src, delta, num_buckets);
}
