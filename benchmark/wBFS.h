#pragma once

#include <cmath>
#include "bucket.h"
#include "lib/index_map.h"
#include "ligra.h"

namespace wbfs {
constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;

struct Visit_F {
  array_imap<uintE>& dists;
  Visit_F(array_imap<uintE>& _dists) : dists(_dists) {}

  inline Maybe<uintE> update(const uintE& s, const uintE& d, const intE& w) {
    uintE oval = dists.s[d];
    uintE dist = oval | TOP_BIT, n_dist = (dists.s[s] | TOP_BIT) + w;
    if (n_dist < dist) {
      if (!(oval & TOP_BIT)) {  // First visitor
        dists.s[d] = n_dist;
        return Maybe<uintE>(oval);
      }
      dists.s[d] = n_dist;
    }
    return Maybe<uintE>();
  }

  inline Maybe<uintE> updateAtomic(const uintE& s, const uintE& d,
                                   const intE& w) {
    uintE oval = dists.s[d];
    uintE dist = oval | TOP_BIT;
    uintE n_dist = (dists.s[s] | TOP_BIT) + w;
    if (n_dist < dist) {
      if (!(oval & TOP_BIT) &&
          CAS(&(dists.s[d]), oval, n_dist)) {  // First visitor
        return Maybe<uintE>(oval);
      }
      writeMin(&(dists.s[d]), n_dist);
    }
    return Maybe<uintE>();
  }

  inline bool cond(const uintE& d) const { return true; }
};

}  // namespace wbfs

template <
    template <typename W> class vertex, class W,
    typename std::enable_if<std::is_same<W, int32_t>::value, int>::type = 0>
auto wBFS(graph<vertex<W>>& G, uintE src, size_t num_buckets = 128,
          bool largemem = false, bool no_blocked = false) {
  auto before_state = get_pcm_state();
  timer t; t.start();

  timer init; init.start();
  auto V = G.V;
  size_t n = G.n, m = G.m;

  auto dists = array_imap<uintE>(n, [&](size_t i) { return INT_E_MAX; });
  dists[src] = 0;

  auto get_bkt = [&](const uintE& dist) -> const uintE {
    return (dist == INT_E_MAX) ? UINT_E_MAX : dist;
  };
  auto get_ring = [&](const size_t& v) -> const uintE {
    auto d = dists[v];
    return (d == INT_E_MAX) ? UINT_E_MAX : d;
  };
  auto b = make_buckets(n, get_ring, increasing, num_buckets);

  auto apply_f = [&](const uintE v, uintE& oldDist) -> void {
    uintE newDist = dists[v] & wbfs::VAL_MASK;
    dists[v] = newDist;  // Remove the TOP_BIT in the distance.
    // Compute the previous bucket and new bucket for the vertex.
    uintE prev_bkt = get_bkt(oldDist), new_bkt = get_bkt(newDist);
    bucket_dest dest = b.get_bucket(prev_bkt, new_bkt);
    oldDist = dest;  // write back
  };

  init.stop(); init.reportTotal("init time");
  timer bt, emt;
  auto bkt = b.next_bucket();
  size_t rd = 0;
  flags fl = dense_forward;
  if (!largemem) fl |= no_dense;
  if (!no_blocked) fl |= sparse_blocked;
  while (bkt.id != b.null_bkt) {
    auto active = bkt.identifiers;
    emt.start();
    // The output of the edgeMap is a vertexSubsetData<uintE> where the value
    // stored with each vertex is its original distance in this round
    auto res =
        edgeMapData<uintE>(G, active, wbfs::Visit_F(dists), G.m / 20, fl);
    vertexMap(res, apply_f);
    // update buckets with vertices that just moved
    emt.stop(); bt.start();
    if (res.dense()) {
      b.update_buckets(res.get_fn_repr(), n);
    } else {
      b.update_buckets(res.get_fn_repr(), res.size());
    }
    res.del(); active.del();
    bkt = b.next_bucket();
    bt.stop();
    rd++;
  }
  bt.reportTotal("bucket time");
  emt.reportTotal("edge map time");
  auto dist_im = make_in_imap<size_t>(
      n, [&](size_t i) { return (dists[i] == INT_E_MAX) ? 0 : dists[i]; });
  cout << "max dist = " << pbbs::reduce_max(dist_im) << endl;
  cout << "n rounds = " << rd << endl;

  double time_per_iter = t.stop();
  auto after_state = get_pcm_state();
  print_pcm_stats(before_state, after_state, 1, time_per_iter);

  return dists;
}

template <
    template <typename W> class vertex, class W,
    typename std::enable_if<!std::is_same<W, int32_t>::value, int>::type = 0>
auto wBFS(graph<vertex<W>>& G, uintE src, size_t num_buckets = 128,
          bool largemem = false, bool no_blocked = false) {
  auto dists = array_imap<uintE>(G.n, [&](size_t i) { return INT_E_MAX; });
  cout << "Unimplemented for unweighted graphs; use a regular BFS." << endl;
  exit(0);
  return dists;
}
