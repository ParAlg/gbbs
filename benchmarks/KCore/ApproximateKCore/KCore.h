#pragma once

#include <math.h>

#include "gbbs/gbbs.h"
#include "gbbs/julienne.h"

namespace gbbs {
namespace approximate_kcore {

template <class Graph>
inline sequence<uintE> KCore(Graph& G, size_t num_buckets = 16,
  double eps = 0.2, double delta = 0.1, bool use_pow = false) {
  const size_t n = G.n;

  auto Degrees =
      sequence<std::pair<uintE, bool>>::from_function(n, [&](size_t i) {
          return std::make_pair(G.get_vertex(i).out_degree(), false); });
  double one_plus_delta = log(1 + delta);
  auto get_bucket = [&](size_t deg) -> uintE {
    return ceil(log(1 + deg) / one_plus_delta);
  };
  auto D = sequence<uintE>::from_function(G.n, [&] (size_t i) {
    return get_bucket(Degrees[i].first); });

  auto em = hist_table<uintE, uintE>(std::make_tuple(UINT_E_MAX, 0), (size_t)G.m / 50);
  auto b = make_vertex_buckets(n, D, increasing, num_buckets);
  timer bt;

  size_t finished = 0, rho = 0, k_max = 0;

  size_t cur_inner_rounds = 0;
  size_t max_inner_rounds = log(G.n) / log(1.0 + eps);
  uintE prev_bkt = UINT_E_MAX;

  while (finished != n) {
    bt.start();
    auto bkt = b.next_bucket();
    bt.stop();
    auto active = vertexSubset(n, std::move(bkt.identifiers));
    uintE k = bkt.id;
    finished += active.size();
    k_max = std::max(k_max, bkt.id);
    if (k != prev_bkt) {
      prev_bkt = k;
      cur_inner_rounds = 0;
    }

    // Check if we hit the threshold for inner peeling rounds.
    if (cur_inner_rounds == max_inner_rounds) {
      // new re-insertions will go to at least bucket k (one greater than before).
      k++;
      cur_inner_rounds = 0;
    }

    // Mark peeled vertices as done.
    parallel_for(0, active.size(), [&] (size_t i) {
      uintE vtx = active.s[i];
      assert(!Degrees[vtx].second);  // not yet peeled
      Degrees[vtx].second = true;  // set to peeled
    });

    uintE lower_bound = ceil(pow((1 + delta), k-1));
    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const std::optional<std::tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edges_removed = std::get<1>(p);
      if (!Degrees[v].second) {
        uintE deg = Degrees[v].first;
        uintE new_deg = std::max(deg - edges_removed, lower_bound);
        assert(new_deg >= 0);
        Degrees[v].first = new_deg;
        uintE new_bkt = std::max(get_bucket(new_deg), k);
        uintE prev_bkt = D[v];
        if (prev_bkt != new_bkt) {
          D[v] = new_bkt;
          return wrap(v, b.get_bucket(new_bkt));
        }
      }
      return std::nullopt;
    };

    auto cond_f = [] (const uintE& u) { return true; };
    vertexSubsetData<uintE> moved = nghCount(G, active, cond_f, apply_f, em,
        no_dense);

    bt.start();
    if (moved.dense()) {
      b.update_buckets(moved.get_fn_repr(), n);
    } else {
      b.update_buckets(moved.get_fn_repr(), moved.size());
    }

    bt.stop();
    rho++;
    cur_inner_rounds++;
  }

  b.del();
  em.del();
  debug(bt.reportTotal("bucket time"););

  parallel_for(0, n, [&] (size_t i) {
    if (use_pow) {  // use 2^{peeled_bkt} as the coreness estimate
      D[i] = (D[i] == 0) ? 0 : 1 << D[i];
    } else {
      D[i] = Degrees[i].first;  // use capped induced degree when peeled as the coreness estimate
    }
  });
  return D;
}

}  // namespace approximate_kcore
}  // namespace gbbs
