#pragma once

#include "bucket.h"
#include "edge_map_reduce.h"
#include "ligra.h"

#include "lib/index_map.h"

template <template <typename W> class vertex, class W>
array_imap<uintE> KCore(graph<vertex<W> >& GA, size_t num_buckets = 16) {
  const size_t n = GA.n;
  const size_t m = GA.m;
  auto D =
      array_imap<uintE>(n, [&](size_t i) { return GA.V[i].getOutDegree(); });

  auto em = EdgeMap<uintE, vertex, W>(GA, make_tuple(UINT_E_MAX, 0),
                                      (size_t)GA.m / 50);
  auto b = make_buckets(n, D, increasing, num_buckets);
  timer bt;

  size_t finished = 0, rho = 0, k_max = 0;
  while (finished != n) {
    bt.start();
    auto bkt = b.next_bucket();
    bt.stop();
    auto active = bkt.identifiers;
    uintE k = bkt.id;
    finished += active.size();
    k_max = std::max(k_max, bkt.id);

    auto apply_f =
        [&](const tuple<uintE, uintE>& p) -> const Maybe<tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
      uintE deg = D[v];
      if (deg > k) {
        uintE new_deg = max(deg - edgesRemoved, k);
        D[v] = new_deg;
        uintE bkt = b.get_bucket(deg, new_deg);
        return wrap(v, bkt);
      }
      return Maybe<tuple<uintE, uintE> >();
    };

    vertexSubsetData<uintE> moved =
        em.template edgeMapCount<uintE>(active, apply_f);
    bt.start();
    b.update_buckets(moved.get_fn_repr(), moved.size());
    bt.stop();
    moved.del();
    active.del();
    rho++;
  }
  cout << "rho = " << rho << " k_{max} = " << k_max << endl;
  bt.reportTotal("bucket time");
  return std::move(D);
}

template <class W>
struct kcore_fetch_add {
  uintE* er;
  uintE* D;
  uintE k;
  kcore_fetch_add(uintE* _er, uintE* _D, uintE _k) : er(_er), D(_D), k(_k) {}
  inline auto update(const uintE& s, const uintE& d, const W& w) {
    er[d]++;
    if (er[d] == 1) {
      return Maybe<uintE>((uintE)0);
    }
    return Maybe<uintE>();
  }
  inline auto updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    if (writeAdd(&er[d], (uintE)1) == 1) {
      return Maybe<uintE>((uintE)0);
    }
    return Maybe<uintE>();
  }
  inline bool cond(uintE d) { return D[d] > k; }
};

template <template <typename W> class vertex, class W>
auto KCore_FA(graph<vertex<W> >& GA, size_t num_buckets = 16) {
  const size_t n = GA.n;
  const size_t m = GA.m;
  auto D =
      array_imap<uintE>(n, [&](size_t i) { return GA.V[i].getOutDegree(); });
  auto ER = array_imap<uintE>(n, [&](size_t i) { return 0; });

  auto b = make_buckets(n, D, increasing, num_buckets);

  size_t finished = 0;
  size_t rho = 0;
  size_t k_max = 0;
  while (finished != n) {
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    uintE k = bkt.id;
    finished += active.size();
    k_max = std::max(k_max, bkt.id);

    auto apply_f = [&](const uintE v, uintE& bkt) -> void {
      uintE deg = D[v];
      uintE edgesRemoved = ER[v];
      ER[v] = 0;
      uintE new_deg = max(deg - edgesRemoved, k);
      D[v] = new_deg;
      bkt = b.get_bucket(deg, new_deg);
    };

    auto moved = edgeMapData<uintE>(
        GA, active, kcore_fetch_add<W>(ER.start(), D.start(), k));
    vertexMap(moved, apply_f);

    if (moved.dense()) {
      b.update_buckets(moved.get_fn_repr(), n);
    } else {
      b.update_buckets(moved.get_fn_repr(), moved.size());
    }
    moved.del();
    active.del();
    rho++;
  }
  cout << "rho = " << rho << " k_{max} = " << k_max << endl;
  return std::move(D);
}
