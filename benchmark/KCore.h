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

#include "bucket.h"
#include "edge_map_reduce.h"
#include "ligra.h"


template <template <typename W> class vertex, class W>
inline sequence<uintE> KCore(graph<vertex<W> >& GA, size_t num_buckets = 16) {
  const size_t n = GA.n;
  auto D =
      sequence<uintE>(n, [&](size_t i) { return GA.V[i].getOutDegree(); });

  auto em = EdgeMap<uintE, vertex, W>(GA, std::make_tuple(UINT_E_MAX, 0),
                                      (size_t)GA.m / 50);
  auto b = make_vertex_buckets(n, D, increasing, num_buckets);
  timer bt;

  size_t finished = 0, rho = 0, k_max = 0;
  while (finished != n) {
    bt.start();
    auto bkt = b.next_bucket();
    bt.stop();
    auto active = vertexSubset(n, bkt.identifiers);
    uintE k = bkt.id;
    finished += active.size();
    k_max = std::max(k_max, bkt.id);

    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const Maybe<std::tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
      uintE deg = D[v];
      if (deg > k) {
        uintE new_deg = std::max(deg - edgesRemoved, k);
        D[v] = new_deg;
        uintE bkt = b.get_bucket(new_deg);
        return wrap(v, bkt);
      }
      return Maybe<std::tuple<uintE, uintE> >();
    };

    vertexSubsetData<uintE> moved =
        em.template edgeMapCount<uintE>(active, apply_f);
    bt.start();
//    b.update_buckets(moved.get_fn_repr(), moved.size());
    if (moved.dense()) {
      b.update_buckets(moved.get_fn_repr(), n);
    } else {
      b.update_buckets(moved.get_fn_repr(), moved.size());
    }

    bt.stop();
    moved.del();
    active.del();
    rho++;
  }
  std::cout << "rho = " << rho << " k_{max} = " << k_max << "\n";
  bt.reportTotal("bucket time");
  return D;
}

template <class W>
struct kcore_fetch_add {
  uintE* er;
  uintE* D;
  uintE k;
  kcore_fetch_add(uintE* _er, uintE* _D, uintE _k) : er(_er), D(_D), k(_k) {}
  inline Maybe<uintE> update(const uintE& s, const uintE& d, const W& w) {
    er[d]++;
    if (er[d] == 1) {
      return Maybe<uintE>((uintE)0);
    }
    return Maybe<uintE>();
  }
  inline Maybe<uintE> updateAtomic(const uintE& s, const uintE& d,
                                   const W& wgh) {
    if (pbbslib::fetch_and_add(&er[d], (uintE)1) == 1) {
      return Maybe<uintE>((uintE)0);
    }
    return Maybe<uintE>();
  }
  inline bool cond(uintE d) { return D[d] > k; }
};

template <template <typename W> class vertex, class W>
inline sequence<uintE> KCore_FA(graph<vertex<W> >& GA,
                                  size_t num_buckets = 16) {
  const size_t n = GA.n;
  auto D =
      sequence<uintE>(n, [&](size_t i) { return GA.V[i].getOutDegree(); });
  auto ER = sequence<uintE>(n, [&](size_t i) { return 0; });

  auto b = make_vertex_buckets(n, D, increasing, num_buckets);

  size_t finished = 0;
  size_t rho = 0;
  size_t k_max = 0;
  while (finished != n) {
    auto bkt = b.next_bucket();
    auto active = vertexSubset(n, bkt.identifiers);
    uintE k = bkt.id;
    finished += active.size();
    k_max = std::max(k_max, bkt.id);

    auto apply_f = [&](const uintE v, uintE& bkt) -> void {
      uintE deg = D[v];
      uintE edgesRemoved = ER[v];
      ER[v] = 0;
      uintE new_deg = std::max(deg - edgesRemoved, k);
      D[v] = new_deg;
      bkt = b.get_bucket(deg, new_deg);
    };

    auto moved = edgeMapData<uintE>(
        GA, active, kcore_fetch_add<W>(ER.begin(), D.begin(), k));
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
  std::cout << "rho = " << rho << " k_{max} = " << k_max << "\n";
  return D;
}
