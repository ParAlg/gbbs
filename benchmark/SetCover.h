// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.


#pragma once

#include "bucket.h"
#include "edge_map_reduce.h"
#include "lib/index_map.h"
#include "lib/random.h"
#include "lib/random_shuffle.h"
#include "ligra.h"

namespace sc {
constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE COVERED = ((uintE)INT_E_MAX) - 1;
constexpr double epsilon = 0.01;
const double x = 1.0 / log(1.0 + sc::epsilon);

template <class W>
struct Visit_Elms {
  uintE* elms;
  uintE* perm;
  Visit_Elms(uintE* _elms, uintE* _perm) : elms(_elms), perm(_perm) {}
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    uintE oval = elms[d];
    uintE p_s = perm[s];
    writeMin(&(elms[d]), p_s);
    return false;
  }
  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    return updateAtomic(s, d, wgh);
  }
  inline bool cond(const uintE& d) const { return elms[d] != sc::COVERED; }
};
}  // namespace sc

template <template <class W> class vertex, class W>
dyn_arr<uintE> SetCover(graph<vertex<W>>& G, size_t num_buckets = 128) {
  auto Elms = array_imap<uintE>(G.n, [&](size_t i) { return UINT_E_MAX; });
  auto D =
      array_imap<uintE>(G.n, [&](size_t i) { return G.V[i].getOutDegree(); });
  auto get_bucket_clamped = [&](size_t deg) -> uintE {
    return (deg == 0) ? UINT_E_MAX : (uintE)floor(sc::x * log((double)deg));
  };
  auto bucket_f = [&](size_t i) { return get_bucket_clamped(D(i)); };
  auto b = make_buckets(G.n, bucket_f, decreasing, num_buckets);

  auto perm = array_imap<uintE>(G.n);
  timer bktt, packt, permt, emt;

  size_t rounds = 0;
  dyn_arr<uintE> cover = dyn_arr<uintE>();
  auto r = pbbs::default_random;
  while (true) {
    bktt.start();
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers;
    size_t cur_bkt = bkt.id;
    if (cur_bkt == b.null_bkt) {
      break;
    }
    bktt.stop();

    cout << "packing, active.size = " << active.size() << endl;
    packt.start();
    // 1. sets -> elements (Pack out sets and update their degree)
    auto pack_predicate = [&](const uintE& u, const uintE& ngh, const W& wgh) {
      return Elms[ngh] != sc::COVERED;
    };
    auto pack_apply = [&](uintE v, size_t ct) { D[v] = ct; };
    auto packed_vtxs = edgeMapFilter(G, active, pack_predicate, pack_edges);
    cout << "packed" << endl;
    vertexMap(packed_vtxs, pack_apply);
    packt.stop();

    // Find sets which still have sufficient degree (degree >= threshold)
    size_t threshold = ceil(pow(1.0 + sc::epsilon, cur_bkt));
    auto above_threshold = [&](const uintE& v, const uintE& deg) {
      return deg >= threshold;
    };
    auto still_active = vertexFilter2<uintE>(packed_vtxs, above_threshold);
    packed_vtxs.del();

    permt.start();
    // Update the permutation for the sets that are active in this round.
    still_active.toSparse();
    auto P = pbbs::random_permutation<uintE>(still_active.size(), r);
    parallel_for_bc(i, 0, still_active.size(), (still_active.size() > pbbs::kSequentialForThreshold), {
      uintE v = still_active.vtx(i);
      uintE pv = P[i];
      perm[v] = pv;
    });
    P.del();
    permt.stop();

    cout << "Round = " << rounds << " bkt = " << cur_bkt
         << " active = " << active.size()
         << " stillactive = " << still_active.size() << endl;

    emt.start();
    // 2. sets -> elements (writeMin to acquire neighboring elements)
    edgeMap(G, still_active, sc::Visit_Elms<W>(Elms.start(), perm.start()), -1,
            no_output | dense_forward);

    // 3. sets -> elements (count and add to cover if enough elms were won)
    const size_t low_threshold =
        std::max((size_t)ceil(pow(1.0 + sc::epsilon, cur_bkt - 1)), (size_t)1);
    auto won_ngh_f = wrap_f<W>([&](const uintE& u, const uintE& v) -> bool {
      return Elms[v] == perm[u];
    });
    auto threshold_f = [&](const uintE& v, const uintE& numWon) {
      if (numWon >= low_threshold) D[v] |= sc::TOP_BIT;
    };
    auto activeAndCts = edgeMapFilter(G, still_active, won_ngh_f);
    vertexMap(activeAndCts, threshold_f);
    auto inCover =
        vertexFilter2(activeAndCts, [&](const uintE& v, const uintE& numWon) {
          return numWon >= low_threshold;
        });
    cover.copyInF([&](uintE i) { return inCover.vtx(i); }, inCover.size());
    inCover.del();
    activeAndCts.del();

    // 4. sets -> elements (Sets that joined the cover mark their neighboring
    // elements as covered. Sets that didn't reset any acquired elements)
    auto reset_f = [&](const uintE& u, const uintE& v, const W& w) -> bool {
      if (Elms[v] == perm[u]) {
        if (D(u) & sc::TOP_BIT)
          Elms[v] = sc::COVERED;
        else
          Elms[v] = UINT_E_MAX;
      }
      return false;
    };
    edgeMap(G, still_active, ligra_utils::EdgeMap_F<W, decltype(reset_f)>(reset_f), -1,
            no_output | dense_forward);
    emt.stop();

    bktt.start();
    // Rebucket the active sets. Ignore those that joined the cover.
    active.toSparse();
    auto f = [&](size_t i) -> Maybe<tuple<uintE, uintE>> {
      const uintE v = active.vtx(i);
      const uintE dv = D(v);
      uintE bkt = UINT_E_MAX;
      if (!(dv & sc::TOP_BIT))
        bkt = b.get_bucket(cur_bkt, get_bucket_clamped(dv));
      return Maybe<tuple<uintE, uintE>>(make_tuple(v, bkt));
    };
    cout << "cover.size = " << cover.size << endl;
    b.update_buckets(f, active.size());
    active.del();
    still_active.del();
    rounds++;
    bktt.stop();
    r = r.next();
  }

  bktt.reportTotal("bucket");
  packt.reportTotal("pack");
  permt.reportTotal("perm");
  emt.reportTotal("emap");
  auto elm_cov = make_in_imap<uintE>(
      G.n, [&](uintE v) { return (uintE)(Elms[v] == sc::COVERED); });
  size_t elms_cov = pbbs::reduce_add(elm_cov);
  cout << "|V| = " << G.n << " |E| = " << G.m << endl;
  cout << "|cover|: " << cover.size << endl;
  cout << "Rounds: " << rounds << endl;
  cout << "Num_uncovered = " << (G.n - elms_cov) << endl;
  return std::move(cover);
}
