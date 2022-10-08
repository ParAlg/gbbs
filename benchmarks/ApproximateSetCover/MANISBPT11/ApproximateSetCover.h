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

#include "gbbs/gbbs.h" /* includes core gbbs libraries (graphs, graph operators, etc) */
#include "gbbs/julienne.h" /* includes bucketing data structure */

namespace gbbs {

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
    uintE p_s = perm[s];
    gbbs::write_min(&(elms[d]), p_s);
    return false;
  }
  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    return updateAtomic(s, d, wgh);
  }
  inline bool cond(const uintE& d) const { return elms[d] != sc::COVERED; }
};

}  // namespace sc

// Try implementing version with priority array.
// Question is how much does log(..) cost, each time when we unpack?

// build default bucket approach
// update_bucket should not take new bucket values. computes new bucket value
// internally.
// reductions are handled by atomic reduction operators; external to bucketing
// interface.

template <class Graph>
inline parlay::sequence<uintE> SetCover(Graph& G, size_t num_buckets = 512) {
  using W = typename Graph::weight_type;
  timer it;
  it.start();
  auto Elms =
      sequence<uintE>::from_function(G.n, [&](size_t i) { return UINT_E_MAX; });
  auto get_bucket_clamped = [&](size_t deg) -> uintE {
    return (deg == 0) ? UINT_E_MAX : (uintE)floor(sc::x * log((double)deg));
  };
  auto D = sequence<uintE>::from_function(G.n, [&](size_t i) {
    return get_bucket_clamped(G.get_vertex(i).out_degree());
  });
  auto b = make_vertex_buckets(G.n, D, decreasing, num_buckets);

  auto perm = sequence<uintE>::uninitialized(G.n);
  timer bktt, packt, permt, emt;

  timer nbt;
  size_t rounds = 0;
  parlay::sequence<uintE> cover;
  auto r = parlay::random();
  it.stop();
  it.next("initialization time");
  while (true) {
    nbt.start();
    auto bkt = b.next_bucket();
    auto active = vertexSubset(G.n, std::move(bkt.identifiers));
    size_t cur_bkt = bkt.id;
    if (cur_bkt == b.null_bkt) {
      break;
    }
    nbt.stop();

    packt.start();
    // 1. sets -> elements (Pack out sets and update their degree)
    auto pack_predicate = [&](const uintE& u, const uintE& ngh, const W& wgh) {
      return Elms[ngh] != sc::COVERED;
    };
    auto pack_apply = [&](uintE v, size_t ct) {
      D[v] = get_bucket_clamped(ct);
    };
    auto packed_vtxs = srcPack(G, active, pack_predicate, pack_edges);
    vertexMap(packed_vtxs, pack_apply);
    packt.stop();

    // Find sets which still have sufficient degree (degree >= threshold)
    size_t threshold = ceil(pow(1.0 + sc::epsilon, cur_bkt));
    auto above_threshold = [&](const uintE& v, const uintE& deg) {
      return deg >= threshold;
    };
    // auto still_active = vertexFilter_sparse(packed_vtxs, above_threshold);
    auto still_active = vertexFilter(packed_vtxs, above_threshold, no_dense);

    permt.start();
    // Update the permutation for the sets that are active in this round.
    still_active.toSparse();
    auto P = parlay::random_permutation<uintE>(still_active.size(), r);
    parallel_for(0, still_active.size(), kDefaultGranularity, [&](size_t i) {
      uintE v = still_active.vtx(i);
      uintE pv = P[i];
      perm[v] = pv;
    });
    P.clear();
    permt.stop();

    debug(std::cout << "Round = " << rounds << " bkt = " << cur_bkt
                    << " active = " << active.size()
                    << " stillactive = " << still_active.size() << "\n";);

    emt.start();
    // 2. sets -> elements (write_min to acquire neighboring elements)
    edgeMap(G, still_active, sc::Visit_Elms<W>(Elms.begin(), perm.begin()), -1,
            no_output | dense_forward);

    // 3. sets -> elements (count and add to cover if enough elms were won)
    const size_t low_threshold =
        std::max((size_t)ceil(pow(1.0 + sc::epsilon, cur_bkt - 1)), (size_t)1);
    auto won_ngh_f = [&](const uintE& u, const uintE& v, const W& wgh) -> bool {
      return Elms[v] == perm[u];
    };
    auto threshold_f = [&](const uintE& v, const uintE& numWon) {
      if (numWon >= low_threshold) D[v] = UINT_E_MAX;
    };
    auto activeAndCts = srcCount(G, still_active, won_ngh_f);
    vertexMap(activeAndCts, threshold_f);
    auto inCover = vertexFilter(activeAndCts,
                                [&](const uintE& v, const uintE& numWon) {
                                  // vertexFilter_sparse(activeAndCts, [&](const
                                  // uintE& v, const uintE& numWon) {
                                  return numWon >= low_threshold;
                                },
                                no_dense);

//     cover.copyInF([&](uintE i) { return inCover.vtx(i); }, inCover.size());

     auto to_copy = parlay::delayed_seq<uintE>(inCover.size(), [&](uintE i) { return inCover.vtx(i); });
     cover.append(to_copy);

    // 4. sets -> elements (Sets that joined the cover mark their neighboring
    // elements as covered. Sets that didn't reset any acquired elements)
    auto reset_f = [&](const uintE& u, const uintE& v, const W& w) -> bool {
      if (Elms[v] == perm[u]) {
        if (D[u] == UINT_E_MAX)
          Elms[v] = sc::COVERED;
        else
          Elms[v] = UINT_E_MAX;
      }
      return false;
    };
    edgeMap(G, still_active, EdgeMap_F<W, decltype(reset_f)>(reset_f), -1,
            no_output | dense_forward);
    emt.stop();

    bktt.start();
    // Rebucket the active sets. Ignore those that joined the cover.
    active.toSparse();
    auto f = [&](size_t i) -> std::optional<std::tuple<uintE, uintE>> {
      const uintE v = active.vtx(i);
      const uintE v_bkt = D[v];
      uintE bucket = UINT_E_MAX;
      if (!(v_bkt == UINT_E_MAX)) bucket = b.get_bucket(v_bkt);
      return std::optional<std::tuple<uintE, uintE>>(
          std::make_tuple(v, bucket));
    };
    debug(std::cout << "cover.size = " << cover.size() << "\n");
    b.update_buckets(f, active.size());
    rounds++;
    bktt.stop();
    r = r.next();
  }

  bktt.next("bucket");
  nbt.next("next bucket time");
  packt.next("pack");
  permt.next("perm");
  emt.next("emap");
  auto elm_cov_f = [&](uintE v) { return (uintE)(Elms[v] == sc::COVERED); };
  auto elm_cov = parlay::delayed_seq<uintE>(G.n, elm_cov_f);
  size_t elms_cov = parlay::reduce(elm_cov);
  std::cout << "|V| = " << G.n << " |E| = " << G.m << "\n";
  std::cout << "|cover|: " << cover.size() << "\n";
  std::cout << "Rounds: " << rounds << "\n";
  std::cout << "Num_uncovered = " << (G.n - elms_cov) << "\n";
  return cover;
}

}  // namespace gbbs
