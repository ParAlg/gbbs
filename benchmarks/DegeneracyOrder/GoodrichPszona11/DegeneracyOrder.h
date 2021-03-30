#pragma once

#include "gbbs/bucket.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"
#include "gbbs/pbbslib/dyn_arr.h"
#include "pbbslib/integer_sort.h"
#include "pbbslib/kth_smallest.h"
#include "pbbslib/random.h"

namespace gbbs {
namespace goodrichpszona_degen {

// Goodrich (2+epsilon) approx for degeneracy ordering where epsilon > 0
// Returns vertice sorted in degeneracy order
template<class Graph>
inline sequence<uintE> DegeneracyOrder(Graph& GA, double epsilon=0.1) {
  const size_t n = GA.n;
  const size_t ns = std::max((size_t) (ceil((n*epsilon) / (2+epsilon))), (size_t) 1);

  auto active = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });

  /* induced degrees sequence */
  auto D =
      sequence<uintE>(n, [&](size_t i) { return GA.get_vertex(i).out_degree(); });

  auto em = EdgeMap<uintE, Graph>(GA, std::make_tuple(UINT_E_MAX, 0),
                                      (size_t)GA.m / 20);

  auto ret = pbbslib::dyn_arr<uintE>(n);

  pbbs::random r;
  timer kt, ft;
  while (active.size() > 0) {
    /* compute cutoff using kth-smallest */

    auto active_degs = pbbslib::make_sequence<uintE>(active.size(), [&] (size_t i) {
      uintE v = active[i];
      return D[v];
    });
    debug(
    std::cout << "Kth smallesting w ns = " << ns << std::endl;
    std::cout << "num remaining = " << active_degs.size() << std::endl;);
    kt.start();
    uintE threshold = pbbs::approximate_kth_smallest(active_degs, ns, std::less<uintE>(), r);
    r = r.next();
    kt.stop();


    auto lte_threshold = [&] (const uintE& v) {
      return D[v] <= threshold;
    };
    auto gt_threshold = [&] (const uintE& v) {
      return D[v] > threshold;
    };

    ft.start();
    auto this_round = pbbs::filter(active, lte_threshold);
    active = pbbs::filter(active, gt_threshold);
    ft.stop();

    ret.copyInF([&] (size_t i) { return this_round[i]; }, this_round.size());

    // least ns, from start to min(ns+start, n), is in order
    // update degrees based on peeled vert
    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const std::optional<std::tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
      D[v] -= edgesRemoved;
      return std::nullopt;
    };
    auto this_round_vs = vertexSubset(n, std::move(this_round));
    auto moved = em.template edgeMapCount_sparse<uintE>(this_round_vs, apply_f);
  }
  pbbs::sequence<uintE> ret_seq(ret.A, n);
  ret.A = nullptr; ret.alloc = false; /* sketchy */
  debug(
  kt.reportTotal("kth time");
  ft.reportTotal("filter time"););
  return ret_seq;
}

// Goodrich (2+epsilon) approx for degeneracy ordering where epsilon > 0
// Returns vertice sorted in degeneracy order
template<class Graph>
inline sequence<uintE> DegeneracyOrder_intsort(Graph& GA, double epsilon=0.001) {
  const size_t n = GA.n;
  const size_t ns = std::max((size_t) (ceil((n*epsilon) / (2+epsilon))), (size_t) 1);

  auto sortD = sequence<uintE>(n, [&](size_t i) {
    return i;
  });
  auto D =
      sequence<uintE>(n, [&](size_t i) { return GA.get_vertex(i).out_degree(); });
  auto em = EdgeMap<uintE, Graph>(GA, std::make_tuple(UINT_E_MAX, 0),
                                      (size_t)GA.m / 20);
  auto get_deg =
      [&](uintE p) -> uintE { return D[p]; };
  for (size_t start = 0; start < n; start += ns) {
    // sort vertices in GA by degree, from start to n
    integer_sort_inplace(sortD.slice(start, n), get_deg);
    //radix::parallelIntegerSort(sortD.begin() + start, n - start, get_deg);
    // uintE deg_max = D[sortD[std::min(ns + start, n)]];

    // least ns, from start to min(ns+start, n), is in order
    // update degrees based on peeled vert
    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const std::optional<std::tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
      //if (D[v] >= deg_max) {
        D[v] -= edgesRemoved;
      //  return wrap(v, D[v]);
      //}
      return std::nullopt;
    };

    size_t num_removed = std::min(ns + start, n) - start;
    auto removed = sequence<uintE>::no_init(num_removed);
    parallel_for(0, num_removed, [&] (size_t i) {
      removed[i] = sortD[start + i];
    });
    auto active =
        vertexSubset(n, std::move(removed));
    auto moved = em.template edgeMapCount_sparse<uintE>(active, apply_f);
  }
  auto ret = sequence<uintE>::no_init(n);
  parallel_for (0,n,[&] (size_t j) { ret[sortD[j]] = j; });
  return ret;
}

}  // namespace goodrichpszona_degen
}  // namespace gbbs
