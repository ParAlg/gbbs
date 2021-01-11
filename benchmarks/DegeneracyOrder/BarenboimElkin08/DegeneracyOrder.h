#pragma once

#include "gbbs/bucket.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"
#include "gbbs/pbbslib/dyn_arr.h"
#include "pbbslib/integer_sort.h"

#include "benchmarks/ApproximateDensestSubgraph/GreedyCharikar/DensestSubgraph.h"
#include "benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12/DensestSubgraph.h"

namespace gbbs {
namespace barenboimelkin_degen {

template<class Graph>
inline sequence<uintE> DegeneracyOrder(Graph& GA, double epsilon=0.1, bool approx=false) {
  double alpha = approx ? CharikarAppxDensestSubgraph(GA) : WorkEfficientDensestSubgraph(GA, epsilon);
  const size_t n = GA.n;
  const size_t deg_cutoff = std::max((size_t) (ceil(alpha * epsilon)), (size_t) 1);
  auto sortD = sequence<uintE>(n, [&](size_t i) {
    return i;
  });
  auto D =
      sequence<uintE>(n, [&](size_t i) { return GA.get_vertex(i).out_degree(); });
  auto em = EdgeMap<uintE, Graph>(GA, std::make_tuple(UINT_E_MAX, 0),
                                      (size_t)GA.m / 50);
  auto get_deg =
      [&](uintE p) -> uintE { return D[p] < deg_cutoff; };
  size_t start = 0;
  while (start < n) {
    // move all vert with deg < deg_cutoff in the front
    integer_sort_inplace(sortD.slice(start, n), get_deg);
    //radix::parallelIntegerSort(sortD.begin() + start, n - start, get_deg);
    auto BS = pbbs::delayed_seq<size_t>(n - start, [&] (size_t i) -> size_t {
      return D[sortD[i + start]] < deg_cutoff ? i + start : 0;});
    size_t end = pbbs::reduce(BS, pbbs::maxm<size_t>());
    if (end == start) end++; //TODO step?
    // least ns, from start to min(ns+start, n), is in order
    // update degrees based on peeled vert
    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const std::optional<std::tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
        D[v] -= edgesRemoved;
      return std::nullopt;
    };
    auto active =
        vertexSubset(n, end - start, sortD.begin() + start);
    auto moved = em.template edgeMapCount_sparse<uintE>(active, apply_f);

    moved.del();
    start = end;
  }
  auto ret = sequence<uintE>::no_init(n);
  parallel_for (0,n,[&] (size_t j) { ret[sortD[j]] = j; });
  return ret;
}

}  // namespace barenboimelkin_degen
}  // namespace gbbs
