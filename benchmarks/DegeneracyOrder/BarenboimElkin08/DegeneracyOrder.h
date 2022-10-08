#pragma once

#include "gbbs/bucket.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"

#include "benchmarks/ApproximateDensestSubgraph/ApproxPeelingBKV12/DensestSubgraph.h"
#include "benchmarks/ApproximateDensestSubgraph/GreedyCharikar/DensestSubgraph.h"

namespace gbbs {
namespace barenboimelkin_degen {

template <class Graph>
inline sequence<uintE> DegeneracyOrder(Graph& GA, double epsilon = 0.1,
                                       bool approx = false) {
  double alpha = approx ? CharikarAppxDensestSubgraph(GA)
                        : WorkEfficientDensestSubgraph(GA, epsilon);
  const size_t n = GA.n;
  const size_t deg_cutoff =
      std::max((size_t)(ceil(alpha * epsilon)), (size_t)1);
  auto sortD = sequence<uintE>::from_function(n, [&](size_t i) { return i; });
  auto D = sequence<uintE>::from_function(
      n, [&](size_t i) { return GA.get_vertex(i).out_degree(); });
  auto em = EdgeMap<uintE, Graph>(GA, std::make_tuple(UINT_E_MAX, 0),
                                  (size_t)GA.m / 50);
  auto get_deg = [&](uintE p) -> uintE { return D[p] < deg_cutoff; };
  size_t start = 0;
  while (start < n) {
    // move all vert with deg < deg_cutoff in the front
    parlay::integer_sort_inplace(sortD.cut(start, n), get_deg);
    // radix::parallelIntegerSort(sortD.begin() + start, n - start, get_deg);
    auto BS = parlay::delayed_seq<size_t>(n - start, [&](size_t i) -> size_t {
      return D[sortD[i + start]] < deg_cutoff ? i + start : 0;
    });
    size_t end = parlay::reduce(BS, parlay::maxm<size_t>());
    if (end == start) end++;  // TODO step?

    auto num_removed = end - start;
    auto removed = sequence<uintE>::uninitialized(num_removed);
    parallel_for(0, num_removed,
                 [&](size_t i) { removed[i] = sortD[start + i]; });
    auto active = vertexSubset(n, std::move(removed));

    // least ns, from start to min(ns+start, n), is in order
    // update degrees based on peeled vert
    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const std::optional<std::tuple<uintE, uintE> > {
          uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
          D[v] -= edgesRemoved;
          return std::nullopt;
        };
    auto moved = em.template edgeMapCount_sparse<uintE>(active, apply_f);

    start = end;
  }
  auto ret = sequence<uintE>::uninitialized(n);
  parallel_for(0, n, [&](size_t j) { ret[sortD[j]] = j; });
  return ret;
}

}  // namespace barenboimelkin_degen
}  // namespace gbbs
