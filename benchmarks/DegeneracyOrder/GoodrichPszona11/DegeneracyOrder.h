#include "ligra/bucket.h"
#include "ligra/edge_map_reduce.h"
#include "ligra/ligra.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "pbbslib/integer_sort.h"

namespace goodrichpszona_degen {
// Goodrich (2+epsilon) approx for degeneracy ordering where epsilon > 0
// Returns vertice sorted in degeneracy order
template<class Graph>
inline sequence<uintE> DegeneracyOrder(Graph& GA, double epsilon=0.001) {
  const size_t n = GA.n;
  const size_t ns = std::max((size_t) (ceil((n*epsilon) / (2+epsilon))), (size_t) 1);

  auto sortD = sequence<uintE>(n, [&](size_t i) {
    return i;
  });
  auto D =
      sequence<uintE>(n, [&](size_t i) { return GA.get_vertex(i).getOutDegree(); });
  auto em = EdgeMap<uintE, Graph>(GA, std::make_tuple(UINT_E_MAX, 0),
                                      (size_t)GA.m / 50);
  auto get_deg =
      [&](uintE& p) -> uintE { return D[p]; };
  for (size_t start = 0; start < n; start += ns) {
    // sort vertices in GA by degree, from start to n
    integer_sort_inplace(sortD.slice(start, n), get_deg);
    //radix::parallelIntegerSort(sortD.begin() + start, n - start, get_deg);
    uintE deg_max = D[sortD[std::min(ns + start, n)]];
    
    // least ns, from start to min(ns+start, n), is in order
    // update degrees based on peeled vert
    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const Maybe<std::tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
      //if (D[v] >= deg_max) {
        D[v] -= edgesRemoved;
      //  return wrap(v, D[v]);
      //}
      return Maybe<std::tuple<uintE, uintE> >();
    };
    auto active =
        vertexSubset(n, std::min(ns + start, n) - start, sortD.begin() + start);
    auto moved = em.template edgeMapCount_sparse<uintE>(active, apply_f);

    moved.del();
  }
  auto ret = sequence<uintE>::no_init(n);
  parallel_for (0,n,[&] (size_t j) { ret[sortD[j]] = j; });
  return ret;
} 
}