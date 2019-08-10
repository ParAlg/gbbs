// Usage:
// numactl -i all ./UnionFindCC -rounds 3 -s -m twitter_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc

#include <algorithm>
#include "ligra.h"
#include "utils/union_find.h"
#include "utils/stats.h"

// Assumes root is negative
// Not making parent array volatile improves
// performance and doesn't affect correctness
inline uintE uf_find(uintE i, const pbbs::sequence<uintE>& parent) {
  uintE j = i;
  if (parent[j] == UINT_E_MAX) return j;
  do j = parent[j];
  while (parent[j] < UINT_E_MAX);
  //note: path compression can happen in parallel in the same tree, so
  //only link from smaller to larger to avoid cycles
  uintE tmp;
  while((tmp=parent[i])<j){ parent[i]=j; i=tmp;}
  return j;
}

// returns a component labeling
template <class G>
pbbs::sequence<uintE> UnionFindCC(G& GA) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;
  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return UINT_E_MAX; });
  auto IDs = pbbs::sequence<uintE>(n, [&] (size_t i) { return UINT_E_MAX; });

  parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (uintE u, uintE v, const W& wgh) {
      if (u < v) {
        while(1) {
          u = uf_find(u, parents);
          v = uf_find(v, parents);
          if(u == v) break;
          if(u > v) std::swap(u,v);
          if(IDs[u] == UINT_E_MAX &&
              pbbs::atomic_compare_and_swap(&IDs[u],UINT_E_MAX,v)) {
            parents[u] = v;
            break;
          }
        }
      }
    };
    GA.V[i].mapOutNgh(i, map_f); // in parallel
  });

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    IDs[i] = uf_find(i,parents);
  });
  ft.stop(); ft.reportTotal("find time");
  return IDs;
}

template <class G>
double UnionFindCC_runner(G& GA, commandLine P) {
  auto beta = P.getOptionDoubleValue("-beta", 0.2);
  std::cout << "### Application: UnionFindCC (Connectivity)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: " << std::endl;
  std::cout << "### -----------------------------------" << endl;

  timer t;
  t.start();
  auto components = UnionFindCC(GA);
  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;

  if (P.getOption("-stats")) {
    num_cc(components);
    largest_cc(components);
  }
  return tt;
}

generate_main(UnionFindCC_runner, false);
