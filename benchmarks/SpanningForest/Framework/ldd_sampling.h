#pragma once

#include "benchmarks/LowDiameterDecomposition/MPX13/LowDiameterDecomposition.h"

namespace gbbs {
template <class W>
struct LDD_Edges_Fn {
  sequence<uintE>& Parents;
  sequence<edge>& Edges;

  LDD_Edges_Fn(sequence<uintE>& Parents, sequence<edge>& Edges)
      : Parents(Parents), Edges(Edges) {}

  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    Parents[d] = Parents[s];
    Edges[d] = std::make_pair(s, d);
    return true;
  }

  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    if (gbbs::atomic_compare_and_swap(&Parents[d], UINT_E_MAX, Parents[s])) {
      Edges[d] = std::make_pair(s, d);
      return true;
    }
    return false;
  }

  inline bool cond(uintE d) { return Parents[d] == UINT_E_MAX; }
};

// Returns a pair containing the clusters and edges.
template <class Graph>
inline std::pair<sequence<uintE>, sequence<edge>> LDD_sample_edges(
    Graph& G, double beta, bool permute = true, bool pack = false) {
  using W = typename Graph::weight_type;
  size_t n = G.n;

  sequence<uintE> vertex_perm;
  if (permute) {
    vertex_perm = parlay::random_permutation<uintE>(n);
  }
  auto shifts = ldd_utils::generate_shifts(n, beta);
  auto Parents = sequence<uintE>(n);
  parallel_for(0, n, kDefaultGranularity,
               [&](size_t i) { Parents[i] = UINT_E_MAX; });

  auto Edges = sequence<edge>(n, empty_edge);

  size_t round = 0, num_visited = 0;
  vertexSubset frontier(n);  // Initially empty
  size_t num_added = 0;
  while (num_visited < n) {
    size_t start = shifts[round];
    size_t end = std::min(static_cast<size_t>(shifts[round + 1]), n);
    size_t num_to_add = end - start;
    if (num_to_add > 0) {
      assert((num_added + num_to_add) <= n);
      auto candidates_f = [&](size_t i) {
        if (permute)
          return vertex_perm[num_added + i];
        else
          return static_cast<uintE>(num_added + i);
      };
      auto candidates = parlay::delayed_seq<uintE>(num_to_add, candidates_f);
      auto pred = [&](uintE v) { return Parents[v] == UINT_E_MAX; };
      auto new_centers = parlay::filter(candidates, pred);
      add_to_vsubset(frontier, new_centers.begin(), new_centers.size());
      parallel_for(0, new_centers.size(), kDefaultGranularity, [&](size_t i) {
        uintE v = new_centers[i];
        Parents[v] = v;
      });
      num_added += num_to_add;
    }

    num_visited += frontier.size();
    if (num_visited >= n) break;

    auto ldd_f = LDD_Edges_Fn<W>(Parents, Edges);
    vertexSubset next_frontier =
        edgeMap(G, frontier, ldd_f, -1, sparse_blocked);
    frontier = std::move(next_frontier);

    round++;
  }
  return std::make_pair(Parents, Edges);
}

template <class G>
struct LDDSamplingTemplate {
  G& GA;

  LDDSamplingTemplate(G& GA, commandLine& P) : GA(GA) {}

  auto initial_spanning_forest() {
    size_t n = GA.n;

    timer lddt;
    lddt.start();
    auto[Parents, Edges] = LDD_sample_edges(GA, 0.2, /* permute = */ false);
    lddt.stop();
    lddt.next("## ldd time");

    return std::make_pair(Parents, Edges);
  }
};
}  // namespace gbbs
