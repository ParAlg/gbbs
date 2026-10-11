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

#include <cstddef>
#include <cstdint>  // Added for uint64_t
#include <vector>

#include "absl/types/span.h"
#include "gbbs/bridge.h"
#include "gbbs/edge_map_data.h"
#include "gbbs/flags.h"
#include "gbbs/macros.h"
#include "gbbs/vertex_subset.h"
#include "parlay/random.h"

namespace gbbs {

// =============================================================================
// Parallel Graph Coloring using Jones-Plassmann Algorithm
// =============================================================================
//
// This file implements parallel graph coloring using the Jones-Plassmann (JP)
// algorithm with two ordering heuristics: LF (Largest-degree-First) and LLF
// (Largest-Log-degree-First).
//
// Reference:
//   Hasenplaugh, W., Kaler, T., Schardl, T. B., & Leiserson, C. E. (2014).
//   "Ordering heuristics for parallel graph coloring."
//   Proceedings of the 26th ACM Symposium on Parallelism in Algorithms and
//   Architectures (SPAA '14), pp. 166-177.
//   https://dl.acm.org/doi/10.1145/2612669.2612697
//
// -----------------------------------------------------------------------------
// Algorithm Overview
// -----------------------------------------------------------------------------
//
// The Jones-Plassmann algorithm is a parallel greedy graph coloring algorithm.
// It works by assigning each vertex a unique priority, then iteratively
// coloring vertices that have higher priority than all their uncolored
// neighbors.
//
// In each round:
//   1. Identify "roots" - uncolored vertices with priority higher than all
//      uncolored neighbors. Roots form an independent set (no two roots are
//      adjacent), since if both u and v are neighbors, at most one can have
//      higher priority than the other.
//   2. Color each root with the smallest color not used by its neighbors.
//      Since roots are independent, they can be colored in parallel without
//      conflicting with each other.
//
// -----------------------------------------------------------------------------
// LF vs LLF Heuristics
// -----------------------------------------------------------------------------
//
// The choice of priority ordering affects both the span (depth) of the overall
// algorithm and the quality of the coloring. Below we give some high-level
// intuition for the two ordering heuristics implemented in this file:
//
// LF (Largest-degree-First):
//   - A vertex v has priority over neighbor u if degree(v) > degree(u), with
//     ties broken by a random permutation (`rank`): among vertices of equal
//     degree, the one with smaller rank has higher priority.
//   - Produces high-quality colorings (fewer colors) in practice, but is
//     vulnerable to adversarial inputs that limit parallelism. The LLF
//     heuristic, described next, avoids this issue while only using marginally
//     more colors in practice.
//
// LLF (Largest-Log-degree-First):
//   - A vertex v has priority over neighbor u if round_up(log2(1 +
//     degree(v))) > round_up(log2(1 + degree(u))), with ties broken by a
//     random permutation (`rank`) as in LF. The rounding groups vertices of
//     similar degree into the same priority class, increasing the number of
//     ties broken by the random permutation and providing better parallelism
//     guarantees.
//   - Has provable bounds on the span of the coloring algorithm (see below).
//   - Produces colorings of competitive quality to LF in practice.
//
// Theoretical guarantees proven in Hasenplaugh et al.:
//   - JP-LLF colors a graph G = (V, E) with degree Δ using:
//       * Work: Θ(|V| + |E|)
//       * Span: O(log |V| + log Δ · min{sqrt(|E|), Δ + log Δ · log |V| / log
//       log |V|}).
//   The span term is usually dominated by the sqrt(|E|) term, so this algorithm
//   provably gets polynomial speedup (the speedup is >= sqrt(|E|) when |E| >=
//   |V|).
//
// -----------------------------------------------------------------------------
// Input Requirements
// -----------------------------------------------------------------------------
//
// The input graph MUST be symmetric (undirected). If the graph is asymmetric
// (directed), the algorithm may fail to terminate. This happens because a
// vertex v becomes a "root" (and gets colored) only when its priority counter
// reaches 0, which requires receiving a decrement from each higher-priority
// neighbor when that neighbor is colored. If edge (u, v) exists but (v, u)
// does not, then u will never send a decrement to v, so v's counter may never
// reach 0, causing the algorithm to hang.
//
// The algorithm does NOT validate that the input graph is symmetric for
// performance reasons.
//
// Self-loops, if present, are silently ignored.
//
// -----------------------------------------------------------------------------
// Output Properties
// -----------------------------------------------------------------------------
//
// **Color Range**: Colors are non-negative integers starting from 0.
//   - The minimum color assigned is always 0 (assuming the graph is non-empty)
//   - Colors are assigned greedily as the smallest available color
//
// **Color Contiguity**: The set of colors used is always {0, 1, ..., k} for
//   some k. This follows from the greedy rule: if vertex v receives color c,
//   then colors 0 through c-1 must each appear on at least one of v's
//   already-colored neighbors.
//
// **Number of Colors**: The algorithm uses at most Δ + 1 colors, where Δ is
// the maximum degree. In practice, it often uses far fewer colors.
//
// -----------------------------------------------------------------------------
// Determinism and Seed Stability
// -----------------------------------------------------------------------------
//
// The algorithm supports deterministic execution via the `ColoringWithSeed`
// function, which accepts a `coloring_seed` parameter.
//
// **Seed Stability Guarantees:**
//
// Within a given process, if the same input graph, `use_lf` setting, and
// `coloring_seed` are used, then the same coloring will be produced.
//
// The default `Coloring` function uses a fixed default seed (0), so it will
// produce deterministic results when called with the same graph and settings.
//
// =============================================================================

namespace coloring {

// Assigns the smallest color not used by any neighbor of vertex v.
// Uses thread-local storage for efficiency in parallel execution.
template <class Graph>
inline uintE AssignSmallestAvailableColor(Graph& G, uintE v,
                                          sequence<uintE>& colors) {
  using W = typename Graph::weight_type;
  uintE deg = G.get_vertex(v).out_degree();
  if (deg == 0) {
    return 0;
  }

  thread_local std::vector<bool> neighbor_colors_used;
  neighbor_colors_used.resize(deg);
  neighbor_colors_used.assign(deg, false);

  auto mark_neighbor_color = [&](uintE src, uintE ngh, const W& wgh) {
    if (src == ngh) return;
    uintE neighbor_color = colors[ngh];
    if (neighbor_color < deg) {
      neighbor_colors_used[neighbor_color] = true;
    }
  };
  G.get_vertex(v).out_neighbors().map(mark_neighbor_color, /*parallel=*/false);

  uintE smallest_available = UINT_E_MAX;
  for (size_t i = 0; i < deg; ++i) {
    if (!neighbor_colors_used[i]) {
      smallest_available = static_cast<uintE>(i);
      break;
    }
  }

  return (smallest_available == UINT_E_MAX) ? deg : smallest_available;
}

// Edge function for propagating coloring updates via edgeMap.
// Decrements the count of higher-priority uncolored neighbors and returns
// true when a vertex becomes a root (all higher-priority neighbors colored).
template <class W>
struct ColoringEdgeFunction {
  absl::Span<intE> remaining_higher_priority_neighbors;

  explicit ColoringEdgeFunction(absl::Span<intE> counts)
      : remaining_higher_priority_neighbors(counts) {}

  inline bool update(const uintE& src, const uintE& dst, const W& w) {
    if (src == dst) return false;
    --remaining_higher_priority_neighbors[dst];
    return remaining_higher_priority_neighbors[dst] == 0;
  }

  inline bool updateAtomic(const uintE& src, const uintE& dst, const W& wgh) {
    if (src == dst) return false;
    return (gbbs::fetch_and_add(&remaining_higher_priority_neighbors[dst],
                                -1) == 1);
  }

  inline bool cond(uintE d) {
    return remaining_higher_priority_neighbors[d] > 0;
  }
};

// Computes the number of higher-priority neighbors for each vertex.
// Priority is determined by degree (LF/LLF heuristic) with rank as tiebreaker.
template <class Graph>
inline sequence<intE> ComputeVertexPrioritiesWithRanks(
    Graph& G, bool use_lf, absl::Span<const uintE> ranks) {
  using W = typename Graph::weight_type;
  const size_t num_vertices = G.n;
  auto priorities = sequence<intE>(num_vertices);

  auto ordering_keys =
      sequence<uintE>::from_function(num_vertices, [&](size_t v) -> uintE {
        if (use_lf) {
          return G.get_vertex(v).out_degree();
        } else {
          return parlay::log2_up(1 + G.get_vertex(v).out_degree());
        }
      });

  parallel_for(0, num_vertices, /*granularity=*/1, [&](size_t v) {
    auto has_higher_priority = [&](uintE /*src*/, uintE ngh, const W& wgh) {
      if (ngh == v) return false;
      return (ordering_keys[ngh] > ordering_keys[v]) ||
             ((ordering_keys[ngh] == ordering_keys[v]) &&
              ranks[ngh] < ranks[v]);
    };

    priorities[v] = G.get_vertex(v).out_neighbors().count(has_higher_priority);
  });

  return priorities;
}

template <class Graph>
inline sequence<intE> ComputeVertexPriorities(Graph& G, bool use_lf) {
  const size_t num_vertices = G.n;
  auto random_ranks = parlay::random_permutation<uintE>(num_vertices);
  return ComputeVertexPrioritiesWithRanks(G, use_lf, random_ranks);
}

template <class Graph>
inline void ColorVertexSet(Graph& G, vertexSubset& vertices,
                           sequence<uintE>& colors) {
  vertices.toSparse();
  parallel_for(0, vertices.size(), /*granularity=*/1, [&](size_t i) {
    uintE v = vertices.vtx(i);
    colors[v] = coloring::AssignSmallestAvailableColor(G, v, colors);
  });
}

// Computes a valid vertex coloring of the input graph using the provided
// ranks for tie-breaking. This overload is useful for deterministic testing.
//
// Args:
//   G: A symmetric (undirected) graph.
//   use_lf: If true, use LF heuristic. If false, use LLF heuristic.
//   ranks: A span of length num_vertices where ranks[v] is the tie-breaking
//          rank for vertex v. Lower ranks have higher priority when degrees
//          are equal.
//
// Returns:
//   A sequence of length num_vertices where the i-th element is the color
//   assigned to vertex i. Colors are non-negative integers starting from 0.
template <class Graph>
inline sequence<uintE> ColoringWithRanks(Graph& G, bool use_lf,
                                         absl::Span<const uintE> ranks) {
  using W = typename Graph::weight_type;
  const size_t num_vertices = G.n;

  auto priorities =
      coloring::ComputeVertexPrioritiesWithRanks(G, use_lf, ranks);

  auto colors = sequence<uintE>::from_function(
      num_vertices, [](size_t i) { return UINT_E_MAX; });

  auto is_root = [&](size_t i) { return priorities[i] == 0; };
  auto root_flags = parlay::delayed_seq<bool>(num_vertices, is_root);
  auto roots =
      vertexSubset(num_vertices, parlay::pack_index<uintE>(root_flags));

  size_t num_colored = 0;
  while (num_colored != num_vertices) {
    num_colored += roots.size();

    coloring::ColorVertexSet(G, roots, colors);

    auto next_roots =
        edgeMap(G, roots, ColoringEdgeFunction<W>(absl::MakeSpan(priorities)),
                -1, sparse_blocked);
    roots = std::move(next_roots);
  }

  return colors;
}

// Computes a valid vertex coloring of the input graph using a specified seed
// for the random number generator.
//
// Args:
//   G: A symmetric (undirected) graph. See documentation above for
//      requirements.
//   use_lf: If true, use LF (Largest-degree-First) heuristic.
//           If false, use LLF (Largest-Log-degree-First) heuristic.
//   coloring_seed: The seed for the random number generator used for
//                  tie-breaking. See the "Determinism and Seed Stability"
//                  section above for stability guarantees.
//
// Returns:
//   A sequence of length num_vertices where the i-th element is the color
//   assigned to vertex i. Colors are non-negative integers starting from 0.
template <class Graph>
inline sequence<uintE> ColoringWithSeed(Graph& G, bool use_lf,
                                        int32_t coloring_seed) {
  const size_t num_vertices = G.n;
  auto random_ranks = parlay::random_permutation<uintE>(
      num_vertices, parlay::random(static_cast<size_t>(coloring_seed)));
  return ColoringWithRanks(G, use_lf, random_ranks);
}

// Counts the number of coloring violations, i.e., vertices that have at least
// one neighbor with the same color. Returns 0 if the coloring is valid.
template <class Graph>
inline size_t CountColoringViolations(Graph& G, sequence<uintE>& colors) {
  using W = typename Graph::weight_type;
  size_t num_vertices = G.n;
  auto has_violation = sequence<bool>(num_vertices);

  parallel_for(0, num_vertices, [&](size_t i) {
    uintE src_color = colors[i];
    auto has_same_color_neighbor = [&](const uintE& src, const uintE& ngh,
                                       const W& wgh) {
      if (src == ngh) return false;
      return src_color == colors[ngh];
    };
    size_t violation_count =
        G.get_vertex(i).out_neighbors().count(has_same_color_neighbor);
    has_violation[i] = (violation_count > 0);
  });

  auto violation_counts = parlay::delayed_seq<size_t>(
      num_vertices,
      [&](size_t i) { return static_cast<size_t>(has_violation[i]); });
  return parlay::reduce(violation_counts);
}

}  // namespace coloring

// Computes a valid vertex coloring of the input graph.
//
// This function uses a default seed of 0 for the random number generator,
// providing deterministic results for the same input graph and settings.
// Use `ColoringWithSeed` to specify a custom seed.
//
// Args:
//   G: A symmetric (undirected) graph. See documentation above for
//      requirements.
//   use_lf: If true, use LF (Largest-degree-First) heuristic.
//           If false (default), use LLF (Largest-Log-degree-First) heuristic.
//
// Returns:
//   A sequence of length n where the i-th element is the color assigned to
//   vertex i. Colors are non-negative integers starting from 0.
template <class Graph>
inline sequence<uintE> Coloring(Graph& G, bool use_lf = false) {
  constexpr int32_t kDefaultColoringSeed = 0;
  return coloring::ColoringWithSeed(G, use_lf, kDefaultColoringSeed);
}

}  // namespace gbbs
