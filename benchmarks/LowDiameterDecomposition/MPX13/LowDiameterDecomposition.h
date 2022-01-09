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

#include "gbbs/gbbs.h"

#include <cmath>

namespace gbbs {
namespace ldd_utils {

inline size_t total_rounds(size_t n, double beta) {
  return std::min<uintE>(n + 1, 2 + ceil(log(n) / beta));
}

// Shifts[i] is the start of the vertices to take on round i
inline sequence<size_t> generate_shifts(size_t n, double beta) {
  // Create (ln n)/beta levels
  uintE last_round = total_rounds(n, beta);
  auto shifts = sequence<size_t>(last_round + 1);
  parallel_for(0, last_round, kDefaultGranularity,
               [&](size_t i) { shifts[i] = floor(exp(i * beta)); });
  shifts[last_round] = 0;
  parlay::scan_inplace(shifts);
  return shifts;
}

template <class Seq>
inline void num_clusters(Seq& s) {
  size_t n = s.size();
  auto flags =
      sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    if (!flags[s[i]]) {
      flags[s[i]] = 1;
    }
  });
  std::cout << "num. clusters = " << parlay::reduce(flags) << "\n";
}

template <class Seq>
inline void cluster_sizes(Seq& s) {
  size_t n = s.size();
  auto flags =
      sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    gbbs::write_add(&flags[s[i]], 1);
    //    if (!flags[s[i]]) {
    //      flags[s[i]] = 1;
    //    }
  });
  for (size_t i = 0; i < n; i++) {
    if (flags[i]) {
      std::cout << "Found cluster with size : " << flags[i] << std::endl;
    }
  }
}

template <class Graph, class Seq>
inline void num_intercluster_edges(Graph& G, Seq& s) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  auto ic_edges =
      sequence<size_t>::from_function(n, [&](size_t i) { return 0; });
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      return s[src] != s[ngh];
    };
    size_t ct = G.get_vertex(i).out_neighbors().count(pred);
    ic_edges[i] = ct;
  });
  std::cout << "num. intercluster edges = " << parlay::reduce(ic_edges) << "\n";
}
}  // namespace ldd_utils

template <class W, class EO>
struct LDD_F {
  uintE* cluster_ids;
  const EO& oracle;

  LDD_F(uintE* _cluster_ids, const EO& _oracle)
      : cluster_ids(_cluster_ids), oracle(_oracle) {}

  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    if (oracle(s, d, wgh)) {
      cluster_ids[d] = cluster_ids[s];
      return true;
    }
    return false;
  }

  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    if (oracle(s, d, wgh)) {
      return gbbs::atomic_compare_and_swap(&cluster_ids[d], UINT_E_MAX,
                                           cluster_ids[s]);
    }
    return false;
  }

  inline bool cond(uintE d) { return cluster_ids[d] == UINT_E_MAX; }
};

template <class Graph, class EO>
inline sequence<uintE> LDD_impl(Graph& G, const EO& oracle, double beta,
                                bool permute = true) {
  // Implementation based on "A Simple and Practical Linear-Work Parallel
  // Algorithm for Connectivity" by Shun, Dhulipala, and Blelloch, which is in
  // turn based on "Parallel Graph Decompositions Using Random Shifts" by
  // Miller, Peng, and Xu.
  timer gs;
  gs.start();
  using W = typename Graph::weight_type;
  size_t n = G.n;

  sequence<uintE> vertex_perm;
  if (permute) {
    vertex_perm = parlay::random_permutation<uintE>(n);
  }
  auto shifts = ldd_utils::generate_shifts(n, beta);
  gs.stop();
  debug(gs.next("generate shifts time"););
  auto cluster_ids = sequence<uintE>(n, UINT_E_MAX);

  timer add_t;
  timer vt;

  size_t round = 0, num_visited = 0;
  vertexSubset frontier(n);  // Initially empty
  size_t num_added = 0;
  while (num_visited < n) {
    size_t start = shifts[round];
    size_t end = std::min(static_cast<size_t>(shifts[round + 1]), n);
    size_t num_to_add = end - start;
    if (num_to_add > 0) {
      add_t.start();
      assert((num_added + num_to_add) <= n);
      auto candidates_f = [&](size_t i) {
        if (permute)
          return vertex_perm[num_added + i];
        else
          return static_cast<uintE>(num_added + i);
      };
      auto candidates = parlay::delayed_seq<uintE>(num_to_add, candidates_f);
      auto pred = [&](uintE v) { return cluster_ids[v] == UINT_E_MAX; };
      auto new_centers = parlay::filter(candidates, pred);
      add_to_vsubset(frontier, new_centers.begin(), new_centers.size());
      parallel_for(0, new_centers.size(), [&](size_t i) {
        uintE new_center = new_centers[i];
        cluster_ids[new_center] = new_center;
      });
      num_added += num_to_add;
      add_t.stop();
    }

    num_visited += frontier.size();
    if (num_visited >= n) break;

    vt.start();
    auto ldd_f = LDD_F<W, EO>(cluster_ids.begin(), oracle);
    frontier = edgeMap(G, frontier, ldd_f, -1, sparse_blocked);
    vt.stop();

    round++;
  }
  debug(add_t.next("add vertices time"); vt.next("edge map time"););
  return cluster_ids;
}

// Partitions vertices of an n-vertex, m-edge graph into subsets where each
// subset has O((log n) / beta) diameter and at most O(beta * m) edges exiting
// the subset.
//
// Returns:
//   Length n-length sequence `S` such that `S[i]` is the subset ID of vertex i.
//   The IDs are in the range [0, n) but are not necessarily contiguous.
template <class Graph>
sequence<uintE> LDD(Graph& G, double beta, bool permute = true) {
  using W = typename Graph::weight_type;
  debug(std::cout << "permute = " << permute << std::endl;);
  auto oracle = [&](const uintE& u, const uintE& v, const W& wgh) {
    return true;
  };
  return LDD_impl(G, oracle, beta, permute);
}

template <class Graph, class EO>
sequence<uintE> LDD_oracle(Graph& G, EO& oracle, double beta,
                           bool permute = true) {
  return LDD_impl(G, oracle, beta, permute);
}

}  // namespace gbbs
