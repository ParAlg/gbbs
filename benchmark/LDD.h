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

#include "pbbslib/random_shuffle.h"
#include "ligra.h"

#include <cmath>

namespace ldd_utils {
inline size_t total_rounds(size_t n, double beta) {
  return std::min<uintE>(n + 1, 2 + ceil(log(n) / beta));
}

// Shifts[i] is the start of the vertices to take on round i
inline sequence<uintE> generate_shifts(size_t n, double beta) {
  // Create (ln n)/beta levels
  uintE last_round = total_rounds(n, beta);
  auto shifts = sequence<uintE>(last_round + 1);
  par_for(0, last_round, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { shifts[i] = floor(exp(i * beta)); });
  shifts[last_round] = 0;
  pbbslib::scan_add_inplace(shifts);
  return shifts;
}

template <class Seq>
inline void num_clusters(Seq& s) {
  size_t n = s.size();
  auto flags = sequence<uintE>(n + 1, [&](size_t i) { return 0; });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    if (!flags[s[i]]) {
      flags[s[i]] = 1;
    }
  });
  std::cout << "num. clusters = " << pbbslib::reduce_add(flags) << "\n";
}

template <template <typename W> class vertex, class W, class Seq>
inline void num_intercluster_edges(graph<vertex<W> >& GA, Seq& s) {
  size_t n = GA.n;
  auto ic_edges = sequence<size_t>(n, [&](size_t i) { return 0; });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      return s[src] != s[ngh];
    };
    size_t ct = GA.V[i].countOutNgh(i, pred);
    ic_edges[i] = ct;
  });
  std::cout << "num. intercluster edges = " << pbbslib::reduce_add(ic_edges)
            << "\n";
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
      return CAS(&cluster_ids[d], UINT_E_MAX, cluster_ids[s]);
    }
    return false;
  }

  inline bool cond(uintE d) { return cluster_ids[d] == UINT_E_MAX; }
};

template <template <typename W> class vertex, class W, class EO>
inline sequence<uintE> LDD_impl(graph<vertex<W> >& GA, const EO& oracle,
                                  double beta, bool permute = true,
                                  bool pack = false) {
  size_t n = GA.n;

  sequence<uintE> vertex_perm;
  if (permute) {
    vertex_perm = pbbslib::random_permutation<uintE>(n);
  }
  auto shifts = ldd_utils::generate_shifts(n, beta);
  auto cluster_ids = sequence<uintE>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { cluster_ids[i] = UINT_E_MAX; });

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
      auto candidates = pbbslib::make_sequence<uintE>(num_to_add, candidates_f);
      auto pred = [&](uintE v) { return cluster_ids[v] == UINT_E_MAX; };
      auto new_centers = pbbslib::filter(candidates, pred);
      add_to_vsubset(frontier, new_centers.begin(), new_centers.size());
      par_for(0, new_centers.size(), pbbslib::kSequentialForThreshold, [&] (size_t i)
                      { cluster_ids[new_centers[i]] = new_centers[i]; });
      num_added += num_to_add;
    }

    num_visited += frontier.size();
    if (num_visited >= n) break;

    auto ldd_f = LDD_F<W, EO>(cluster_ids.begin(), oracle);
    vertexSubset next_frontier =
        edgeMap(GA, frontier, ldd_f, -1, sparse_blocked);
    if (pack) {
      auto pred = [&](const uintE& src, const uintE& dest, const W& w) {
        return oracle(src, dest, w) && (cluster_ids[src] != cluster_ids[dest]);
      };
//     timer t; t.start();
      edgeMapFilter(GA, frontier, pred, pack_edges | no_output);
//      t.stop(); t.reportTotal("pack time");
    }
    frontier.del();
    frontier = next_frontier;

    round++;
  }
  return cluster_ids;
}

template <template <typename W> class vertex, class W>
sequence<uintE> LDD(graph<vertex<W> >& GA, double beta, bool permute = true,
                      bool pack = false) {
  auto oracle = [&](const uintE& u, const uintE& v, const W& wgh) {
    return true;
  };
  return LDD_impl(GA, oracle, beta, permute, pack);
}

template <template <typename W> class vertex, class W, class EO>
sequence<uintE> LDD_oracle(graph<vertex<W> >& GA, EO& oracle, double beta,
                             bool permute = true, bool pack = false) {
  return LDD_impl(GA, oracle, beta, permute, pack);
}
