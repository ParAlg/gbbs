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

#include "lib/index_map.h"
#include "lib/random_shuffle.h"
#include "ligra.h"

#include <cmath>

namespace ldd_utils {
size_t total_rounds(size_t n, double beta) {
  return min<uintE>(n + 1, 2 + ceil(log(n) / beta));
}

// Shifts[i] is the start of the vertices to take on round i
auto generate_shifts(size_t n, double beta) {
  // Create (ln n)/beta levels
  uintE last_round = total_rounds(n, beta);
  auto shifts = array_imap<uintE>(last_round + 1);
  parallel_for_bc(i, 0, last_round,
                  (last_round > pbbs::kSequentialForThreshold),
                  { shifts[i] = floor(exp(i * beta)); });
  shifts[last_round] = 0;
  pbbs::scan_add(shifts, shifts);
  return shifts;
}

template <class Seq>
void num_clusters(Seq& s) {
  size_t n = s.size();
  auto flags = array_imap<uintE>(n + 1, [&](size_t i) { return 0; });
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    if (!flags[s[i]]) {
      flags[s[i]] = 1;
    }
  });
  size_t n_clusters = pbbs::reduce_add(flags);
  cout << "num. clusters = " << pbbs::reduce_add(flags) << endl;
}

template <template <typename W> class vertex, class W, class Seq>
void num_intercluster_edges(graph<vertex<W> >& GA, Seq& s) {
  size_t n = GA.n;
  auto ic_edges = array_imap<size_t>(n, [&](size_t i) { return 0; });
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      return s[src] != s[ngh];
    };
    size_t ct = GA.V[i].countOutNgh(i, pred);
    ic_edges[i] = ct;
  });
  cout << "num. intercluster edges = " << pbbs::reduce_add(ic_edges) << endl;
}
}  // namespace ldd_utils

template <class W>
struct LDD_F {
  uintE* cluster_ids;

  LDD_F(uintE* _cluster_ids) : cluster_ids(_cluster_ids) {}

  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    cluster_ids[d] = cluster_ids[s];
    return true;
  }

  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    return CAS(&cluster_ids[d], UINT_E_MAX, cluster_ids[s]);
  }

  inline bool cond(uintE d) { return cluster_ids[d] == UINT_E_MAX; }
};

template <template <typename W> class vertex, class W>
auto LDD(graph<vertex<W> >& GA, double beta, bool permute = true,
         bool pack = false) {
  using w_vertex = vertex<W>;
  size_t n = GA.n;
  size_t m = GA.m;
  auto V = GA.V;
  timer lddt;
  lddt.start();

  uintE* vertex_perm;
  if (permute) {
    auto perm = pbbs::random_permutation<uintE>(n);
    vertex_perm = perm.get_array();
  }
  auto shifts = ldd_utils::generate_shifts(n, beta);
  auto cluster_ids = array_imap<uintE>(n);
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                  { cluster_ids[i] = UINT_E_MAX; });

  size_t last_round = ldd_utils::total_rounds(n, beta);
  size_t round = 0, num_visited = 0;
  vertexSubset frontier(n);  // Initially empty
  size_t num_added = 0;
  cout << "last_round = " << last_round << endl;
  while (num_visited < n) {
    assert(round < last_round);
    size_t start = shifts[round];
    size_t end = min(static_cast<size_t>(shifts[round + 1]), n);
    size_t num_to_add = end - start;
    if (num_to_add > 0) {
      assert((num_added + num_to_add) <= n);
      auto candidates = make_in_imap<uintE>(num_to_add, [&](size_t i) {
        if (permute)
          return vertex_perm[num_added + i];
        else
          return static_cast<uintE>(num_added + i);
      });
      auto pred = [&](uintE v) { return cluster_ids[v] == UINT_E_MAX; };
      auto new_centers = pbbs::filter(candidates, pred);
      add_to_vsubset(frontier, new_centers.start(), new_centers.size());
      parallel_for_bc(i, 0, new_centers.size(), (new_centers.size() > pbbs::kSequentialForThreshold),
                      { cluster_ids[new_centers[i]] = new_centers[i]; });
      num_added += num_to_add;
    }

    num_visited += frontier.size();
    if (num_visited >= n) break;

    auto ldd_f = LDD_F<W>(cluster_ids.start());
    vertexSubset next_frontier =
        edgeMap(GA, frontier, ldd_f, -1, sparse_blocked);
    if (pack) {
      auto pred = [&](const uintE& src, const uintE& dest, const W& w) {
        return cluster_ids[src] != cluster_ids[dest];
      };
      edgeMapFilter(GA, frontier, pred, pack_edges | no_output);
    }
    frontier.del();
    frontier = next_frontier;

    round++;
  }
  double tt = lddt.stop();
  return std::move(cluster_ids);
}
