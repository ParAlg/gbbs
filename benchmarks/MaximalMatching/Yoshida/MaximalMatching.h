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

#include "gbbs/bridge.h"
#include "gbbs/gbbs.h"
#include "reorder.h"
#include "yoshida_matching_utils.h"

namespace gbbs {

enum mm_status { in, out, unknown };

/* returns pair of (in_matching, query_work) */
template <class Graph>
std::pair<mm_status, size_t> mm_query(uintE u, uintE v, Graph& G,
                                      size_t work_so_far, size_t query_cutoff) {
  assert(u < v);
  auto vtx_u = G.get_vertex(u);
  auto vtx_v = G.get_vertex(v);
  uintE deg_u = vtx_u.out_degree();
  uintE deg_v = vtx_v.out_degree();
  auto nghs_u = (uintE*)vtx_u.neighbors;
  auto nghs_v = (uintE*)vtx_v.neighbors;

  size_t our_pri = get_edge_pri(u, v);

  auto it = sorted_it(nghs_u, nghs_v, deg_u, deg_v, u, v);

  size_t work = work_so_far;
  [[maybe_unused]] size_t last_pri = 0;
  while (it.has_next()) {
    if (work == query_cutoff) {
      return std::make_pair(unknown, query_cutoff);
    }
    auto[a, n_a] = it.get_next();
    work++;

    auto ep = get_edge_pri(a, n_a);
    assert(ep >= last_pri); /* no inversions */
    last_pri = ep;
    if (our_pri < ep) {
      break; /* done */
    } else if (our_pri == ep) {
      if (a == u && n_a == v) {
        break;
      } /* not equal, but hashes equal */
      if (std::make_pair(u, v) < std::make_pair(a, n_a)) {
        break; /* done */
      }
    }

    /* recursively check (a, n_a) */
    auto[status, rec_work] = mm_query(a, n_a, G, work, query_cutoff);
    work = rec_work;
    if (status == in) { /* neighboring edge in mm */
      return std::make_pair(out, work);
    } else if (status == unknown) { /* hit query limit---return */
      return std::make_pair(unknown, work);
    }
  }
  /* made it! return in_matching */
  return std::make_pair(in, work);
}

// Runs a constant number of filters on the whole graph, and runs the
// prefix-based algorithm on them. Finishes off the rest of the graph with the
// prefix-based algorithm.
template <class Graph>
auto MaximalMatching(Graph& G, size_t query_cutoff) {
  using W = typename Graph::weight_type;

  size_t m = G.m;
  size_t n = G.n;

  auto edge_to_priority = [&](const uintE& u, const std::tuple<uintE, W>& ngh) {
    auto v = std::get<0>(ngh);
    return get_edge_pri(u, v);
  };
  auto RG = reorder_graph(G, edge_to_priority);

  size_t max_query_length = 0;
  auto matching_cts = sequence<uintE>(n, (uintE)0);
  auto total_work = sequence<size_t>(n, (size_t)0);
  auto got_answer = sequence<size_t>(n, (size_t)0);
  auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
    if (u < v) {
      auto[in_mm, work] = mm_query(u, v, RG, 0, query_cutoff);
      gbbs::write_add(&total_work[v], work);
      gbbs::write_max(&max_query_length, work, std::less<size_t>());
      if (in_mm == in) {
        gbbs::write_add(&matching_cts[u], 1);
        gbbs::write_add(&matching_cts[v], 1);
      }
      if (in_mm == in || in_mm == out) {
        gbbs::write_add(&got_answer[u], 1);
      }
    }
  };
  RG.mapEdges(map_f);
  size_t tot_work = parlay::reduce(make_slice(total_work));
  double fraction_covered =
      (static_cast<double>(parlay::reduce(make_slice(got_answer))) /
       (static_cast<double>(m) / 2));
  std::cout << "# Max query length = " << max_query_length << std::endl;
  std::cout << "# Total work = " << tot_work << std::endl;
  std::cout << "# Answers for: " << fraction_covered << " fraction of edges"
            << std::endl;

  if (query_cutoff == std::numeric_limits<size_t>::max()) {
    // Verify that we found a matching.
    for (size_t i = 0; i < n; i++) {
      if (matching_cts[i] > 1) {
        std::cout << "mct = " << matching_cts[i] << " i = " << i
                  << " deg = " << G.get_vertex(i).out_degree() << std::endl;
      }
      assert(matching_cts[i] <= 1);
    }

    // Verify that the matching is maximal.
    auto verify_f = [&](const uintE u, const uintE& v, const W& wgh) {
      assert(!(matching_cts[u] == 0 &&
               matching_cts[v] ==
                   0));  // !(we could add this edge to the matching)
    };
    G.mapEdges(verify_f);
  }

  return std::make_tuple(max_query_length, tot_work, fraction_covered);
}
}  // namespace gbbs
