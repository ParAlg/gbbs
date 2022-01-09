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
#include "reorder.h"

namespace gbbs {
enum mis_status { in, out, unknown };

size_t get_vertex_pri(uintE v) { return parlay::hash64(v); }

/* returns pair of (in_matching, query_work) */
template <class Graph>
std::pair<mis_status, size_t> mis_query(uintE u, Graph& G, size_t work_so_far,
                                        size_t query_cutoff) {
  auto vtx_u = G.get_vertex(u);
  uintE deg_u = vtx_u.out_degree();
  auto nghs_u = (uintE*)vtx_u.neighbors;
  size_t our_pri = get_vertex_pri(u);

  size_t work = work_so_far;
  size_t last_pri = 0;
  size_t ct_u = 0;
  while (ct_u < deg_u) {
    if (work == query_cutoff) {
      return std::make_pair(unknown, work);
    }
    uintE ngh = nghs_u[ct_u++];
    work++;

    auto pri = get_vertex_pri(ngh);
    assert(pri >= last_pri); /* no inversions */
    if (pri < last_pri) {
      exit(-1);
    }
    last_pri = pri;
    if (our_pri < pri) {
      break; /* done */
    } else if (our_pri == pri) {
      if (u <= ngh) {
        break; /* done */
      }
    }

    /* recursively check ngh */
    auto[status, rec_work] = mis_query(ngh, G, work, query_cutoff);
    work = rec_work;
    if (status == in) { /* neighboring edge in mm */
      return std::make_pair(out, work);
    } else if (status == unknown) {
      return std::make_pair(unknown, work);
    }
  }
  /* made it! return in_mis */
  return std::make_pair(in, work);
}

// Runs a constant number of filters on the whole graph, and runs the
// prefix-based algorithm on them. Finishes off the rest of the graph with the
// prefix-based algorithm.
template <class Graph>
auto MaximalIndependentSet(
    Graph& G, size_t query_cutoff = std::numeric_limits<size_t>::max()) {
  using W = typename Graph::weight_type;

  size_t n = G.n;

  auto RG = reorder_graph(G, get_vertex_pri);

  size_t max_query_length = 0;
  auto mis = sequence<bool>(n, false);
  auto answered = sequence<bool>(n, false);
  auto total_work = sequence<size_t>(n, (size_t)0);
  parallel_for(0, n, [&](size_t i) {
    auto[status, work] = mis_query(i, RG, 0, query_cutoff);
    mis[i] = (status == in);
    total_work[i] = work;
    answered[i] = (status != unknown);
    gbbs::write_max(&max_query_length, work, std::less<size_t>());
  });
  size_t tot_work = parlay::reduce(make_slice(total_work));
  auto answered_seq = parlay::delayed_seq<size_t>(
      n, [&](size_t i) { return static_cast<size_t>(answered[i]); });
  size_t num_answered = parlay::reduce(answered_seq);
  double fraction_covered =
      static_cast<double>(num_answered) / static_cast<double>(n);
  std::cout << "# Max query length = " << max_query_length << std::endl;
  std::cout << "# Total work = " << tot_work << std::endl;
  std::cout << "# Fraction covered = " << fraction_covered << std::endl;
  std::cout << "# Num answered = " << fraction_covered << std::endl;

  // Verify that the set is independent and maximal (if no cutoff)
  if (query_cutoff == std::numeric_limits<size_t>::max()) {
    parallel_for(0, n, [&](size_t i) {
      auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
        return (size_t)(mis[v]);
      };
      auto mon = parlay::addm<size_t>();
      size_t nghs_ct = G.get_vertex(i).out_neighbors().reduce(map_f, mon);
      if (mis[i]) {  // if in, ensure no neighbors are in
        assert(nghs_ct == 0);
        if (nghs_ct != 0) {
          exit(-1);
        }
      } else {  // if out, ensure some neighbor is in
        assert(nghs_ct > 0);
        if (nghs_ct == 0) {
          exit(-1);
        }
      }
    });
  }
  return std::make_tuple(max_query_length, tot_work, fraction_covered);
}
}  // namespace gbbs
