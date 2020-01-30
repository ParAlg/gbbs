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

#include "ligra/bridge.h"
#include "ligra/ligra.h"
#include "ligra/speculative_for.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "pbbslib/random_shuffle.h"
#include "reorder.h"

size_t n = 0;

size_t get_edge_pri(uintE u, uintE v) {
  uintE min_u = std::min(u, v);
  uintE max_u = std::max(u, v);
  size_t key = (static_cast<size_t>(min_u) << 32UL) | static_cast<size_t>(max_u);
  return pbbs::hash64(key);
}

template <class LT>
struct sorted_it {
  uintE* ngh_a;
  uintE* ngh_b;
  size_t deg_a;
  size_t deg_b;
  uintE id_a;
  uintE id_b;
  size_t ct_a;
  size_t ct_b;
  size_t tot;
  LT& lt;
  sorted_it(uintE* ngh_a, uintE* ngh_b, uintE deg_a, uintE deg_b, uintE id_a, uintE id_b, LT& lt) : ngh_a(ngh_a), ngh_b(ngh_b), deg_a(deg_a), deg_b(deg_b), id_a(id_a), id_b(id_b), lt(lt) {
    ct_a = 0;
    ct_b = 0;
    tot = deg_a + deg_b;
  }

  std::pair<uintE, uintE> get_next() {
    if (ct_a < deg_a) {
      if (ct_b < deg_b) {
        uintE a = ngh_a[ct_a];
        uintE b = ngh_b[ct_b];
        auto pair_a = std::make_pair(std::min(id_a, a), std::max(id_a, a));
        auto pair_b = std::make_pair(std::min(id_b, b), std::max(id_b, b));
        auto ret = lt(pair_a.first, pair_a.second, pair_b.first, pair_b.second);
        if (ret == pair_a) {
          ct_a++;
        } else {
          ct_b++;
        }
        assert(ret.first < n && ret.second < n);
        return ret;
      } else { /* exhausted b */
        assert(ct_a < deg_a);
        assert(id_a < n && ngh_a[ct_a] < n);
        return std::make_pair(id_a, ngh_a[ct_a++]);
      }
    } else { /* exhausted a */
      assert(ct_b < deg_b);
      assert(id_b < n && ngh_b[ct_b] < n);
      return std::make_pair(id_b, ngh_b[ct_b++]);
    }
  }

  bool has_next() {
    return ct_a + ct_b < tot;
  }
};

/* returns pair of (in_matching, query_work) */
template <class Graph>
std::pair<bool, size_t> mm_query(uintE u, uintE v, Graph& G) {
  assert(u < v);
  auto vtx_u = G.get_vertex(u); auto vtx_v = G.get_vertex(v);
  uintE deg_u = vtx_u.getOutDegree(); uintE deg_v = vtx_v.getOutDegree();
  auto nghs_u = (uintE*)vtx_u.getOutNeighbors(); auto nghs_v = (uintE*)vtx_v.getOutNeighbors();

  size_t our_pri = get_edge_pri(u,v);

  auto lt = [&] (const uintE& aa, const uintE& n_aa, const uintE& bb, const uintE& n_bb) {
    uintE a = std::min(aa, n_aa); uintE n_a = std::max(aa, n_aa);
    uintE b = std::min(bb, n_bb); uintE n_b = std::max(bb, n_bb);

    auto a_p = get_edge_pri(a,n_a);
    auto b_p = get_edge_pri(b,n_b);
    if (a_p < b_p) {
      return std::make_pair(a, n_a);
    } else if (b_p < a_p) {
      return std::make_pair(b, n_b);
    } else {
      assert(a != n_a && b != n_b);
      if (std::make_pair(a, n_a) < std::make_pair(b, n_b)) {
        return std::make_pair(a, n_a);
      } else {
        return std::make_pair(b, n_b);
      }
    }
  };
  auto it = sorted_it(nghs_u, nghs_v, deg_u, deg_v, u, v, lt);

  size_t work = 0;
  size_t last_pri = 0;
  while (it.has_next()) {
    auto [aa, n_aa] = it.get_next();
    auto a = std::min(aa, n_aa); auto n_a = std::max(aa, n_aa);
    work++;

    auto e_p = get_edge_pri(a, n_a);
    assert(e_p >= last_pri); /* no inversions */
    last_pri = e_p;
    if (our_pri < e_p) {
      break; /* done */
    } else if (our_pri == e_p) {
      if (a == u && n_a == v) {
        break;
      } /* not equal, but hashes equal */
      if (std::make_pair(u,v) < std::make_pair(a, n_a)) {
        break; /* done */
      }
    }
    /* have to recursively check (a, n_a) */

    auto [status, rec_work] = mm_query(a, n_a, G);
    work += rec_work;
    if (status) {
      return std::make_pair(false, work);
    } /* else: continue */
  }
  /* made it! return in_matching */
  return std::make_pair(true, work);
}

// Runs a constant number of filters on the whole graph, and runs the
// prefix-based algorithm on them. Finishes off the rest of the graph with the
// prefix-based algorithm.
template <class Graph>
auto MaximalMatching(Graph& G) {
  using W = typename Graph::weight_type;
  using edge = std::tuple<uintE, uintE, W>;

  size_t m = G.m;
  n = G.n;

  auto edge_to_priority = [&] (const uintE& u, const std::tuple<uintE, W>& ngh) {
    auto v = std::get<0>(ngh);
    return get_edge_pri(u, v);
  };
  auto RG = reorder_graph(G, edge_to_priority);

  size_t max_query_length = 0;
  auto matching_cts = pbbs::sequence<uintE>(n, (uintE)0);
  auto total_work = pbbs::sequence<size_t>(n, (size_t)0);
  auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
    if (u < v) {
      auto [in_mm, work] = mm_query(u, v, RG);
      pbbs::write_add(&total_work[v], work);
      pbbs::write_max(&max_query_length, work, std::less<size_t>());
      if (in_mm) {
        pbbs::write_add(&matching_cts[u], 1);
        pbbs::write_add(&matching_cts[v], 1);
      }
    }
  };
  RG.map_edges(map_f);
  std::cout << "Max query length = " << max_query_length << std::endl;
  std::cout << "Total work = " << pbbslib::reduce_add(total_work.slice()) << std::endl;

  // Verify that we found a matching.
  for (size_t i=0; i<n; i++) {
    if (matching_cts[i] > 1) {
      std::cout << "mct = " << matching_cts[i] << " i = " << i << " deg = " << G.get_vertex(i).getOutDegree() << std::endl;
    }
    assert(matching_cts[i] <= 1);
  }

  // Verify that the matching is maximal.
  auto verify_f = [&] (const uintE u, const uintE& v, const W& wgh) {
    assert(!(matching_cts[u] == 0 && matching_cts[v] == 0)); // !(we could add this edge to the matching)
  };
  G.map_edges(verify_f);

  return 1;
}
