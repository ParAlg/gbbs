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
#include "gbbs/helpers/speculative_for.h"

namespace gbbs {
namespace mm {

constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;
constexpr size_t n_filter_steps = 5;

template <class W>
struct matchStep {
  using edge = std::tuple<uintE, uintE, W>;
  edge* E;
  uintE* R;
  bool* matched;
  matchStep(edge* _E, uintE* _R, bool* m) : E(_E), R(_R), matched(m) {}

  bool reserve(uintE i) {
    uintE u = std::get<0>(E[i]);
    uintE v = std::get<1>(E[i]);
    if (matched[u] || matched[v] || (u == v)) return 0;
    reserveLoc<uintE>(&R[u], i);
    reserveLoc<uintE>(&R[v], i);
    return 1;
  }

  bool commit(uintE i) {
    uintE u = std::get<0>(E[i]);
    uintE v = std::get<1>(E[i]);
    if (R[v] == i) {
      R[v] = UINT_E_MAX;
      if (R[u] == i) {
        matched[u] = matched[v] = 1;
        R[u] = UINT_E_MAX;
        // mark edge
        E[i] = std::make_tuple(std::get<0>(E[i]) |= TOP_BIT, std::get<1>(E[i]),
                               std::get<2>(E[i]));
        return 1;
      }
    } else if (R[u] == i)
      R[u] = UINT_E_MAX;
    return 0;
  }
};

inline parlay::sequence<size_t> sample_indices(size_t n, size_t m, parlay::random r = random()) {
  auto out = parlay::tabulate(n, [&](size_t i) { return r.ith_rand(i) % m; });
  r = r.next();
  return out;
}

template <template <class W> class vertex, class W>
inline parlay::sequence<std::tuple<uintE, uintE, W>> get_edges(symmetric_graph<vertex, W>& G, size_t k,
                               bool* matched, size_t round, parlay::random r) {
  using edge = std::tuple<uintE, uintE, W>;

  timer tt;
  tt.start();
  auto live_degree = parlay::tabulate<size_t>(G.num_vertices(), [&](size_t i) {
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      return (round == 0 || !matched[v]) && u < v;
    };
    uintE degree = 0;
    if (!matched[i]) {
      auto monoid = parlay::plus<uintE>();
      degree = G.get_vertex(i).out_neighbors().reduce(map_f, monoid);
    }
    return degree;
  });
  gbbs_debug(tt.next("live degree time"););

  size_t tot = parlay::scan_inplace(parlay::make_slice(live_degree));
  gbbs_debug(std::cout << "tot = " << tot << std::endl;);

  parlay::sequence<size_t> indices;
  parlay::sequence<edge> out;

  if (tot < k) {
    indices = parlay::tabulate(tot, [&](size_t i) { return i; });
    out = parlay::sequence<edge>(tot);
    gbbs_debug(tt.next("permutation time"););
    k = tot;
  } else {
    indices = sample_indices(k, tot, r);
    gbbs_debug(tt.next("permutation time"););
    parlay::sort_inplace(indices);
    gbbs_debug(tt.next("sort time"););
    out = parlay::sequence<edge>(k);
  }

  constexpr size_t kBlockSize = 1 << 16;
  size_t num_blocks = parlay::num_blocks(k, kBlockSize);
  gbbs_debug(std::cout << "k = " << k << " num_blocks = " << num_blocks
                  << std::endl;);
  parlay::parallel_for(0, num_blocks, [&](size_t b) {
    size_t start = b * kBlockSize;
    size_t end = std::min(k, (b + 1) * kBlockSize);

    size_t i = start;
    size_t idx = indices[i];
    auto cur_id = parlay::binary_search(parlay::make_slice(live_degree), idx,
                                        std::less<size_t>());

    while (i < end) {
      idx = indices[i];
      size_t incr = parlay::internal::linear_search(
          parlay::make_slice(live_degree).cut(cur_id, live_degree.size()), idx,
          std::less<size_t>());
      cur_id += incr;
      if (live_degree[cur_id] > idx) cur_id--;

      size_t our_edge_start = live_degree[cur_id];
      size_t our_edge_end =
          ((idx == G.num_vertices() - 1) ? tot : live_degree[cur_id + 1]);
      size_t our_degree = our_edge_end - our_edge_start;
      while (our_degree == 0) {
        cur_id++;
        our_edge_start = live_degree[cur_id];
        our_edge_end =
            ((idx == G.num_vertices() - 1) ? tot : live_degree[cur_id + 1]);
        our_degree = our_edge_end - our_edge_start;
      }

      // detect first edge for vertex
      bool is_first_edge = (i == 0) || (indices[i - 1] < our_edge_start);
      if (is_first_edge) {
        size_t ctr = our_edge_start;
        size_t target = idx;
        size_t offset = i;
        auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
          if (!matched[v] && u < v) {
            if (target == ctr) {
              out[offset] = {u, v, wgh};
              offset++;
              if (offset < k) {
                target = indices[offset];
              } else {
                target = std::numeric_limits<size_t>::max();
              }
            }
            ctr++;
          }
        };
        G.get_vertex(cur_id).out_neighbors().map(map_f, false);

        i = offset;
      } else {
        i++;
      }
    }
  }, 1);

  gbbs_debug(tt.next("extract edges time"););

  return out;
}

};  // namespace mm

// Runs a constant number of filters on the whole graph, and runs the
// prefix-based algorithm on them. Finishes off the rest of the graph with the
// prefix-based algorithm.
template <template <class W> class vertex, class W>
inline sequence<std::tuple<uintE, uintE, W>> MaximalMatching(
    symmetric_graph<vertex, W>& G) {
  using edge = std::tuple<uintE, uintE, W>;

  timer mt;
  mt.start();
  size_t n = G.n;
  auto r = parlay::random();

  auto R =
      sequence<uintE>::from_function(n, [&](size_t i) { return UINT_E_MAX; });
  auto matched =
      sequence<bool>::from_function(n, [&](size_t i) { return false; });

  size_t k = std::max(((3 * G.n) / 2), 1024 * static_cast<size_t>(1024));
  auto matching = parlay::sequence<edge>();

  size_t round = 0;
  timer gete;
  timer eff;
  while (true) {
    gete.start();
    auto edges = mm::get_edges(G, k, matched.begin(), round, r);

    if (edges.size() == 0) break;
    gbbs_debug(gete.next("Get Edges Time"););

    gbbs_debug(std::cout << "Got: " << edges.size() << " edges " << std::endl;);

    mm::matchStep<W> mStep(edges.begin(), R.begin(), matched.begin());
    eff.start();
    eff_for<uintE>(mStep, 0, edges.size(), 50, 0, G.n);
    eff.stop();
    gbbs_debug(gete.next("Match Time"););

    auto e_added = parlay::filter(
        parlay::make_slice(edges), [](edge e) { return std::get<0>(e) & mm::TOP_BIT; });

    matching.append(parlay::make_slice(e_added));

    round++;
    r = r.next();
  }
  std::cout << "matching size = " << matching.size() << "\n";
  auto output = sequence<edge>::from_function(matching.size(), [&](size_t i) {
    uintE u, v;
    W wgh;
    std::tie(u, v, wgh) = matching[i];
    return edge(u & mm::VAL_MASK, v & mm::VAL_MASK, wgh);
  });  // allocated
  mt.stop();
  gbbs_debug(
  eff.next("eff for time");
  gete.next("get edges time");
  mt.next("Matching time"););
  return output;
}

template <template <class W> class vertex, class W, class Seq>
inline void verify_matching(symmetric_graph<vertex, W>& G, Seq& matching) {
  size_t n = G.n;
  auto ok = parlay::sequence<bool>::from_function(n, [](size_t i) { return 1; });
  auto matched = parlay::sequence<uintE>::from_function(n, [](size_t i) { return 0; });

  // Check that this is a valid matching
  parlay::parallel_for(0, matching.size(), [&](size_t i) {
    const auto& edge = matching[i];
    gbbs::write_add(&matched[std::get<0>(edge)], 1);
    gbbs::write_add(&matched[std::get<1>(edge)], 1);
  });

  bool valid = true;
  parlay::parallel_for(0, n, [&](size_t i) {
    if (matched[i] > 1) valid = false;
  });
  assert(valid == true);

  // Check maximality of the matching
  auto map2_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    if (!matched[src] && !matched[ngh]) {
      // could have added this edge, increasing the size of the matching
      ok[src] = 0;
      ok[ngh] = 0;
    }
  };
  parlay::parallel_for(0, n,
               [&](size_t i) { G.get_vertex(i).out_neighbors().map(map2_f); }, 1);

  auto ok_f = [&](size_t i) { return ok[i]; };
  auto ok_im = parlay::delayed_seq<size_t>(n, ok_f);
  size_t n_ok = parlay::reduce(ok_im);
  if (n == n_ok) {
    std::cout << "Matching OK! matching size is: " << matching.size() << "\n";
  } else {
    std::cout << "Matching invalid---" << (n - n_ok)
              << " vertices saw bad neighborhoods."
              << "\n";
  }
}

}  // namespace gbbs
