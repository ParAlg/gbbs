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

inline size_t hash_to_range(size_t hsh, size_t range) { return hsh & range; }

inline size_t key_for_pair(uintE k1, uintE k2, parlay::random rnd) {
  size_t l = std::min(k1, k2);
  size_t r = std::max(k1, k2);
  size_t key = (l << 32) + r;
  return rnd.ith_rand(key);
}

template <template <class W> class vertex, class W>
inline edge_array<W> get_all_edges(symmetric_graph<vertex, W>& G, bool* matched,
                                   parlay::random rnd) {
  auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return !(matched[src] || matched[ngh]) && (src < ngh);
  };
  edge_array<W> E = filterAllEdges(G, pred);

  timer perm_t;
  perm_t.start();

  using edge = std::tuple<uintE, uintE, W>;
  sequence<edge> e_arr = std::move(E.E);

  auto perm = parlay::random_permutation<uintT>(e_arr.size());
  auto out = sequence<edge>(e_arr.size());
  parallel_for(0, e_arr.size(), kDefaultGranularity, [&](size_t i) {
    out[i] = e_arr[perm[i]];  // gather or scatter?
  });
  E.E = std::move(out);
  E.n = G.n;
  perm_t.stop();
  perm_t.next("permutation time");
  return E;
}

template <template <class W> class vertex, class W>
inline edge_array<W> get_edges(symmetric_graph<vertex, W>& G, size_t k,
                               bool* matched, parlay::random r) {
  using edge = std::tuple<uintE, uintE, W>;
  size_t m = G.m / 2;  // assume sym
  bool finish = (m <= k);

  std::cout << "Threshold, using m = " << m << "\n";
  size_t range = parlay::log2_up(G.m);
  range = 1L << range;
  range -= 1;

  size_t threshold = k;
  auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    size_t hash_val = hash_to_range(key_for_pair(src, ngh, r), range);
    if ((src > ngh) || matched[ngh]) {
      return 1;  // pack out, not in edgearr
    } else if (hash_val < threshold || finish) {
      return 2;  // pack out, returned in edgearr
    }
    return 0;  // keep in graph, not in edgearr
  };
  timer fet;
  fet.start();
  edge_array<W> E = filterEdges(G, pred);
  fet.stop();
  fet.next("Filter edges time");

  // permute the retrieved edges

  auto e_arr = E.to_seq();
  timer perm_t;
  perm_t.start();
  auto perm = parlay::random_permutation<uintT>(e_arr.size());
  auto out = sequence<edge>(e_arr.size());
  parallel_for(0, e_arr.size(), [&](size_t i) {
    out[i] = e_arr[perm[i]];  // gather or scatter?
  });
  E.E = std::move(out);
  E.n = G.n;
  perm_t.stop();
  perm_t.next("permutation time");
  return E;
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

  size_t k = ((3 * G.n) / 2);
  auto matching = parlay::sequence<edge>();

  size_t round = 0;
  timer gete;
  timer eff;
  while (G.m > 0) {
    gete.start();
    edge_array<W> e_arr = (round < mm::n_filter_steps)
                              ? mm::get_edges(G, k, matched.begin(), r)
                              : mm::get_all_edges(G, matched.begin(), r);

    auto eim_f = [&](size_t i) { return e_arr.E[i]; };
    auto eim = parlay::delayed_seq<edge>(e_arr.size(), eim_f);
    gete.stop();

    std::cout << "Got: " << e_arr.size() << " edges "
              << " G.m is now: " << G.m << "\n";
    mm::matchStep<W> mStep(e_arr.E.begin(), R.begin(), matched.begin());
    eff.start();
    eff_for<uintE>(mStep, 0, e_arr.size(), 50, 0, G.n);
    eff.stop();

    auto e_added = parlay::filter(
        eim, [](edge e) { return std::get<0>(e) & mm::TOP_BIT; });
    auto sizes = sequence<size_t>(e_added.size());
    parallel_for(0, e_added.size(), [&](size_t i) {
      const auto& e = e_added[i];
      uintE u = std::get<0>(e) & mm::VAL_MASK;
      uintE v = std::get<1>(e) & mm::VAL_MASK;
      uintE deg_u = G.get_vertex(u).out_degree();
      uintE deg_v = G.get_vertex(v).out_degree();
      G.zeroVertexDegree(u);
      G.zeroVertexDegree(v);
      sizes[i] = deg_u + deg_v;
    });
    size_t total_size = parlay::reduce(sizes);
    G.m -= total_size;
    std::cout << "removed: " << total_size << " many edges"
              << "\n";

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
  eff.next("eff for time");
  gete.next("get edges time");
  mt.next("Matching time");
  return output;
}

template <template <class W> class vertex, class W, class Seq>
inline void verify_matching(symmetric_graph<vertex, W>& G, Seq& matching) {
  size_t n = G.n;
  auto ok = sequence<bool>::from_function(n, [](size_t i) { return 1; });
  auto matched = sequence<uintE>::from_function(n, [](size_t i) { return 0; });

  // Check that this is a valid matching
  parallel_for(0, matching.size(), [&](size_t i) {
    const auto& edge = matching[i];
    gbbs::write_add(&matched[std::get<0>(edge)], 1);
    gbbs::write_add(&matched[std::get<1>(edge)], 1);
  });

  bool valid = true;
  parallel_for(0, n, [&](size_t i) {
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
  parallel_for(0, n, 1,
               [&](size_t i) { G.get_vertex(i).out_neighbors().map(map2_f); });

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
