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

#include <cassert>
#include "ligra/ligra.h"
#include "ligra/speculative_for.h"
#include "ligra/union_find.h"
#include "ligra/pbbslib/dyn_arr.h"

#include "pbbslib/binary_search.h"
#include "pbbslib/random.h"
#include "pbbslib/sample_sort.h"

namespace MinimumSpanningForest_spec_for {
constexpr size_t sample_size = 10000;
inline size_t hash_to_range(size_t hsh, size_t range) { return hsh & range; }

inline size_t key_for_pair(uint32_t k1, uintE k2, pbbslib::random rnd) {
  size_t key = (static_cast<size_t>(k1) << 32) + static_cast<size_t>(k2);
  return rnd.ith_rand(key);
}

template <template <class W> class vertex, class W, class UF>
inline edge_array<W> get_remaining(symmetric_graph<vertex, W>& G, size_t k, UF& uf,
                                   pbbslib::random r) {
  auto filter_pred = [&](const uint32_t& src, const uintE& ngh, const W& wgh) {
    if (src < ngh) {
      return 2;  // return in array
    } else {
      return 1;  // filter, don't return in array
    }
    return 0;
  };
  return filter_edges(G, filter_pred);
}

template <template <class W> class vertex, class W, class UF>
inline void pack_shortcut_edges(symmetric_graph<vertex, W>& G, UF& uf) {
  auto filter_pred = [&](const uint32_t& src, const uintE& ngh,
                         const W& wgh) -> int {
    if (src > ngh) {
      return true;
    }
    auto c_src = uf.find(src);
    auto c_ngh = uf.find(ngh);
    return c_src == c_ngh;
  };
  std::cout << "Calling filter, G.m = " << G.m << "\n";
  filter_edges(G, filter_pred, no_output);
  std::cout << "G.m is now " << G.m << "\n";
}

template <template <class W> class vertex, class W, class UF>
inline edge_array<W> get_top_k(symmetric_graph<vertex, W>& G, size_t k, UF& uf,
                               pbbslib::random r, bool first_round = false) {
  if (k == static_cast<size_t>(G.m)) {
    return get_remaining(G, k, uf, r);
  }

  size_t m = (first_round) ? (G.m / 2) : G.m;
  size_t range = (1L << pbbslib::log2_up(G.m)) - 1;
  size_t scaled_size = sample_size * (((double)range) / m);

  // 1. Sample sample_size many edges and sort them.
  auto pred = [&](const uint32_t& src, const uintE& ngh, const W& wgh) -> bool {
    if (src < ngh) {
      size_t hash_val = hash_to_range(key_for_pair(src, ngh, r), range);
      return hash_val < scaled_size;
    }
    return 0;
  };
  auto sampled_e = sample_edges(G, pred);
  if (sampled_e.non_zeros == 0) {
    std::cout << "non_zeros = 0"
              << "\n";
    exit(0);
    return get_remaining(G, k, uf, r);
  }
  auto cmp_by_wgh = [](const std::tuple<uint32_t, uintE, intE>& left,
                       const std::tuple<uintE, uintE, intE>& right) {
    return std::get<2>(left) < std::get<2>(right);
  };
  pbbslib::sample_sort(pbbslib::make_sequence(sampled_e.E, sampled_e.non_zeros), cmp_by_wgh);

  // 2. Get approximate splitter.
  size_t ind = ((double)(k * sampled_e.non_zeros)) / G.m;
  auto splitter = sampled_e.E[ind];
  int32_t split_weight = std::get<2>(splitter);
  sampled_e.del();
  std::cout << "split wgh is: " << split_weight << "\n";

  // 3. Filter edges based on splitter
  auto filter_pred = [&](const uint32_t& src, const uintE& ngh, const W& wgh) {
    if (wgh <= split_weight) {
      if (src < ngh) {
        return 2;  // return in array
      } else {
        return 1;  // filter, but don't return in array
      }
    }
    return 0;
  };
  return filter_edges(G, filter_pred);
}

template <template <class W> class vertex, class W>
inline void MinimumSpanningForest(symmetric_graph<vertex, W>& GA) {
  using res = reservation<uintE>;
  using edge_t = std::tuple<uintE, uintE, W>;

  size_t n = GA.n;
  auto r = pbbslib::random();
  auto uf = UnionFind(n);

  auto mst_edges = pbbslib::dyn_arr<edge_t>(n);

  size_t iter = 0;
  while (GA.m > 0) {
    std::cout << "iter = " << iter << " m = " << GA.m << "\n";
    // 1. get weight of k'th smallest edge, and all edges smaller.
    size_t split_idx = std::min(n, (size_t)GA.m);
    timer get_t;
    get_t.start();
    auto edges = get_top_k(GA, split_idx, uf, r, (iter == 0));
    get_t.stop();
    get_t.reportTotal("get time");
    size_t n_edges = edges.non_zeros;
    auto cmp_by_wgh = [](const std::tuple<uint32_t, uintE, W>& left,
                         const std::tuple<uintE, uintE, W>& right) {
      return std::get<2>(left) < std::get<2>(right);
    };
    pbbslib::sample_sort(pbbslib::make_sequence(edges.E, n_edges), cmp_by_wgh);
    std::cout << "Prefix size = " << split_idx << " #edges = " << n_edges
              << " G.m is now = " << GA.m << "\n";

    // 2. initialize reservations, copy edge info, and run UF step.
    auto R = pbbslib::new_array_no_init<res>(n);
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { R[i] = res(); });
    sequence<bool> mstFlags =
        sequence<bool>(n_edges, [](size_t i) { return 0; });

    auto UFStep = make_uf_step<uintE>(edges, R, mstFlags, uf);
    speculative_for<uintE>(UFStep, 0, n_edges, 8);

    UFStep.clear();
    pbbslib::free_array(R);
    auto edge_imap_f = [&](size_t i) { return edges.E[i]; };
    auto edge_im =
        pbbslib::make_sequence<edge_t>(n_edges, edge_imap_f);
    auto edges_ret = pbbslib::pack(edge_im, mstFlags);
    std::cout << "added " << edges_ret.size() << "\n";
    mst_edges.copyIn(edges_ret, edges_ret.size());
    edges.del();
    mstFlags.clear();

    timer pack_t;
    pack_t.start();
    pack_shortcut_edges(GA, uf);
    pack_t.stop();
    pack_t.reportTotal("pack time");
    iter++;
  }
  std::cout << "n in mst: " << mst_edges.size << "\n";
  auto wgh_imap_f = [&](size_t i) { return std::get<2>(mst_edges.A[i]); };
  auto wgh_imap = pbbslib::make_sequence<size_t>(
      mst_edges.size, wgh_imap_f);
  std::cout << "wgh = " << pbbslib::reduce_add(wgh_imap) << "\n";

  mst_edges.del();
}

}  // namespace MinimumSpanningForest_spec_for
