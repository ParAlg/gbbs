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

#include "benchmarks/Connectivity/common.h"
#include "benchmarks/LowDiameterDecomposition/MPX13/LowDiameterDecomposition.h"
#include "gbbs/contract.h"
#include "gbbs/gbbs.h"

namespace gbbs {
namespace workefficient_cc {

template <class Graph>
inline sequence<parent> CC_impl(Graph& G, double beta, size_t level,
                                bool pack = false, bool permute = false) {
  size_t n = G.n;
  permute |= (level > 0);
  timer ldd_t;
  ldd_t.start();
  sequence<parent> clusters = LDD(G, beta, permute);
  debug(ldd_t.next("ldd time"););

  timer relabel_t;
  relabel_t.start();
  size_t num_clusters = contract::RelabelIds(clusters);
  debug(relabel_t.next("relabel time"););

  timer contract_t;
  contract_t.start();
  auto c_out = contract::contract(G, clusters, num_clusters);
  debug(contract_t.next("contract time"););
  // flags maps from clusters -> no-singleton-clusters
  auto GC = std::move(std::get<0>(c_out));
  auto& flags = std::get<1>(c_out);
  auto& mapping = std::get<2>(c_out);

  if (GC.m == 0) return clusters;

  auto new_labels = CC_impl(GC, beta, level + 1);
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    uintE cluster = clusters[i];
    uintE gc_cluster = flags[cluster];
    if (gc_cluster != flags[cluster + 1]) {  // was not a singleton
      // new_labels[gc_cluster] is the gc vertex that captured the whole
      // component. mapping maps this back to the original label range.
      clusters[i] = mapping[new_labels[gc_cluster]];
    }
  });
  new_labels.clear();
  flags.clear();
  mapping.clear();
  return clusters;
}

template <class Seq>
inline size_t num_cc(Seq& labels) {
  size_t n = labels.size();
  auto flags =
      sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    if (!flags[labels[i]]) {
      flags[labels[i]] = 1;
    }
  });
  parlay::scan_inplace(flags);
  std::cout << "# n_cc = " << flags[n] << "\n";
  return flags[n];
}

template <class Seq>
inline size_t largest_cc(Seq& labels) {
  size_t n = labels.size();
  // could histogram to do this in parallel.
  auto flags =
      sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    flags[labels[i]] += 1;
  }
  size_t sz = parlay::reduce_max(flags);
  std::cout << "# largest_cc has size: " << sz << "\n";
  return sz;
}

// Outputs a sequence `S` of length `G.n` such that the `i`-th vertex is in
// connected component `S[i]`. The component IDs will be in the range `[0, G.n)`
// but are not necessarily contiguous.
template <class Graph>
inline sequence<parent> CC(Graph& G, double beta = 0.2, bool pack = false,
                           bool permute = false) {
  return CC_impl(G, beta, 0, pack, permute);
}

}  // namespace workefficient_cc
}  // namespace gbbs
