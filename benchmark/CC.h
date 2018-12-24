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

#include "LDD.h"
#include "lib/sparse_table.h"
#include "ligra.h"

namespace cc {

// Relabel ids, a length n array with values in the range [0,n), to be in the
// range [0..u-1] where u is the number of unique ids
// Returns u
template <class Seq>
inline size_t RelabelIds(Seq& ids) {
  using T = typename Seq::T;
  size_t n = ids.size();
  auto inverse_map = sequence<T>(n + 1);
  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i)
                  { inverse_map[i] = 0; });
  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i) {
    if (!inverse_map[ids[i]]) inverse_map[ids[i]] = 1;
  });
  pbbs::scan_add(inverse_map, inverse_map);

  size_t new_n = inverse_map[n];
  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i)
                  { ids[i] = inverse_map[ids[i]]; });
  return new_n;
}

template <template <typename W> class vertex, class W, class Seq, class EO>
inline std::tuple<graph<symmetricVertex<pbbs::empty>>, sequence<uintE>,
                  sequence<uintE>>
contract(graph<vertex<W>>& GA, Seq& clusters, size_t num_clusters, EO& oracle) {
  // Remove duplicates by hashing
  using K = std::tuple<uintE, uintE>;
  using V = pbbs::empty;
  using KV = std::tuple<K, V>;

  size_t n = GA.n;

  timer count_t;
  count_t.start();
  auto deg_map = sequence<uintE>(n + 1);
  auto pred = [&](const uintE& src, const uintE& ngh, const W& w) {
    if (oracle(src, ngh, w)) {
      uintE c_src = clusters[src];
      uintE c_ngh = clusters[ngh];
      return c_src < c_ngh;
    }
    return false;
  };
  parallel_for_bc(i, 0, n, true,
                  { deg_map[i] = GA.V[i].countOutNgh(i, pred); });
  deg_map[n] = 0;
  pbbs::scan_add(deg_map, deg_map);
  count_t.stop();
  count_t.reportTotal("count time");

  timer ins_t;
  ins_t.start();
  KV empty =
      std::make_tuple(std::make_tuple(UINT_E_MAX, UINT_E_MAX), pbbs::empty());
  auto hash_pair = [](const std::tuple<uintE, uintE>& t) {
    return pbbs::hash64(std::get<0>(t)) ^ pbbs::hash64(std::get<1>(t));
  };
  auto edge_table = make_sparse_table<K, V>(deg_map[n], empty, hash_pair);
  deg_map.del();

  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
    if (oracle(src, ngh, w)) {
      uintE c_src = clusters[src];
      uintE c_ngh = clusters[ngh];
      if (c_src < c_ngh) {
        edge_table.insert(
            std::make_tuple(std::make_tuple(c_src, c_ngh), pbbs::empty()));
      }
    }
  };
  parallel_for_bc(i, 0, n, true, { GA.V[i].mapOutNgh(i, map_f); });
  auto edges = edge_table.entries();
  edge_table.del();
  ins_t.stop();
  ins_t.reportTotal("ins time");

  // Pack out singleton clusters
  auto flags = sequence<uintE>(num_clusters + 1, [](size_t i) { return 0; });

  parallel_for_bc(i, 0, edges.size(),
                  (edges.size() > pbbs::kSequentialForThreshold), {
                    auto e = std::get<0>(edges[i]);
                    uintE u = std::get<0>(e);
                    uintE v = std::get<1>(e);
                    if (!flags(u)) flags[u] = 1;
                    if (!flags(v)) flags[v] = 1;
                  });
  pbbs::scan_add(flags, flags);

  size_t num_ns_clusters = flags[num_clusters];  // num non-singleton clusters
  auto mapping = sequence<uintE>(num_ns_clusters);
  parallel_for_bc(i, 0, num_clusters,
                  (num_clusters > pbbs::kSequentialForThreshold), {
                    if (flags[i] != flags[i + 1]) {
                      mapping[flags[i]] = i;
                    }
                  });

  auto sym_edges = sequence<K>(2 * edges.size(), [&](size_t i) {
    size_t src_edge = i / 2;
    if (i % 2) {
      return std::make_tuple(flags[std::get<0>(std::get<0>(edges[src_edge]))],
                             flags[std::get<1>(std::get<0>(edges[src_edge]))]);
    } else {
      return std::make_tuple(flags[std::get<1>(std::get<0>(edges[src_edge]))],
                             flags[std::get<0>(std::get<0>(edges[src_edge]))]);
    }
  });

  auto EA = edge_array<pbbs::empty>(
      (std::tuple<uintE, uintE, pbbs::empty>*)sym_edges.start(),
      num_ns_clusters, num_ns_clusters, sym_edges.size());

  auto GC = sym_graph_from_edges<pbbs::empty>(EA);
  return std::make_tuple(GC, std::move(flags), std::move(mapping));
}

template <template <class W> class vertex, class W>
inline sequence<uintE> CC_impl(graph<vertex<W>>& GA, double beta,
                                 size_t level, bool pack = false,
                                 bool permute = false) {
  size_t n = GA.n;
  permute |= (level > 0);
  timer ldd_t;
  ldd_t.start();
  auto clusters = LDD(GA, beta, permute, pack);
  ldd_t.stop();
  ldd_t.reportTotal("ldd time");

  timer relabel_t;
  relabel_t.start();
  size_t num_clusters = RelabelIds(clusters);
  relabel_t.stop();
  relabel_t.reportTotal("relabel time");

  timer contract_t;
  contract_t.start();
  auto oracle = [&](const uintE& u, const uintE& v, const W& wgh) {
    return true;  // noop oracle
  };
  auto c_out = contract(GA, clusters, num_clusters, oracle);
  contract_t.stop();
  contract_t.reportTotal("contract time");
  // flags maps from clusters -> no-singleton-clusters
  auto GC = std::get<0>(c_out);
  auto flags = std::get<1>(c_out);
  auto mapping = std::get<2>(c_out);

  if (GC.m == 0) return clusters;

  auto new_labels = CC_impl(GC, beta, level + 1);
  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i) {
    uintE cluster = clusters[i];
    uintE gc_cluster = flags[cluster];
    if (gc_cluster != flags[cluster + 1]) {  // was not a singleton
      // new_labels[gc_cluster] is the gc vertex that captured the whole
      // component. mapping maps this back to the original label range.
      clusters[i] = mapping[new_labels[gc_cluster]];
    }
  });
  GC.del();
  new_labels.del();
  flags.del();
  mapping.del();
  return clusters;
}

template <class Seq>
inline size_t num_cc(Seq& labels) {
  size_t n = labels.size();
  auto flags = sequence<uintE>(n + 1, [&](size_t i) { return 0; });
  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i) {
    if (!flags[labels[i]]) {
      flags[labels[i]] = 1;
    }
  });
  pbbs::scan_add(flags, flags);
  std::cout << "n_cc = " << flags[n] << "\n";
}

template <class Seq>
inline size_t largest_cc(Seq& labels) {
  size_t n = labels.size();
  // could histogram to do this in parallel.
  auto flags = sequence<uintE>(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    flags[labels[i]] += 1;
  }
  std::cout << "largest_cc has size: " << pbbs::reduce_max(flags) << "\n";
}

template <class vertex>
inline sequence<uintE> CC(graph<vertex>& GA, double beta = 0.2,
                            bool pack = false) {
  return CC_impl(GA, beta, 0, pack, true);
}

template <class vertex, class EO>
inline sequence<uintE> CC_oracle_impl(graph<vertex>& GA, EO& oracle,
                                        double beta, size_t level,
                                        bool pack = false,
                                        bool permute = false) {
  size_t n = GA.n;
  if (level > 0) {
    permute = true;
  }
  auto clusters = LDD_impl(GA, oracle, beta, permute, pack);

  timer relabel_t;
  relabel_t.start();
  size_t num_clusters = RelabelIds(clusters);
  relabel_t.stop();
  relabel_t.reportTotal("relabel time");

  timer contract_t;
  contract_t.start();
  auto c_out = contract(GA, clusters, num_clusters, oracle);
  contract_t.stop();
  contract_t.reportTotal("contract time");
  // flags maps from clusters -> no-singleton-clusters
  auto GC = std::get<0>(c_out);
  auto flags = std::get<1>(c_out);
  auto mapping = std::get<2>(c_out);

  if (GC.m == 0) return clusters;

  auto new_labels = cc::CC_impl(GC, beta, level + 1);
  par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i) {
    uintE cluster = clusters[i];
    uintE gc_cluster = flags[cluster];
    if (gc_cluster != flags[cluster + 1]) {  // was not a singleton
      // new_labels[gc_cluster] is the gc vertex that captured the whole
      // component. mapping maps this back to the original label range.
      clusters[i] = mapping[new_labels[gc_cluster]];
    }
  });
  GC.del();
  new_labels.del();
  flags.del();
  mapping.del();
  return clusters;
}

// Connectivity with an edge oracle
template <class vertex, class EO>
inline sequence<uintE> CC_oracle(graph<vertex>& GA, EO& oracle,
                                   double beta = 0.2, bool pack = false) {
  return CC_oracle_impl(GA, oracle, beta, 0, pack, true);
}

}  // namespace cc
