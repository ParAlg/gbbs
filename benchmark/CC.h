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
size_t RelabelIds(Seq& ids) {
  using T = typename Seq::T;
  size_t n = ids.size();
  auto inverse_map = array_imap<T>(n + 1);
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                  { inverse_map[i] = 0; });
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    if (!inverse_map[ids[i]]) inverse_map[ids[i]] = 1;
  });
  pbbs::scan_add(inverse_map, inverse_map);

  size_t new_n = inverse_map[n];
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                  { ids[i] = inverse_map[ids[i]]; });
  return new_n;
}

template <template <typename W> class vertex, class W, class Seq>
auto contract(graph<vertex<W> >& GA, Seq& clusters, size_t num_clusters) {
  // Remove duplicates by hashing
  using K = tuple<uintE, uintE>;
  using V = pbbs::empty;
  using KV = tuple<K, V>;

  size_t n = GA.n;
  size_t m = GA.m;

  timer count_t;
  count_t.start();
  auto deg_map = array_imap<uintE>(n + 1);
  auto pred = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    return c_src < c_ngh;
  };
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                  { deg_map[i] = GA.V[i].countOutNgh(i, pred); });
  deg_map[n] = 0;
  pbbs::scan_add(deg_map, deg_map);
  count_t.stop();
  count_t.reportTotal("count time");

  timer ins_t;
  ins_t.start();
  KV empty = make_tuple(make_tuple(UINT_E_MAX, UINT_E_MAX), pbbs::empty());
  auto hash_pair = [](const tuple<uintE, uintE>& t) {
    return pbbs::hash64(get<0>(t)) ^ pbbs::hash64(get<1>(t));
  };
  auto edge_table = make_sparse_table<K, V>(deg_map[n], empty, hash_pair);

  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    if (c_src < c_ngh) {
      edge_table.insert(make_tuple(make_tuple(c_src, c_ngh), pbbs::empty()));
    }
  };
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                  { GA.V[i].mapOutNgh(i, map_f, true); });
  auto e2 = edge_table.entries();
  auto edges = make_array_imap<K>((K*)e2.start(), e2.size());
  ins_t.stop();
  ins_t.reportTotal("ins time");

  // Pack out singleton clusters
  auto flags = array_imap<uintE>(num_clusters + 1, [](size_t i) { return 0; });
  parallel_for_bc(i, 0, edges.size(),
                  (edges.size() > pbbs::kSequentialForThreshold), {
                    auto e = edges[i];
                    uintE u = get<0>(e) COMMA v = get<1>(e);
                    if (!flags(u)) flags[u] = 1;
                    if (!flags(v)) flags[v] = 1;
                  }) pbbs::scan_add(flags, flags);

  size_t num_ns_clusters = flags[num_clusters];  // num non-singleton clusters
  auto mapping = array_imap<uintE>(num_ns_clusters);
  parallel_for_bc(i, 0, num_clusters,
                  (num_clusters > pbbs::kSequentialForThreshold), {
                    if (flags[i] != flags[i + 1]) {
                      mapping[flags[i]] = i;
                    }
                  });

  auto sym_edges = array_imap<K>(2 * edges.size(), [&](size_t i) {
    size_t src_edge = i / 2;
    if (i % 2) {
      return make_tuple(flags[get<0>(edges[src_edge])],
                        flags[get<1>(edges[src_edge])]);
    } else {
      return make_tuple(flags[get<1>(edges[src_edge])],
                        flags[get<0>(edges[src_edge])]);
    }
  });

  auto EA = edge_array<pbbs::empty>(
      (tuple<uintE, uintE, pbbs::empty>*)sym_edges.start(), num_ns_clusters,
      num_ns_clusters, sym_edges.size());

  auto GC = sym_graph_from_edges<pbbs::empty>(EA);
  return make_tuple(GC, std::move(flags), std::move(mapping));
}

template <class vertex>
auto CC_impl(graph<vertex>& GA, double beta, size_t level, bool pack = false,
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
  auto c_out = contract(GA, clusters, num_clusters);
  contract_t.stop();
  contract_t.reportTotal("contract time");
  // flags maps from clusters -> no-singleton-clusters
  auto GC = get<0>(c_out);
  auto flags = get<1>(c_out);
  auto mapping = get<2>(c_out);

  if (GC.m == 0) return clusters;

  auto new_labels = CC_impl(GC, beta, level + 1);
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    uintE cluster = clusters[i];
    uintE gc_cluster = flags[cluster];
    if (gc_cluster != flags[cluster + 1]) {  // was not a singleton
      // new_labels[gc_cluster] is the gc vertex that captured the whole
      // component. mapping maps this back to the original label range.
      clusters[i] = mapping[new_labels[gc_cluster]];
    }
  });
  return clusters;
}

template <class Seq>
size_t num_cc(Seq& labels) {
  size_t n = labels.size();
  auto flags = array_imap<uintE>(n + 1, [&](size_t i) { return 0; });
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    if (!flags[labels[i]]) {
      flags[labels[i]] = 1;
    }
  });
  pbbs::scan_add(flags, flags);
  size_t n_cc = flags[n];
  cout << "n_cc = " << flags[n] << endl;
}

template <class Seq>
size_t largest_cc(Seq& labels) {
  size_t n = labels.size();
  // could histogram to do this in parallel.
  auto flags = array_imap<uintE>(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    flags[labels[i]] += 1;
  }
  cout << "largest_cc has size: " << pbbs::reduce_max(flags) << endl;
}

template <class vertex>
uintE* CC(graph<vertex>& GA, double beta = 0.2, bool pack = false) {
  long n = GA.n;
  timer t;
  t.start();
  auto labels = CC_impl(GA, beta, 0, pack, true);
  double tt = t.stop();
  return labels.get_array();
}
}  // namespace cc
