#pragma once

#include "pbbslib/sparse_table.h"

namespace contract {

// Relabel ids, a length n array with values in the range [0,n), to be in the
// range [0..u-1] where u is the number of unique ids
// Returns u
template <class Seq>
inline size_t RelabelIds(Seq& ids) {
  using T = typename Seq::value_type;
  size_t n = ids.size();
  auto inverse_map = sequence<T>(n + 1);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { inverse_map[i] = 0; });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    if (!inverse_map[ids[i]]) inverse_map[ids[i]] = 1;
  });
  pbbslib::scan_add_inplace(inverse_map);

  size_t new_n = inverse_map[n];
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { ids[i] = inverse_map[ids[i]]; });
  return new_n;
}

// Contract when the numbers of clusters is < 2048
template <template <typename W> class vertex, class W>
inline std::tuple<graph<symmetricVertex<pbbslib::empty>>, sequence<uintE>,
                  sequence<uintE>>
contract_small(graph<vertex<W>>& GA, sequence<uintE>& clusters, size_t num_clusters) {
  debug(cout << "Running contract small" << endl;);
  using K = uint32_t; // combine both endpoints into a single key
  using V = pbbslib::empty;

  size_t n = GA.n;

  // Worst case is if the cluster-graph is a clique
  size_t m_upper_bound = 2048*2048;

  auto hash_key = [](const uintE& t) {
    return pbbslib::hash32_2(t);
  };
  std::tuple<uintE, pbbs::empty> empty = std::make_tuple(UINT_E_MAX, pbbs::empty());
  auto edge_table = make_sparse_table<K, V>(m_upper_bound, empty, hash_key);

  timer ins_t; ins_t.start();
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
      uintE c_src = clusters[src];
      uintE c_ngh = clusters[ngh];
      if (c_src < c_ngh) {
        uintE edge_key = (c_ngh << 10) + c_src;
        edge_table.insert(std::make_tuple(edge_key, pbbslib::empty()));
      }
  };
  par_for(0, n, 1, [&] (size_t i) { GA.V[i].mapOutNgh(i, map_f); });
  auto edges = edge_table.entries();
  edge_table.del();
  ins_t.stop(); debug(ins_t.reportTotal("insertion time"););
  debug(cout << "edges.size = " << edges.size() << endl);

  uintE mask = (1 << 10) - 1;
  // Pack out singleton clusters
  auto flags = sequence<uintE>(num_clusters + 1, [](size_t i) { return 0; });

  par_for(0, edges.size(), pbbslib::kSequentialForThreshold, [&] (size_t i) {
                    uintE e = std::get<0>(edges[i]);
                    uintE u = e & mask;
                    uintE v = e >> 10;
                    if (!flags[u]) flags[u] = 1;
                    if (!flags[v]) flags[v] = 1;
                  });
  pbbslib::scan_add_inplace(flags);

  size_t num_ns_clusters = flags[num_clusters];  // num non-singleton clusters
  auto mapping = sequence<uintE>(num_ns_clusters);
  par_for(0, num_clusters, pbbslib::kSequentialForThreshold, [&] (size_t i) {
                    if (flags[i] != flags[i + 1]) {
                      mapping[flags[i]] = i;
                    }
                  });

  auto sym_edges = sequence<std::tuple<uintE, uintE>>(2 * edges.size(), [&](size_t i) {
    size_t src_edge = i / 2;
    uintE e = std::get<0>(edges[src_edge]);
    uintE u = e & mask;
    uintE v = e >> 10;
    if (i % 2) {
      return std::make_tuple(flags[u], flags[v]);
    } else {
      return std::make_tuple(flags[v], flags[u]);
    }
  });

  auto EA = edge_array<pbbslib::empty>(
      (std::tuple<uintE, uintE, pbbslib::empty>*)sym_edges.begin(),
      num_ns_clusters, num_ns_clusters, sym_edges.size());

  auto GC = sym_graph_from_edges<pbbslib::empty>(EA);
  return std::make_tuple(GC, std::move(flags), std::move(mapping));
}

template <template <typename W> class vertex, class W>
inline std::tuple<graph<symmetricVertex<pbbslib::empty>>, sequence<uintE>,
                  sequence<uintE>>
contract_te(graph<vertex<W>>& GA, sequence<uintE>& clusters, size_t num_clusters) {
  // Remove duplicates by hashing
  using K = std::tuple<uintE, uintE>;
  using V = pbbslib::empty;
  using KV = std::tuple<K, V>;

  size_t n = GA.n;

  if (num_clusters < 1024) {
    return contract_small(GA, clusters, num_clusters);
  }

  debug(cout << "num_clusters = " << num_clusters << endl;);
  timer count_t;
  count_t.start();
  auto deg_map = sequence<uintE>(n + 1);
  auto pred = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    return c_src < c_ngh;
  };
  par_for(0, n, 1, [&] (size_t i)
                  { deg_map[i] = GA.V[i].countOutNgh(i, pred); });
  deg_map[n] = 0;
  pbbslib::scan_add_inplace(deg_map);
  count_t.stop();
  debug(count_t.reportTotal("count time"););

  timer ins_t;
  ins_t.start();
  KV empty =
      std::make_tuple(std::make_tuple(UINT_E_MAX, UINT_E_MAX), pbbslib::empty());
  auto hash_pair = [](const std::tuple<uintE, uintE>& t) {
    size_t l = std::min(std::get<0>(t), std::get<1>(t));
    size_t r = std::max(std::get<0>(t), std::get<1>(t));
    size_t key = (l << 32) + r;
    return pbbslib::hash64_2(key);
//    return pbbslib::hash64(std::get<0>(t)) ^ pbbslib::hash64(std::get<1>(t));
  };
  auto edge_table = make_sparse_table<K, V>(deg_map[n], empty, hash_pair);
  debug(cout << "sizeof table = " << edge_table.m << endl;);
  deg_map.clear();

  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    if (c_src < c_ngh) {
      edge_table.insert(
          std::make_tuple(std::make_tuple(c_src, c_ngh), pbbslib::empty()));
    }
  };
  par_for(0, n, 1, [&] (size_t i) { GA.V[i].mapOutNgh(i, map_f); });
  auto edges = edge_table.entries();
  edge_table.del();
  ins_t.stop();
  debug(ins_t.reportTotal("ins time"););
  debug(cout << "edges.size = " << edges.size() << endl);

  // Pack out singleton clusters
  auto flags = sequence<uintE>(num_clusters + 1, [](size_t i) { return 0; });

  par_for(0, edges.size(), pbbslib::kSequentialForThreshold, [&] (size_t i) {
                    auto e = std::get<0>(edges[i]);
                    uintE u = std::get<0>(e);
                    uintE v = std::get<1>(e);
                    if (!flags[u]) flags[u] = 1;
                    if (!flags[v]) flags[v] = 1;
                  });
  pbbslib::scan_add_inplace(flags.slice());

  size_t num_ns_clusters = flags[num_clusters];  // num non-singleton clusters
  auto mapping = sequence<uintE>(num_ns_clusters);
  par_for(0, num_clusters, pbbslib::kSequentialForThreshold, [&] (size_t i) {
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

  auto EA = edge_array<pbbslib::empty>(
      (std::tuple<uintE, uintE, pbbslib::empty>*)sym_edges.begin(),
      num_ns_clusters, num_ns_clusters, sym_edges.size());

  auto GC = sym_graph_from_edges<pbbslib::empty>(EA);
  return std::make_tuple(GC, std::move(flags), std::move(mapping));
}

template <template <typename W> class vertex, class W>
inline std::tuple<graph<symmetricVertex<pbbslib::empty>>, sequence<uintE>,
                  sequence<uintE>>
contract(graph<vertex<W>>& GA, sequence<uintE>& clusters, size_t num_clusters) {
  // Remove duplicates by hashing
  using K = std::tuple<uintE, uintE>;
  using V = pbbslib::empty;
  using KV = std::tuple<K, V>;

  size_t n = GA.n;

  if (num_clusters < 1024) {
    return contract_small(GA, clusters, num_clusters);
  }

  debug(cout << "num_clusters = " << num_clusters << endl;);

  size_t estimated_edges = num_clusters*5;

  timer ins_t;
  ins_t.start();
  KV empty =
      std::make_tuple(std::make_tuple(UINT_E_MAX, UINT_E_MAX), pbbslib::empty());
  auto hash_pair = [](const std::tuple<uintE, uintE>& t) {
    size_t l = std::min(std::get<0>(t), std::get<1>(t));
    size_t r = std::max(std::get<0>(t), std::get<1>(t));
    size_t key = (l << 32) + r;
    return pbbslib::hash64_2(key);
  };
  auto edge_table = make_sparse_table<K, V>(estimated_edges, empty, hash_pair);
  debug(cout << "sizeof table = " << edge_table.m << endl;);

  bool abort = false;
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    if (c_src < c_ngh) {
      edge_table.insert_check(
          std::make_tuple(std::make_tuple(c_src, c_ngh), pbbslib::empty()), &abort);
    }
  };
  parallel_for(0, n, [&] (size_t i) { GA.V[i].mapOutNgh(i, map_f); }, 1);
  if (abort) {
    debug(cout << "calling contract_te" << endl;);
    return contract_te(GA, clusters, num_clusters);
  }
  auto edges = edge_table.entries();
  edge_table.del();
  ins_t.stop();
  debug(ins_t.reportTotal("ins time"););
  debug(cout << "edges.size = " << edges.size() << endl);

  // Pack out singleton clusters
  auto flags = sequence<uintE>(num_clusters + 1, [](size_t i) { return 0; });

  par_for(0, edges.size(), pbbslib::kSequentialForThreshold, [&] (size_t i) {
                    auto e = std::get<0>(edges[i]);
                    uintE u = std::get<0>(e);
                    uintE v = std::get<1>(e);
                    if (!flags[u]) flags[u] = 1;
                    if (!flags[v]) flags[v] = 1;
                  });
  pbbslib::scan_add_inplace(flags.slice());

  size_t num_ns_clusters = flags[num_clusters];  // num non-singleton clusters
  auto mapping = sequence<uintE>(num_ns_clusters);
  par_for(0, num_clusters, pbbslib::kSequentialForThreshold, [&] (size_t i) {
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

  auto EA = edge_array<pbbslib::empty>(
      (std::tuple<uintE, uintE, pbbslib::empty>*)sym_edges.begin(),
      num_ns_clusters, num_ns_clusters, sym_edges.size());

  auto GC = sym_graph_from_edges<pbbslib::empty>(EA);
  return std::make_tuple(GC, std::move(flags), std::move(mapping));
}

} //namespace contract
