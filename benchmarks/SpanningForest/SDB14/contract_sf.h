#pragma once

#include <tuple>
#include "gbbs/graph.h"
#include "gbbs/helpers/sparse_table.h"

namespace gbbs {
namespace contract_sf {

using edge = std::tuple<uintE, uintE>;

constexpr size_t small_cluster_size = 2048;
constexpr size_t m_upper_bound =
    small_cluster_size * small_cluster_size;  // clique on small clusters

struct hash_pair {
  inline size_t operator()(const std::tuple<uintE, uintE>& t) const {
    size_t l = std::min(std::get<0>(t), std::get<1>(t));
    size_t r = std::max(std::get<0>(t), std::get<1>(t));
    size_t key = (l << 32) + r;
    return parlay::hash64_2(key);
  }
};

// Relabel ids, a length n array with values in the range [0,n), to be in the
// range [0..u-1] where u is the number of unique ids
// Returns u
template <class Seq>
inline size_t RelabelIds(Seq& ids) {
  using T = typename Seq::value_type;
  size_t n = ids.size();
  auto inverse_map = sequence<T>(n + 1);
  parallel_for(0, n, kDefaultGranularity,
               [&](size_t i) { inverse_map[i] = 0; });
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    if (!inverse_map[ids[i]]) inverse_map[ids[i]] = 1;
  });
  parlay::scan_inplace(inverse_map);

  size_t new_n = inverse_map[n];
  parallel_for(0, n, kDefaultGranularity,
               [&](size_t i) { ids[i] = inverse_map[ids[i]]; });
  return new_n;
}

// Fetch edges when the numbers of clusters is < small_cluster_size
template <class Graph, class C, class E>
auto fetch_intercluster_small(Graph& GA, C& clusters, size_t num_clusters,
                              E& edge_mapping) {
  debug(std::cout << "Running fetch edges small" << std::endl;);
  using K = std::pair<uintE, uintE>;
  using V = std::pair<uintE, uintE>;
  using KV = std::tuple<K, V>;
  using W = typename Graph::weight_type;

  assert(num_clusters <= small_cluster_size);
  size_t n = GA.n;

  KV empty = std::make_tuple(std::make_pair(UINT_E_MAX, UINT_E_MAX),
                             std::make_pair(UINT_E_MAX, UINT_E_MAX));

  auto edge_table = gbbs::sparse_table<K, V, hash_pair>(small_cluster_size,
                                                        empty, hash_pair());

  timer ins_t;
  ins_t.start();
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    auto orig_edge = edge_mapping(std::make_pair(src, ngh));
    if (c_src < c_ngh) {
      edge_table.insert(
          std::make_pair(std::make_pair(c_src, c_ngh), orig_edge));
    }
  };
  parallel_for(0, n, 1,
               [&](size_t i) { GA.get_vertex(i).out_neighbors().map(map_f); });

  return edge_table;
}

template <class Graph, class C, class E>
auto fetch_intercluster_te(Graph& GA, C& clusters, size_t num_clusters,
                           E& edge_mapping) {
  debug(std::cout << "Running fetch edges te" << std::endl;);
  using K = std::pair<uintE, uintE>;
  using V = std::pair<uintE, uintE>;
  using KV = std::tuple<K, V>;
  using W = typename Graph::weight_type;

  size_t n = GA.n;

  debug(std::cout << "num_clusters = " << num_clusters << std::endl;);
  timer count_t;
  count_t.start();
  auto deg_map = sequence<uintE>(n + 1);
  auto pred = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    return c_src < c_ngh;
  };
  parallel_for(0, n, 1, [&](size_t i) {
    deg_map[i] = GA.get_vertex(i).out_neighbors().count(pred);
  });
  deg_map[n] = 0;
  parlay::scan_inplace(deg_map);
  count_t.stop();
  debug(count_t.next("count time"););

  timer ins_t;
  ins_t.start();
  KV empty = std::make_tuple(std::make_pair(UINT_E_MAX, UINT_E_MAX),
                             std::make_pair(UINT_E_MAX, UINT_E_MAX));

  auto edge_table =
      gbbs::sparse_table<K, V, hash_pair>(deg_map[n], empty, hash_pair());
  debug(std::cout << "sizeof table = " << edge_table.m << std::endl;);
  deg_map.clear();

  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    auto orig_edge = edge_mapping(std::make_pair(src, ngh));
    if (c_src < c_ngh) {
      edge_table.insert(
          std::make_pair(std::make_pair(c_src, c_ngh), orig_edge));
    }
  };
  parallel_for(0, n, 1,
               [&](size_t i) { GA.get_vertex(i).out_neighbors().map(map_f); });
  return edge_table;
}

template <class Graph, class C, class E>
auto fetch_intercluster(Graph& GA, C& clusters, size_t num_clusters,
                        E& edge_mapping) {
  using K = std::pair<uintE, uintE>;
  using V = std::pair<uintE, uintE>;
  using KV = std::tuple<K, V>;
  using W = typename Graph::weight_type;
  size_t n = GA.n;
  debug(std::cout << "num_clusters = " << num_clusters << std::endl;);
  size_t estimated_edges = num_clusters * 5;

  timer ins_t;
  ins_t.start();
  KV empty = std::make_tuple(std::make_pair(UINT_E_MAX, UINT_E_MAX),
                             std::make_pair(UINT_E_MAX, UINT_E_MAX));

  auto edge_table =
      gbbs::sparse_table<K, V, hash_pair>(estimated_edges, empty, hash_pair());
  debug(std::cout << "sizeof table = " << edge_table.m << std::endl;);

  bool abort = false;
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    auto orig_edge = edge_mapping(std::make_pair(src, ngh));
    if (c_src < c_ngh) {
      edge_table.insert_check(
          std::make_pair(std::make_pair(c_src, c_ngh), orig_edge), &abort);
    }
  };
  parallel_for(0, n, 1,
               [&](size_t i) { GA.get_vertex(i).out_neighbors().map(map_f); });
  if (abort) {
    debug(std::cout << "calling fetch_intercluster_te" << std::endl;);
    return fetch_intercluster_te(GA, clusters, num_clusters, edge_mapping);
  }
  return edge_table;
}

template <class Graph, class E>
inline auto contract(Graph& GA, sequence<uintE>& clusters, size_t num_clusters,
                     E& edge_mapping) {
  // Remove duplicates by hashing
  using K = std::pair<uintE, uintE>;
  using V = std::pair<uintE, uintE>;

  auto table =
      (num_clusters < small_cluster_size)
          ? fetch_intercluster_small(GA, clusters, num_clusters, edge_mapping)
          : fetch_intercluster(GA, clusters, num_clusters, edge_mapping);

  auto edges = table.entries();  // sequence
  size_t edges_size = edges.size();

  // Pack out singleton clusters
  auto flags = sequence<uintE>(num_clusters + 1, static_cast<uintE>(0));

  parallel_for(0, edges_size, kDefaultGranularity, [&](size_t i) {
    auto e = std::get<0>(edges[i]);
    uintE u = e.first;
    uintE v = e.second;
    if (!flags[u]) flags[u] = 1;
    if (!flags[v]) flags[v] = 1;
  });
  parlay::scan_inplace(make_slice(flags));

  size_t num_ns_clusters = flags[num_clusters];  // num non-singleton clusters
  auto mapping = sequence<uintE>(num_ns_clusters);
  parallel_for(0, num_clusters, kDefaultGranularity, [&](size_t i) {
    if (flags[i] != flags[i + 1]) {
      mapping[flags[i]] = i;
    }
  });

  auto sym_edges =
      sequence<std::tuple<uintE, uintE, gbbs::empty>>::from_function(
          2 * edges.size(), [&](size_t i) {
            size_t src_edge = i / 2;
            auto e0 = std::get<0>(edges[src_edge]);
            if (i % 2) {
              return std::make_tuple(flags[e0.first], flags[e0.second],
                                     gbbs::empty());
            } else {
              return std::make_tuple(flags[e0.second], flags[e0.first],
                                     gbbs::empty());
            }
          });

  auto GC = sym_graph_from_edges<gbbs::empty>(sym_edges, num_ns_clusters);

  debug(std::cout << "table.size = " << table.m << std::endl;);
  auto ret_table =
      gbbs::sparse_table<K, V, hash_pair>(table.m, table.empty, hash_pair());
  // Go through the edge table and map edges to their new ids
  parallel_for(0, table.m, [&](size_t i) {
    auto& e = table.table[i];
    if (e != table.empty) {
      auto& e0 = std::get<0>(e);
      auto& e1 = std::get<1>(e);
      uintE u = flags[e0.first];
      uintE v = flags[e0.second];
      uintE fst = std::min(u, v);
      uintE snd = std::max(u, v);
      ret_table.insert(std::make_tuple(std::make_pair(fst, snd), e1));
    }
  });

  table.clear_table();

  return std::make_pair(std::move(GC), std::move(ret_table));
}

}  // namespace contract
}  // namespace gbbs
