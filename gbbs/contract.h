#pragma once

#include <tuple>

#include "gbbs/graph.h"
#include "gbbs/pbbslib/sparse_table.h"

namespace gbbs {
namespace contract {

  using edge = std::tuple<uintE, uintE>;
  constexpr size_t small_cluster_size = 2048;
  constexpr size_t m_upper_bound = small_cluster_size*small_cluster_size; // clique on small clusters

  // Relabel ids, a length n array with values in the range [0,n), to be in the
  // range [0..u-1] where u is the number of unique ids
  // Returns u
  template <class Seq>
  inline size_t RelabelIds(Seq& ids) {
    using T = typename Seq::value_type;
    size_t n = ids.size();
    auto inverse_map = parlay::sequence<T>::uninitialized(n + 1);
    parallel_for(0, n, [&] (size_t i) { inverse_map[i] = 0; });
    parallel_for(0, n, [&] (size_t i) {
      if (!inverse_map[ids[i]]) inverse_map[ids[i]] = 1;
    });
    parlay::scan_inplace(inverse_map);

    size_t new_n = inverse_map[n];
    parallel_for(0, n, [&] (size_t i)
                    { ids[i] = inverse_map[ids[i]]; });
    return new_n;
  }

  // Fetch edges when the numbers of clusters is < small_cluster_size
  template <class Graph, class C>
  parlay::sequence<std::tuple<std::tuple<uintE, uintE>, gbbs::empty>> fetch_intercluster_small(Graph& GA, C& clusters, size_t num_clusters) {
    debug(std::cout << "# Running fetch edges small" << std::endl;);
    using K = std::tuple<uintE, uintE>;
    using V = gbbs::empty;
    using KV = std::tuple<K, V>;
    using W = typename Graph::weight_type;

    assert(num_clusters <= small_cluster_size);
    size_t n = GA.n;

    KV empty =
        std::make_tuple(std::make_tuple(UINT_E_MAX, UINT_E_MAX), gbbs::empty());
    auto hash_pair = [](const std::tuple<uintE, uintE>& t) {
      size_t l = std::min(std::get<0>(t), std::get<1>(t));
      size_t r = std::max(std::get<0>(t), std::get<1>(t));
      size_t key = (l << 32) + r;
      return parlay::hash64_2(key);
    };
    auto edge_table = pbbslib::make_sparse_table<K, V>(small_cluster_size, empty, hash_pair);

    timer ins_t; ins_t.start();
    auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
      uintE c_src = clusters[src];
      uintE c_ngh = clusters[ngh];
      if (c_src < c_ngh) {
        edge_table.insert(
            std::make_tuple(std::make_tuple(c_src, c_ngh), gbbs::empty()));
      }
    };
    parallel_for(0, n, [&] (size_t i) { GA.get_vertex(i).out_neighbors().map(map_f); }, 1);
    auto edges = edge_table.entries();
    edge_table.del();
    ins_t.stop(); debug(ins_t.reportTotal("insertion time"););
    debug(std::cout << "# edges.size = " << edges.size() << std::endl);

    return edges;
  }

  template <class Graph, class C>
  parlay::sequence<std::tuple<std::tuple<uintE, uintE>, gbbs::empty>> fetch_intercluster_te(Graph& GA, C& clusters, size_t num_clusters) {
    debug(std::cout << "# Running fetch edges te" << std::endl;);
    using K = std::tuple<uintE, uintE>;
    using V = gbbs::empty;
    using KV = std::tuple<K, V>;
    using W = typename Graph::weight_type;

    size_t n = GA.n;

    debug(std::cout << "# num_clusters = " << num_clusters << std::endl;);
    timer count_t;
    count_t.start();
    auto deg_map = sequence<uintE>::uninitialized(n + 1);
    auto pred = [&](const uintE& src, const uintE& ngh, const W& w) {
      uintE c_src = clusters[src];
      uintE c_ngh = clusters[ngh];
      return c_src < c_ngh;
    };
    parallel_for(0, n, [&] (size_t i)
                    { deg_map[i] = GA.get_vertex(i).out_neighbors().count(pred); }, 1);
    deg_map[n] = 0;
    parlay::scan_inplace(deg_map);
    count_t.stop();
    debug(count_t.reportTotal("count time"););

    timer ins_t;
    ins_t.start();
    KV empty =
        std::make_tuple(std::make_tuple(UINT_E_MAX, UINT_E_MAX), gbbs::empty());
    auto hash_pair = [](const std::tuple<uintE, uintE>& t) {
      size_t l = std::min(std::get<0>(t), std::get<1>(t));
      size_t r = std::max(std::get<0>(t), std::get<1>(t));
      size_t key = (l << 32) + r;
      return parlay::hash64_2(key);
    };
    auto edge_table = pbbslib::make_sparse_table<K, V>(deg_map[n], empty, hash_pair);
    debug(std::cout << "# sizeof table = " << edge_table.m << std::endl;);
    deg_map.clear();

    auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
      uintE c_src = clusters[src];
      uintE c_ngh = clusters[ngh];
      if (c_src < c_ngh) {
        edge_table.insert(
            std::make_tuple(std::make_tuple(c_src, c_ngh), gbbs::empty()));
      }
    };
    parallel_for(0, n, [&] (size_t i) { GA.get_vertex(i).out_neighbors().map(map_f); }, 1);
    auto edges = edge_table.entries();
    edge_table.del();
    ins_t.stop();
    debug(ins_t.reportTotal("ins time"););
    debug(std::cout << "# edges.size = " << edges.size() << std::endl);
    return edges;
  }

  template <class Graph, class C>
  parlay::sequence<std::tuple<std::tuple<uintE, uintE>, gbbs::empty>> fetch_intercluster(Graph& GA, C& clusters, size_t num_clusters) {
    using K = std::tuple<uintE, uintE>;
    using V = gbbs::empty;
    using KV = std::tuple<K, V>;
    using W = typename Graph::weight_type;
    size_t n = GA.n;
    debug(std::cout << "# num_clusters = " << num_clusters << std::endl;);
    size_t estimated_edges = num_clusters*5;

    timer ins_t;
    ins_t.start();
    KV empty =
        std::make_tuple(std::make_tuple(UINT_E_MAX, UINT_E_MAX), gbbs::empty());
    auto hash_pair = [](const std::tuple<uintE, uintE>& t) {
      size_t l = std::min(std::get<0>(t), std::get<1>(t));
      size_t r = std::max(std::get<0>(t), std::get<1>(t));
      size_t key = (l << 32) + r;
      return parlay::hash64_2(key);
    };
    auto edge_table = pbbslib::make_sparse_table<K, V>(estimated_edges, empty, hash_pair);
    debug(std::cout << "# sizeof table = " << edge_table.m << std::endl;);

    bool abort = false;
    auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
      uintE c_src = clusters[src];
      uintE c_ngh = clusters[ngh];
      if (c_src < c_ngh) {
        edge_table.insert_check(
            std::make_tuple(std::make_tuple(c_src, c_ngh), gbbs::empty()), &abort);
      }
    };
    parallel_for(0, n, [&] (size_t i) { GA.get_vertex(i).out_neighbors().map(map_f); }, 1);
    if (abort) {
      debug(std::cout << "# calling fetch_intercluster_te" << std::endl;);
      return fetch_intercluster_te(GA, clusters, num_clusters);
    }
    auto edges = edge_table.entries();
    edge_table.del();
    ins_t.stop();
    debug(ins_t.reportTotal("ins time"););
    debug(std::cout << "# edges.size = " << edges.size() << std::endl);
    return edges;
  }

  // Given a graph and a vertex partitioning of the graph, returns a contracted
  // graph where each vertex grouping/cluster in the partition is contracted
  // into a single vertex. Self-loop edges and duplicate edges are removed, and
  // singleton vertices (vertices with degree zero) in the contracted graph also
  // removed.
  //
  // Arguments:
  //   GA
  //     The graph to contract.
  //   clusters
  //     A `GA.n`-length sequence where `clusters[i]` is the cluster ID of the
  //     i-th vertex. Each cluster ID must be in the range [0, `num_clusters`).
  //   num_clusters
  //     The number of clusters.
  //
  // Returns:
  // (0) The contracted graph.
  // (1) A sequence `S` of length `num_clusters + 1` such that
  //   - if cluster i does not contract into a singleton, then it is vertex
  //   `S[i]` in the contracted graph.
  //   - `S[i] == S[i + 1]` iff the cluster `i` contracts into a singleton.
  // (2) A sequence `T` that is the inverse to (1). Vertex i in the contracted
  //   graph is the contraction of cluster `T[i]`.
  template <class Graph>
  inline std::tuple<symmetric_graph<symmetric_vertex, gbbs::empty>, sequence<uintE>, sequence<uintE>>
  contract(Graph& GA, parlay::sequence<uintE>& clusters, size_t num_clusters) {
    // Remove duplicates by hashing
    using K = std::tuple<uintE, uintE, gbbs::empty>;

    auto edges = (num_clusters < small_cluster_size) ?
      fetch_intercluster_small(GA, clusters, num_clusters) :
      fetch_intercluster(GA, clusters, num_clusters);
    auto edges_begin = (std::tuple<uintE, uintE>*)edges.begin();


    // Pack out singleton clusters
    auto flags = sequence<uintE>(num_clusters + 1, static_cast<uintE>(0));

    parallel_for(0, edges.size(), [&] (size_t i) {
                      auto e = edges_begin[i];
                      uintE u = std::get<0>(e);
                      uintE v = std::get<1>(e);
                      if (!flags[u]) flags[u] = 1;
                      if (!flags[v]) flags[v] = 1;
                    });
    parlay::scan_inplace(parlay::make_slice(flags));

    size_t num_ns_clusters = flags[num_clusters];  // num non-singleton clusters
    auto mapping = sequence<uintE>(num_ns_clusters);
    parallel_for(0, num_clusters, [&] (size_t i) {
                      if (flags[i] != flags[i + 1]) {
                        mapping[flags[i]] = i;
                      }
                    });

    auto sym_edges = parlay::sequence<K>::from_function(2 * edges.size(), [&](size_t i) {
      size_t src_edge = i / 2;
      if (i % 2) {
        return std::make_tuple(flags[std::get<0>(edges_begin[src_edge])],
                               flags[std::get<1>(edges_begin[src_edge])],
                               gbbs::empty());
      } else {
        return std::make_tuple(flags[std::get<1>(edges_begin[src_edge])],
                               flags[std::get<0>(edges_begin[src_edge])],
                               gbbs::empty());
      }
    });

//    std::cout << "Creating graph from " << sym_edges.size() << " many edges. Num clusters = " << num_ns_clusters << std::endl;
//    for (size_t i=0; i<sym_edges.size(); i++) {
//      auto [u, v, _] = sym_edges[i];
//      if (u > num_ns_clusters || v > num_ns_clusters) {
//        std::cout << std::get<0>(sym_edges[i]) << " " << std::get<1>(sym_edges[i]) << std::endl;
//
//      }
//    }

    auto GC = sym_graph_from_edges<gbbs::empty>(
        sym_edges, num_ns_clusters);
    return std::make_tuple(GC, std::move(flags), std::move(mapping));
  }

}  // namespace contract
}  // namespace gbbs
