#pragma once

#include "pbbslib/sequence_ops.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/macros.h"
#include "ligra/graph.h"
#include <tuple>
#include "benchmarks/Connectivity/common.h"

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

  // Fetch edges when the numbers of clusters is < small_cluster_size
  template <class Graph, class C>
  std::pair<edge*, size_t> fetch_intercluster_small(Graph& GA, C& clusters, size_t num_clusters) {
    debug(cout << "Running fetch edges small" << endl;);
    using K = std::tuple<uintE, uintE>;
    using V = pbbslib::empty;
    using KV = std::tuple<K, V>;
    using W = typename Graph::weight_type;

    assert(num_clusters <= small_cluster_size);
    size_t n = GA.n;

    KV empty =
        std::make_tuple(std::make_tuple(UINT_E_MAX, UINT_E_MAX), pbbslib::empty());
    auto hash_pair = [](const std::tuple<uintE, uintE>& t) {
      size_t l = std::min(std::get<0>(t), std::get<1>(t));
      size_t r = std::max(std::get<0>(t), std::get<1>(t));
      size_t key = (l << 32) + r;
      return pbbslib::hash64_2(key);
    };
    auto edge_table = make_sparse_table<K, V>(small_cluster_size, empty, hash_pair);

    timer ins_t; ins_t.start();
    auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
      uintE c_src = clusters[src];
      uintE c_ngh = clusters[ngh];
      if (c_src < c_ngh) {
        edge_table.insert(
            std::make_tuple(std::make_tuple(c_src, c_ngh), pbbslib::empty()));
      }
    };
    par_for(0, n, 1, [&] (size_t i) { GA.get_vertex(i).mapOutNgh(i, map_f); });
    auto edges = edge_table.entries();
    edge_table.del();
    ins_t.stop(); debug(ins_t.reportTotal("insertion time"););
    debug(cout << "edges.size = " << edges.size() << endl);

    size_t edge_size = edges.size();
    edge* edge_ret = (edge*)edges.to_array();
    return std::make_pair(edge_ret, edge_size);
  }

  template <class Graph, class C>
  std::pair<edge*, size_t> fetch_intercluster_te(Graph& GA, C& clusters, size_t num_clusters) {
    debug(cout << "Running fetch edges te" << endl;);
    using K = std::tuple<uintE, uintE>;
    using V = pbbslib::empty;
    using KV = std::tuple<K, V>;
    using W = typename Graph::weight_type;

    size_t n = GA.n;

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
                    { deg_map[i] = GA.get_vertex(i).countOutNgh(i, pred); });
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
    par_for(0, n, 1, [&] (size_t i) { GA.get_vertex(i).mapOutNgh(i, map_f); });
    auto edges = edge_table.entries();
    edge_table.del();
    ins_t.stop();
    debug(ins_t.reportTotal("ins time"););
    debug(cout << "edges.size = " << edges.size() << endl);
    size_t edge_size = edges.size();
    edge* edge_ret = (edge*)edges.to_array();
    return std::make_pair(edge_ret, edge_size);
  }

  template <class Graph, class C>
  std::pair<edge*, size_t> fetch_intercluster(Graph& GA, C& clusters, size_t num_clusters) {
    using K = std::tuple<uintE, uintE>;
    using V = pbbslib::empty;
    using KV = std::tuple<K, V>;
    using W = typename Graph::weight_type;
    size_t n = GA.n;
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
    parallel_for(0, n, [&] (size_t i) { GA.get_vertex(i).mapOutNgh(i, map_f); }, 1);
    if (abort) {
      debug(cout << "calling fetch_intercluster_te" << endl;);
      return fetch_intercluster_te(GA, clusters, num_clusters);
    }
    auto edges = edge_table.entries();
    edge_table.del();
    ins_t.stop();
    debug(ins_t.reportTotal("ins time"););
    debug(cout << "edges.size = " << edges.size() << endl);
    size_t edge_size = edges.size();
    edge* edge_ret = (edge*)edges.to_array();
    return std::make_pair(edge_ret, edge_size);
  }

  // Given a graph and a vertex partitioning of the graph, returns a contracted
  // graph where each vertex subset in the partition is contracted into a single
  // vertex. Self-loop edges and duplicate edges are removed, and singleton
  // vertices in the contracted graph also removed.
  //
  // If graph `GA` has `n` vertices, `clusters` should be an `n`-length sequence
  // where `clusters[i]` represents the subset ID of the `i`-th vertex.
  // `num_clusters` should be greater than any value in `clusters`.
  //
  // Returns:
  // (0) The contracted graph.
  // (1) A sequence `S` of length `num_clusters + 1` such that
  //   - if the `i`-th vertex subset does not contract into a singleton vertex,
  //   then it is the `S[i]`-th vertex subset after removing all subsets that
  //   contract into singleton vertices.
  //   - `S[i] == S[i + 1]` iff the `i`-th vertex subset contracts into a
  //   singleton vertex.
  // (2) This sequence `T` is the inverse to (1). If we have a vertex subset
  //   that is the `i`-th subset after removing all subsets that contract to a
  //   singleton, then it is the `T[i]`-th subset before removing those
  //   singleton subsets.
  template <class Graph>
  inline std::tuple<symmetric_graph<symmetric_vertex, pbbslib::empty>, sequence<uintE>, sequence<uintE>>
  contract(Graph& GA, sequence<parent>& clusters, size_t num_clusters) {
    // Remove duplicates by hashing
    using K = std::tuple<uintE, uintE, pbbs::empty>;

    edge* edges;
    size_t edges_size;
    std::tie(edges, edges_size) = (num_clusters < small_cluster_size) ?
      fetch_intercluster_small(GA, clusters, num_clusters) :
      fetch_intercluster(GA, clusters, num_clusters);

    // Pack out singleton clusters
    auto flags = sequence<uintE>(num_clusters + 1, static_cast<uintE>(0));

    par_for(0, edges_size, pbbslib::kSequentialForThreshold, [&] (size_t i) {
                      auto e = edges[i];
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

    auto sym_edges = sequence<K>(2 * edges_size, [&](size_t i) {
      size_t src_edge = i / 2;
      if (i % 2) {
        return std::make_tuple(flags[std::get<0>(edges[src_edge])],
                               flags[std::get<1>(edges[src_edge])],
                               pbbs::empty());
      } else {
        return std::make_tuple(flags[std::get<1>(edges[src_edge])],
                               flags[std::get<0>(edges[src_edge])],
                               pbbs::empty());
      }
    });

    pbbs::free_array(edges);

    auto GC = sym_graph_from_edges<pbbslib::empty>(/* edges = */sym_edges, /* n = */num_ns_clusters);
    return std::make_tuple(GC, std::move(flags), std::move(mapping));
  }

} //namespace contract
