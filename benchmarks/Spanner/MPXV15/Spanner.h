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

#include "benchmarks/LowDiameterDecomposition/MPX13/LowDiameterDecomposition.h"
#include "gbbs/gbbs.h"
#include "gbbs/pbbslib/dyn_arr.h"
#include "gbbs/pbbslib/sparse_table.h"

namespace gbbs {
namespace spanner {

using edge = std::pair<uintE, uintE>;

// cluster, parent
struct cluster_and_parent {
  uintE cluster;
  uintE parent;
  cluster_and_parent(uintE _cluster, uintE _parent) : cluster(_cluster), parent(_parent) { }
};

template <class Graph, class C>
pbbs::sequence<edge> fetch_intercluster_te(Graph& G, C& clusters, size_t num_clusters) {
  using W = typename Graph::weight_type;
  debug(cout << "Running fetch edges te" << std::endl;);
  using K = edge;
  using V = edge;
  using KV = std::tuple<K, V>;

  size_t n = G.n;

  debug(cout << "num_clusters = " << num_clusters << std::endl;);
  timer count_t;
  count_t.start();
  auto deg_map = pbbs::sequence<uintE>(n + 1);
  auto pred = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    return c_src < c_ngh;
  };
  par_for(0, n, 1, [&] (size_t i)
                  { deg_map[i] = G.get_vertex(i).countOutNgh(i, pred); });
  deg_map[n] = 0;
  pbbslib::scan_add_inplace(deg_map);
  count_t.stop();
  debug(count_t.reportTotal("count time"););

  timer ins_t;
  ins_t.start();
  auto empty_edge = std::make_pair(UINT_E_MAX, UINT_E_MAX);
  KV empty = std::make_tuple(empty_edge, empty_edge);
  auto hash_pair = [](const edge& t) {
    size_t l = std::min(t.first, t.second);
    size_t r = std::max(t.first, t.second);
    size_t key = (l << 32) + r;
    return pbbslib::hash64_2(key);
  };
  auto edge_table = pbbslib::make_sparse_table<K, V>(deg_map[n], empty, hash_pair);
  debug(cout << "sizeof table = " << edge_table.m << std::endl;);
  deg_map.clear();

  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    if (c_src < c_ngh) {
      edge_table.insert(
          std::make_tuple(std::make_pair(c_src, c_ngh), std::make_pair(src, ngh)));
    }
  };
  par_for(0, n, 1, [&] (size_t i) { G.get_vertex(i).mapOutNgh(i, map_f); });
  auto edge_pairs = edge_table.entries();
  edge_table.del();
  ins_t.stop();
  debug(ins_t.reportTotal("ins time"););
  debug(cout << "edges.size = " << edge_pairs.size() << endl);

  auto edges = pbbs::sequence<edge>(edge_pairs.size(), [&] (size_t i) {
    return std::get<1>(edge_pairs[i]);
  });
  return edges;
}

template <class Graph, class C>
pbbs::sequence<edge> fetch_intercluster(Graph& G, C& clusters, size_t num_clusters) {
  using K = edge;
  using V = edge;
  using KV = std::tuple<K, V>;
  using W = typename Graph::weight_type;

  size_t n = G.n;
  debug(cout << "num_clusters = " << num_clusters << std::endl;);
  size_t estimated_edges = num_clusters*5;

  timer ins_t;
  ins_t.start();
  auto empty_edge = std::make_pair(UINT_E_MAX, UINT_E_MAX);
  KV empty = std::make_tuple(empty_edge, empty_edge);
  auto hash_pair = [](const edge& t) {
    size_t l = std::min(t.first, t.second);
    size_t r = std::max(t.first, t.second);
    size_t key = (l << 32) + r;
    return pbbslib::hash64_2(key);
  };

  auto edge_table = pbbslib::make_sparse_table<K, V>(estimated_edges, empty, hash_pair);
  debug(cout << "sizeof table = " << edge_table.m << std::endl;);

  bool abort = false;
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
    uintE c_src = clusters[src];
    uintE c_ngh = clusters[ngh];
    if (c_src < c_ngh) {
      edge_table.insert_check(
          std::make_tuple(std::make_pair(c_src, c_ngh), std::make_pair(src, ngh)), &abort);
    }
  };
  parallel_for(0, n, [&] (size_t i) { G.get_vertex(i).mapOutNgh(i, map_f); }, 1);
  if (abort) {
    debug(cout << "calling fetch_intercluster_te" << std::endl;);
    return fetch_intercluster_te(G, clusters, num_clusters);
  }
  auto edge_pairs = edge_table.entries();
  edge_table.del();
  ins_t.stop();
  debug(ins_t.reportTotal("ins time"););
  debug(cout << "edges.size = " << edge_pairs.size() << endl);

  auto edges = pbbs::sequence<edge>(edge_pairs.size(), [&] (size_t i) {
    return std::get<1>(edge_pairs[i]);
  });
  return edges;
}

template <class Graph>
pbbs::sequence<edge> tree_and_intercluster_edges(Graph& G,
    pbbs::sequence<cluster_and_parent>& cluster_and_parents) {
  size_t n = G.n;
  auto edge_list = pbbslib::dyn_arr<edge>(2*n);

  // Compute and add in tree edges.
  auto tree_edges_with_loops = pbbs::delayed_seq<edge>(n, [&] (size_t i) {
      return std::make_pair(i, cluster_and_parents[i].parent); });
  auto tree_edges = pbbs::filter(tree_edges_with_loops, [&] (const edge& e) {
    return e.first != e.second;
  });
  edge_list.copyIn(tree_edges, tree_edges.size());

  // Compute inter-cluster using hashing.
  auto clusters = pbbs::delayed_seq<uintE>(n, [&] (size_t i) {
    return cluster_and_parents[i].cluster;
  });
  pbbs::sequence<bool> flags(n, false);
  parallel_for(0, n, [&] (size_t i) {
    uintE cluster = clusters[i];
    if (!flags[cluster]) {
      flags[cluster] = true;
    }
  });
  auto cluster_size_seq = pbbs::delayed_seq<size_t>(n, [&] (size_t i) {
    return static_cast<size_t>(flags[i]);
  });
  size_t num_clusters = pbbs::reduce(cluster_size_seq, pbbs::addm<size_t>());

  auto intercluster = fetch_intercluster(G, clusters, num_clusters);
  debug(cout << "num_intercluster edges = " << intercluster.size() << std::endl;);
  edge_list.copyIn(intercluster, intercluster.size());
  size_t edge_list_size = edge_list.size;
  return pbbs::sequence<edge>(edge_list.A, edge_list_size);
}

template <class W>
struct LDD_Parents_F {
  cluster_and_parent* clusters;

  LDD_Parents_F(cluster_and_parent* _clusters) : clusters(_clusters) {}

  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    clusters[d].cluster = clusters[s].cluster;
    clusters[d].parent = s;
    return true;
  }

  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    if (pbbslib::atomic_compare_and_swap(&clusters[d].cluster, UINT_E_MAX, clusters[s].cluster)) {
      clusters[d].parent = s;
      return true;
    }
    return false;
  }

  inline bool cond(uintE d) { return clusters[d].cluster == UINT_E_MAX; }
};

template <class Graph>
inline pbbs::sequence<cluster_and_parent> LDD_parents(Graph& G, double beta, bool permute = true) {
  using W = typename Graph::weight_type;
  size_t n = G.n;

  pbbs::sequence<uintE> vertex_perm;
  if (permute) {
    vertex_perm = pbbslib::random_permutation<uintE>(n);
  }
  auto shifts = ldd_utils::generate_shifts(n, beta);
  auto clusters = pbbs::sequence<cluster_and_parent>(n, cluster_and_parent(UINT_E_MAX, UINT_E_MAX));

  size_t round = 0, num_visited = 0;
  vertexSubset frontier(n);  // Initially empty
  size_t num_added = 0;
  while (num_visited < n) {
    size_t start = shifts[round];
    size_t end = std::min(static_cast<size_t>(shifts[round + 1]), n);
    size_t num_to_add = end - start;
    if (num_to_add > 0) {
      assert((num_added + num_to_add) <= n);
      auto candidates_f = [&](size_t i) {
        if (permute)
          return vertex_perm[num_added + i];
        else
          return static_cast<uintE>(num_added + i);
      };
      auto candidates = pbbslib::make_sequence<uintE>(num_to_add, candidates_f);
      auto pred = [&](uintE v) { return clusters[v].cluster == UINT_E_MAX; };
      auto new_centers = pbbslib::filter(candidates, pred);
      add_to_vsubset(frontier, new_centers.begin(), new_centers.size());
      par_for(0, new_centers.size(), pbbslib::kSequentialForThreshold,
        [&] (size_t i) {
            uintE v = new_centers[i];
            clusters[new_centers[i]] = cluster_and_parent(v, v);
        });
      num_added += num_to_add;
    }

    num_visited += frontier.size();
    if (num_visited >= n) break;

    auto ldd_f = LDD_Parents_F<W>(clusters.begin());
    vertexSubset next_frontier =
        edgeMap(G, frontier, ldd_f, -1, sparse_blocked);
    frontier.del();
    frontier = next_frontier;

    round++;
  }
  return clusters;
}

template <class Graph>
inline pbbs::sequence<edge> Spanner_impl(Graph& G, double beta) {
  bool permute = true;
  timer ldd_t;
  ldd_t.start();
  auto clusters_and_parents = LDD_parents(G, beta, permute);
  ldd_t.stop();
  debug(ldd_t.reportTotal("ldd time"););

  timer build_el_t;
  build_el_t.start();
  auto spanner_edges = tree_and_intercluster_edges(G, clusters_and_parents);
  build_el_t.stop();
  debug(build_el_t.reportTotal("build spanner edges time"););

  // return spanner as an edge-list.
  debug(cout << "Spanner size = " << spanner_edges.size() << std::endl;);
  return spanner_edges;
}

template <class Graph>
inline pbbs::sequence<edge> Spanner(Graph& G, double beta) {
  return Spanner_impl(G, beta);
}

}  // namespace cc
}  // namespace gbbs
