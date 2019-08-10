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

#include "../benchmark/LDD.h"
#include "GBBSCC.h"
#include "pbbslib/sparse_table.h"
#include "pbbslib/dyn_arr.h"
#include "ligra.h"
#include "utils/contract_sf.h"

namespace gbbs_sf {

  using edge = std::pair<uintE, uintE>;

  template <class W>
  struct LDD_Edge_F {
    uintE* cluster_ids;
    uintE* parents;

    LDD_Edge_F(uintE* _cluster_ids, uintE* _parents)
        : cluster_ids(_cluster_ids), parents(_parents) {}

    inline bool update(const uintE& s, const uintE& d, const W& wgh) {
      cluster_ids[d] = cluster_ids[s];
      parents[d] = s;
      return true;
    }

    inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
      if (pbbslib::atomic_compare_and_swap(&cluster_ids[d], UINT_E_MAX, cluster_ids[s])) {
        parents[d] = s;
        return true;
      }
      return false;
    }

    inline bool cond(uintE d) { return cluster_ids[d] == UINT_E_MAX; }
  };

  // Returns a pair containing the cluster_ids and parents.
  template <class G>
  inline std::pair<sequence<uintE>, sequence<uintE>> LDD_edges(G& GA,
      double beta, bool permute = true, bool pack = false) {
    using W = typename G::weight_type;
    size_t n = GA.n;

    sequence<uintE> vertex_perm;
    if (permute) {
      vertex_perm = pbbslib::random_permutation<uintE>(n);
    }
    auto shifts = ldd_utils::generate_shifts(n, beta);
    auto cluster_ids = sequence<uintE>(n);
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { cluster_ids[i] = UINT_E_MAX; });

    auto parents = sequence<uintE>(n, UINT_E_MAX);

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
        auto pred = [&](uintE v) { return cluster_ids[v] == UINT_E_MAX; };
        auto new_centers = pbbslib::filter(candidates, pred);
        add_to_vsubset(frontier, new_centers.begin(), new_centers.size());
        par_for(0, new_centers.size(), pbbslib::kSequentialForThreshold, [&] (size_t i) {
          uintE v = new_centers[i];
          cluster_ids[v] = v;
          parents[v] = v;
        });
        num_added += num_to_add;
      }

      num_visited += frontier.size();
      if (num_visited >= n) break;

      auto ldd_f = LDD_Edge_F<W>(cluster_ids.begin(), parents.begin());
      vertexSubset next_frontier =
          edgeMap(GA, frontier, ldd_f, -1, sparse_blocked);
      if (pack) {
        auto pred = [&](const uintE& src, const uintE& dest, const W& w) {
          return (cluster_ids[src] != cluster_ids[dest]);
        };
       timer t; t.start();
        edgeMapFilter(GA, frontier, pred, pack_edges | no_output);
        t.stop(); debug(t.reportTotal("pack time"););
      }
      frontier.del();
      frontier = next_frontier;

      round++;
    }
    return std::make_pair(cluster_ids, parents);
  }

//  template <template <typename W> class vertex, class W, class E>
//  inline auto contract(graph<vertex<W>>& GA, sequence<uintE>& clusters, size_t num_clusters, E& edge_mapping) {
//    // Remove duplicates by hashing
//    using K = std::pair<uintE, uintE>;
//    using V = std::pair<uintE, uintE>;
//    using KV = std::tuple<K, V>;
//
//    size_t n = GA.n;
//
//    size_t et_size = 0;
//    constexpr size_t small_nclusters = 2048;
//
//    if (num_clusters < small_nclusters) {
//      et_size = small_nclusters^2;
//    } else {
//      debug(cout << "num_clusters = " << num_clusters << endl;);
//      timer count_t;
//      count_t.start();
//      auto deg_map = sequence<uintE>(n + 1);
//      auto pred = [&](const uintE& src, const uintE& ngh, const W& w) {
//        uintE c_src = clusters[src];
//        uintE c_ngh = clusters[ngh];
//        return c_src < c_ngh;
//      };
//      par_for(0, n, 1, [&] (size_t i)
//                      { deg_map[i] = GA.V[i].countOutNgh(i, pred); });
//      deg_map[n] = 0;
//      et_size = pbbslib::reduce_add(deg_map.slice());
//      count_t.stop();
//      debug(count_t.reportTotal("count time"););
//    }
//
//    timer ins_t;
//    ins_t.start();
//    KV empty =
//        std::make_tuple(std::make_pair(UINT_E_MAX, UINT_E_MAX), std::make_pair(UINT_E_MAX, UINT_E_MAX));
//    auto hash_pair = [](const edge& t) {
//      size_t l = std::min(t.first,t.second);
//      size_t r = std::max(t.first,t.second);
//      size_t key = (l << 32) + r;
//      return pbbslib::hash64_2(key);
//    };
//    auto edge_table = make_sparse_table<K, V>(et_size, empty, hash_pair);
//    debug(cout << "sizeof table = " << edge_table.m << endl;);
//
//    auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
//      uintE c_src = clusters[src];
//      uintE c_ngh = clusters[ngh];
//      auto orig_edge = edge_mapping(std::make_pair(src,ngh));
//      if (c_src < c_ngh) {
//        edge_table.insert(
//            std::make_pair(std::make_pair(c_src, c_ngh), orig_edge));
//      }
//    };
//    par_for(0, n, 512, [&] (size_t i) { GA.V[i].mapOutNgh(i, map_f); });
//    auto edges = edge_table.entries();
//    ins_t.stop();
//    debug(ins_t.reportTotal("ins time"););
//
//    // Pack out singleton clusters
//    auto flags = sequence<uintE>(num_clusters + 1, [](size_t i) { return 0; });
//
//    par_for(0, edges.size(), pbbslib::kSequentialForThreshold, [&] (size_t i) {
//                      auto e = std::get<0>(edges[i]);
//                      uintE u = e.first;
//                      uintE v = e.second;
//                      if (!flags[u]) flags[u] = 1;
//                      if (!flags[v]) flags[v] = 1;
//                    });
//    pbbslib::scan_add_inplace(flags.slice());
//
//    size_t num_ns_clusters = flags[num_clusters];  // num non-singleton clusters
//    debug(cout << "num ns_clusters = " << num_ns_clusters << " num orig clusters = " << num_clusters << endl;);
//    debug(cout << "#edges in GC = " << edges.size() << endl;);
//
//    auto sym_edges = sequence<std::tuple<uintE, uintE>>(2 * edges.size(), [&](size_t i) {
//      size_t src_edge = i / 2;
//      auto e0 = std::get<0>(edges[src_edge]);
//      if (i % 2) {
//        return std::make_tuple(flags[e0.first], flags[e0.second]);
//      } else {
//        return std::make_tuple(flags[e0.second], flags[e0.first]);
//      }
//    });
//
//    auto EA = edge_array<pbbslib::empty>(
//        (std::tuple<uintE, uintE, pbbslib::empty>*)sym_edges.begin(),
//        num_ns_clusters, num_ns_clusters, sym_edges.size());
//
//
//    auto GC = sym_graph_from_edges<pbbslib::empty>(EA);
//
//    debug(cout << "etable.size = " << edge_table.m << endl;);
//    auto ret_table = make_sparse_table<K, V>(edge_table.m, empty, hash_pair);
//    // Go through the edge table and map edges to their new ids
//    parallel_for(0, edge_table.m, [&] (size_t i) {
//      auto& e = edge_table.table[i];
//      if (e != edge_table.empty) {
//        auto& e0 = std::get<0>(e);
//        auto& e1 = std::get<1>(e);
//        uintE u = flags[e0.first];
//        uintE v = flags[e0.second];
//        uintE fst = std::min(u,v);
//        uintE snd = std::max(u,v);
//        ret_table.insert(std::make_tuple(std::make_pair(fst, snd), e1));
//      }
//    });
//
//    edge_table.clear();
//
//    return std::make_pair(GC, ret_table);
//  }



  // edge_mapping: edge -> edge
  template <class G>
  inline pbbslib::dyn_arr<edge> SpanningForest_Impl(G& GA, double beta,
                                            size_t level, std::function<edge(edge)>& edge_mapping, bool
                                            pack = false, bool permute = false)
  {
    using W = typename G::weight_type;
    permute |= (level > 0);
    timer ldd_t;
    ldd_t.start();
    auto clusters_and_parents = LDD_edges(GA, beta, permute, pack);
    auto clusters = std::move(clusters_and_parents.first);
    auto parents = std::move(clusters_and_parents.second);
    ldd_t.stop();
    debug(ldd_t.reportTotal("ldd time"););

    // Filter out tree edges added this round (ids are in the current level)
    auto delayed_edges = pbbs::delayed_seq<edge>(parents.size(), [&] (size_t i) {
        return std::make_pair(parents[i], i); });
    auto edges = pbbs::filter(delayed_edges, [&] (const edge& e) { return e.first != e.second; });
    // Apply the mapping to map
    par_for(0, edges.size(), [&] (size_t i) {
      auto e_i = edges[i];
      edges[i] = edge_mapping(e_i);
    });

    timer relabel_t;
    relabel_t.start();
    size_t num_clusters = contract::RelabelIds(clusters);
    relabel_t.stop();
    debug(relabel_t.reportTotal("relabel time"););

    timer contract_t;
    contract_t.start();

    // The contraction here also returns a mapping from edge --> edge. This is
    // because edges incident to a single contracted vertex can come from
    // multiple original vertices.
    auto GC_and_new_mapping = contract_sf::contract(GA, clusters, num_clusters, edge_mapping);
    contract_t.stop();
    debug(contract_t.reportTotal("contract time"););
    auto GC = GC_and_new_mapping.first;
    auto& new_mapping = GC_and_new_mapping.second; // sparse_table<edge, edge>

    size_t edges_size = edges.size();
    if (GC.m == 0) return pbbslib::dyn_arr<edge>(edges.to_array(), edges_size, edges_size, true);

    auto empty_val = std::make_pair(UINT_E_MAX, UINT_E_MAX);
    std::function<edge(edge)> new_edge_mapping = [&] (edge e) {
      uintE u = std::min(e.first, e.second);
      uintE v = std::max(e.first, e.second);
      auto ret = new_mapping.find(std::make_pair(u,v), empty_val);
      assert(ret != empty_val);
      return ret;
    };

    auto rec_edge_arr = SpanningForest_Impl(GC, beta, level + 1, new_edge_mapping);
    rec_edge_arr.copyIn(edges, edges.size());
    GC.del();
    return rec_edge_arr;
  }

  // Algorithm maintains a map from edges in the current graph to original
  // edges (initially just identity).
  template <class G>
  inline pbbslib::dyn_arr<edge> SpanningForest(G& GA, double beta = 0.2,
                                        bool pack = false, bool permute = false) {
    std::function<edge(edge)> identity_mapping = [&] (edge e) {
      return e;
    };
    return SpanningForest_Impl(GA, beta, 0, identity_mapping, pack, permute);
  }

}  // namespace gbbs_sf
