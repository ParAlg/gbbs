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

#include <queue>
#include <unordered_set>
#include <vector>

#include "gbbs/gbbs.h"
#include "pam/pam.h"

namespace gbbs {
namespace clustering {

template <class Weights, class IW, template <class W> class w_vertex>
struct clustered_graph {

  using Graph = symmetric_graph<w_vertex, IW>;
  using Graph_vertex = w_vertex<IW>;

  using W = typename Weights::weight_type;
  using edge = std::pair<uintE, W>;

  struct neighbor_entry {
    using key_t = uintE;  // neighbor_id
    using val_t = W;      // weight
    static inline bool comp(key_t a, key_t b) { return a < b; }
  };

  using neighbor_map = pam_map<neighbor_entry>;

  struct clustered_vertex {

    clustered_vertex() {}


    clustered_vertex(uintE vtx_id, Graph_vertex& vertex, const Weights& weights, edge* edges) {
      auto cluster_size = vertex.out_degree();
      auto combine_w = [&] (W l, W r) { return l; };
      if (cluster_size > 0) {
        neighbors = neighbor_map(edges, edges + cluster_size, combine_w);
      } else {
        neighbors = neighbor_map();
      }

      num_in_cluster = 1;
      staleness = 1;
      cas_size = 1;
      current_id = vtx_id;
      active = true;
    }

    struct Add {
      using T = size_t;
      static T identity() { return 0;}
      static T add(T a, T b) { return a + b;}
    };

    // F is a predicate from neighbor -> bool
    template <class F>
    uintE countNeighbors(uintE id, F& f) {
      auto pred = [&] (edge e) -> size_t {
        return f(id, e.first, e.second);
      };
      return neighbor_map::map_reduce(neighbors, pred, Add());
    }

    template <class F>
    void iterate(uintE id, F& f) {
      auto iter = [&] (edge e) -> void {
        f(id, e.first, e.second);
      };
      neighbor_map::foreach_seq(neighbors, iter);
    }

    template <class F>
    void iterate_cond(uintE id, F& f) {
      auto iter = [&] (edge e) -> bool {
        return f(id, e.first, e.second);
      };
      neighbor_map::foreach_cond(neighbors, iter);
    }

    void print_edges() {
      auto f = [&] (const uintE& u, const uintE& v, const auto& wgh) {
        std::cout << u << " " << v << " " << wgh << std::endl;
      };
      iterate(current_id, f);
    }

    template <class F>
    void map_index(uintE id, F& f) {
      auto iter = [&] (edge e, size_t index) -> void {
        f(id, e.first, e.second, index);
      };
      neighbor_map::map_index(neighbors, iter);
    }

    uintE neighbor_size() {
      return neighbors.size();
    }

    // Number of vertices contained in this cluster.
    uintE cluster_size() {
      return num_in_cluster;
    }

    // Tracks the last cluster update size.
    uintE staleness;
    // Initially equal to staleness, but cas'd in parallel rounds.
    uintE cas_size;
    // The "current" id of this cluster, updated upon a merge that keeps this cluster active.
    uintE current_id;
    // Number of vertices contained in this cluster.
    uintE num_in_cluster;
    // Active = false iff this cluster is no longer active.
    bool active;
    // A map storing our neighbors + weights.
    neighbor_map neighbors;
  };

  Graph& G;
  Weights& weights;

  uintE n;
  uintE last_cluster_id;
  uintE num_merges_performed;

  uintE merge_order_idx = 0;
  sequence<clustered_vertex> clusters;
  sequence<std::pair<uintE, float>> merge_order;
  sequence<std::pair<uintE, W>> dendrogram;

  uintE degree(uintE v) {
    return clusters[v].neighbor_size();
  }


  template <class Sim>
  void insert_neighbor_endpoints(sequence<std::tuple<uintE, uintE, Sim>>&& triples) {
    timer nt; nt.start();
    parlay::sort_inplace(make_slice(triples));

    auto all_starts = parlay::delayed_seq<size_t>(triples.size(), [&] (size_t i) {
      if ((i == 0) || std::get<0>(triples[i]) != std::get<0>(triples[i-1])) {
        return i;
      }
      return std::numeric_limits<size_t>::max();
    });
    auto starts = parlay::filter(all_starts, [&] (auto v) { return v != std::numeric_limits<size_t>::max(); });
    std::cout << "Number of neighbor-insertion targets = " << starts.size() << std::endl;

    auto pairs = sequence<std::pair<uintE, Sim>>::from_function(triples.size(), [&] (size_t i) {
      return std::make_pair(std::get<1>(triples[i]), std::get<2>(triples[i]));
    });

    parallel_for(0, starts.size(), [&] (size_t i) {
      size_t start = starts[i];
      size_t end = (i == starts.size() - 1) ? triples.size() : starts[i+1];
      uintE our_id = std::get<0>(triples[start]);

      assert(clusters[our_id].active);

      auto our_map = std::move(clusters[our_id].neighbors);
      auto sl = pairs.cut(start, end);
        auto op = [&] (const uintE& key, const Sim& old_val, const Sim& new_val) {
          return old_val + new_val;  // sum to get the new total weight across the cut
        };
      auto updated = neighbor_map::keyed_multi_insert_sorted(std::move(our_map), sl, op);
      clusters[our_id].neighbors = std::move(updated);
    });
    nt.stop(); nt.reportTotal("Neighbor insertion time");
  }

  // Delete all (u,v) edges from the graph. The "v" endpoint is a newly
  // deactivated vertex, so no need to delete from that side, but may need to
  // delete from the "u" side.
  // - first endpoint: status unknown (either active or newly deactivated)
  // - second endpoint: newly deactivated
  void process_deletions(sequence<std::pair<uintE, uintE>>&& dels) {
    timer dt; dt.start();
    parlay::sort_inplace(make_slice(dels));  // sort by the first endpoint

    auto all_starts = parlay::delayed_seq<size_t>(dels.size(), [&] (size_t i) {
      if ((i == 0) || dels[i].first != dels[i-1].first) {
        return i;
      }
      return std::numeric_limits<size_t>::max();
    });
    auto starts = parlay::filter(all_starts, [&] (auto v) { return v != std::numeric_limits<size_t>::max(); });
    std::cout << "Number of deletion targets = " << starts.size() << std::endl;
    auto del_keys = sequence<uintE>::from_function(dels.size(), [&] (size_t i) { return dels[i].second; });

    parallel_for(0, starts.size(), [&] (size_t i) {
      size_t start = starts[i];
      size_t end = (i == starts.size() - 1) ? dels.size() : starts[i+1];
      uintE our_id = dels[start].first;

      if (clusters[our_id].active) {  // otherwise this vertex is deactivated and we can safely skip.
        // perform deletions.
        auto our_map = std::move(clusters[our_id].neighbors);
        auto sl = del_keys.cut(start, end);
        auto updated = neighbor_map::multi_delete(std::move(our_map), sl);
        clusters[our_id].neighbors = std::move(updated);
      }
    });
    dt.stop(); dt.reportTotal("deletion time");
  }

  // input: sequence of (u, v) pairs representing that v will merge to u
  template <class Sim>
  void unite_merge(sequence<std::tuple<uintE, uintE, float>>&& merge_seq) {
    timer um; um.start();
    std::cout << "Start of unite merge" << std::endl;
    // Sort.
    parlay::sort_inplace(make_slice(merge_seq));

    // For each instance, find largest component.
    auto all_starts = parlay::delayed_seq<uintE>(merge_seq.size(), [&] (size_t i) {
      if ((i == 0) || std::get<0>(merge_seq[i]) != std::get<0>(merge_seq[i-1])) {
        return (uintE)i;
      }
      return UINT_E_MAX;
    });
    auto starts = parlay::filter(all_starts, [&] (uintE v) { return v != UINT_E_MAX; });
    std::cout << "Number of merge targets = " << starts.size() << std::endl;

    auto edge_sizes = sequence<size_t>(merge_seq.size());

    // In parallel over every instance.
    // Make the merge target the cluster with the largest number of out-edges.
    parallel_for(0, starts.size(), [&] (size_t i) {
      size_t start = starts[i];
      size_t end = (i == starts.size() - 1) ? merge_seq.size() : starts[i+1];

      uintE our_id = std::get<0>(merge_seq[start]);
      auto our_size = clusters[our_id].neighbor_size();

      auto sizes_and_id = parlay::delayed_seq<std::pair<uintE, uintE>>(end - start, [&] (size_t i) {
        uintE vtx_id = std::get<1>(merge_seq[start + i]);
        return std::make_pair(clusters[vtx_id].neighbor_size(), vtx_id);
      });
      std::pair<uintE, uintE> id = std::make_pair((uintE)0, (uintE)UINT_E_MAX);
      auto mon = parlay::make_monoid([&] (const auto& l, const auto& r) { return (l.first < r.first) ? r : l;}, id);
      auto [largest_size, largest_id] = parlay::reduce(sizes_and_id, mon);

      if (our_size < largest_size) {  // Relabel merge targets to largest_id.
        for (size_t j=start; j<end; j++) {
          std::get<0>(merge_seq[j]) = largest_id;
          if (std::get<1>(merge_seq[j]) == largest_id) {
            std::get<1>(merge_seq[j]) = our_id;
          }
        }
        our_id = largest_id;
        our_size = largest_size;
      }

      // Write the neighbor size before scan.
      // - Update the current_id for each (being merged) neighbor.
      // - Set the active flag for each such neighbor to false ( Deactivate ).
      size_t total_size = 0;
      for (size_t j=start; j<end; j++) {
        uintE ngh_id = std::get<1>(merge_seq[j]);
        uintE ngh_size = clusters[ngh_id].neighbor_size();
        edge_sizes[j] = ngh_size;

        // Was an active cluster before.
        assert(clusters[ngh_id].current_id == ngh_id);
        assert(clusters[ngh_id].active);

        // Update total_size with the size of the cluster being merge
        total_size += clusters[ngh_id].cluster_size();

        // Update id to point to the merge target, and deactivate.
        clusters[ngh_id].current_id = our_id;
        clusters[ngh_id].active = false;
      }
      // Update our cluster size.
      clusters[our_id].num_in_cluster += total_size;
      // Update the CAS size for the next round.
      clusters[our_id].cas_size = clusters[our_id].num_in_cluster;
    });

//    std::cout << "Edges incident to 2862074" << std::endl;
//    clusters[2862074].print_edges();
//
//    std::cout << "Edges incident to 2619886" << std::endl;
//    clusters[2619886].print_edges();
//
//    std::cout << "Merging neighbors...." << std::endl << std::endl;
//    std::cout << "Edges incident to " << 2862067 << std::endl;
//    clusters[2862067].print_edges();
//
//    std::cout << "Edges incident to " << 2862057 << std::endl;
//    clusters[2862057].print_edges();

//    for (size_t i=0; i<merge_seq.size(); i++) {
//      auto [u, v] = merge_seq[i];
//      if ((u == 2862074) || (u == 2619886)) {
//        std::cout << "Merging: " << u << " and " << v << std::endl;
//      }
//      if ((v == 2862074) || (v == 2619886)) {
//        std::cout << "V merging: " << u << " and " << v << std::endl;
//      }
//    }

    // Scan to compute #edges we need to merge.
    size_t total_edges = parlay::scan_inplace(make_slice(edge_sizes));
    std::cout << "Total edges = " << total_edges << std::endl;
    auto edges = sequence<std::pair<uintE, Sim>>::uninitialized(total_edges);
    // Deletions store deleted edges (incident to newly deactivated vertices) as
    // (v,u) pairs, where u is deactivated, and the status of v is unknown.
    // The last merge_seq.size() many deletions are directly copied from the
    // merges.
    auto deletions = sequence<std::pair<uintE, uintE>>::uninitialized(total_edges + merge_seq.size());

    // Copy edges from trees to edges and deletions.
    parallel_for(0, merge_seq.size(), [&] (size_t i) {
      uintE ngh_id = std::get<1>(merge_seq[i]);
      size_t k = 0;
      size_t off = edge_sizes[i];
      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh, size_t k) {
        uintE cur_ngh_id = clusters[v].current_id;
        // Following symmetry-breaking check is very important to prevent sending a (u,v) edge
        // between two deactivated vertices _twice_ to the activated targets.
        bool ngh_active = clusters[v].active;
        if (ngh_active || (u < v)) {
          edges[off + k] = std::make_pair(cur_ngh_id, wgh);
        } else {
          edges[off + k] = std::make_pair(UINT_E_MAX, 0);
        }
        deletions[off + k] = std::make_pair(v, u);
      };
      clusters[ngh_id].map_index(ngh_id, map_f);
    });

    // Delete every (active_id, deactivated_id) edge incident to the active
    // merge targets.
    parallel_for(0, merge_seq.size(), [&] (size_t i) {
      auto [u, v, wgh] = merge_seq[i];
      deletions[total_edges + i] = std::make_pair(u, v);
    });

    // (1) First perform the deletions. This makes life easier later when we
    // perform the insertions.
    process_deletions(std::move(deletions));

    auto compacted_edge_sizes = sequence<size_t>::uninitialized(starts.size());

    // Sort within each instance, scan to merge weights for identical edges.
    parallel_for(0, starts.size(), [&] (size_t i) {
      size_t sort_start = starts[i];
      size_t sort_end = (i == starts.size() - 1) ? merge_seq.size() : starts[i+1];
      size_t edges_start = edge_sizes[sort_start];
      size_t edges_end = (sort_end == merge_seq.size()) ? total_edges : edge_sizes[sort_end];

      uintE our_id = std::get<0>(merge_seq[sort_start]);
      auto our_size = clusters[our_id].cluster_size();

      auto our_edges = edges.cut(edges_start, edges_end);

      auto comp = [&] (const auto& l, const auto& r) { return l.first < r.first;};
      parlay::sort_inplace(our_edges, comp);

      auto scan_f = [&] (const auto& l, const auto& r) {
        if (r.first == UINT_E_MAX) return l;
        if (r.first != l.first) return r;
        return std::make_pair(r.first, l.second + r.second);
      };
      std::pair<uintE, Sim> id = std::make_pair(UINT_E_MAX, (Sim)0);
      auto scan_mon = parlay::make_monoid(scan_f, id);
      // After the scan, last occurence of each ngh has the sum'd weight.
      parlay::scan_inclusive_inplace(our_edges, scan_mon);

      // Pack-inplace. (optimize later if it seems necessary)
      // Also get rid of self-loop edges.
      assert(our_edges.size() > 0);
      size_t k = 0;
      for (size_t j=0; j<our_edges.size(); j++) {
        if (our_edges[j].first == our_id || our_edges[j].first == UINT_E_MAX) {
          continue;  // skip self-loops and deleted edges
        }
        if ((j == our_edges.size() - 1) || (our_edges[j].first != our_edges[j+1].first)) {
          our_edges[k] = our_edges[j];
          k++;
        }
      }
      compacted_edge_sizes[i] = k;

      if (k > 0) {
        auto our_map = std::move(clusters[our_id].neighbors);
        auto sl = edges.cut(edges_start, edges_start + k);
        auto op = [&] (const uintE& key, const Sim& old_val, const Sim& new_val) {
          return old_val + new_val;  // sum to get the new total weight across the cut
        };
        auto updated = neighbor_map::keyed_multi_insert_sorted(std::move(our_map), sl, op);
        clusters[our_id].neighbors = std::move(updated);
      }
    });

    size_t total_compacted = parlay::scan_inplace(make_slice(compacted_edge_sizes));

    // Map to triples, and reinsert correctly weighted edges into neighbors
    // endpoints (transposed).
    using triple = std::tuple<uintE, uintE, Sim>;
    triple id{UINT_E_MAX, UINT_E_MAX, 0};
    auto triples = sequence<triple>(total_compacted, id);
    std::cout << "Total compacted (triples size) = " << total_compacted << std::endl;

    parallel_for(0, starts.size(), [&] (size_t i) {
      size_t sort_start = starts[i];
      size_t sort_end = (i == starts.size() - 1) ? merge_seq.size() : starts[i+1];
      size_t edges_start = edge_sizes[sort_start];

      size_t num_edges = ((i == starts.size() - 1) ? total_compacted : compacted_edge_sizes[i+1]) - compacted_edge_sizes[i];
      size_t edges_end = edges_start + num_edges;

      uintE our_id = std::get<0>(merge_seq[sort_start]);
      auto our_edges = edges.cut(edges_start, edges_end);

      size_t offset = compacted_edge_sizes[i];

      for (size_t j=0; j<our_edges.size(); j++) {
        auto [ngh_id, sim] = our_edges[j];
        triples[offset + j] = std::make_tuple(ngh_id, our_id, sim);
      }
    });

    std::cout << "Inserting neighbor endpoints" << std::endl;

    insert_neighbor_endpoints(std::move(triples));

    // Finally, destroy all of the trees incident to newly deactivated vertices.
    parallel_for(0, merge_seq.size(), [&] (size_t i) {
      auto [u, v, wgh] = merge_seq[i];
      auto nghs = std::move(clusters[v].neighbors);
      assert(!clusters[v].active);
      if (nghs.root) {
        if (nghs.root->ref_cnt != 1) std::cout << nghs.root->ref_cnt << std::endl;
        assert(nghs.root->ref_cnt == 1);
      }
    });


    // Sort the weights incident to each merged vertex. This ensures that the
    // merge order stores parallel merges into this vertex in descending order
    // of weight (preventing inversions in the dendrogram).
    parallel_for(0, starts.size(), [&] (size_t i) {
      size_t start = starts[i];
      size_t end = (i == starts.size() - 1) ? merge_seq.size() : starts[i+1];
      auto our_merges = merge_seq.cut(start, end);
      auto comp = [&] (const auto& l, const auto& r) {
        return std::get<2>(l) > std::get<2>(r);  // sort by weights in desc. order
      };
      parlay::sort_inplace(our_merges, comp);
      auto prev_wgh = std::numeric_limits<float>::max();
      for (size_t j=0; j< our_merges.size(); j++) {
        if (std::get<2>(our_merges[j]) > prev_wgh) { std::cout << "sort error!" << std::endl; }
        prev_wgh = std::get<2>(our_merges[j]);
      }
    });

    // Lastly for each merged vertex, save it in the merge_order seq.
    parallel_for(0, merge_seq.size(), [&] (size_t i) {
      assert(merge_order[merge_order_idx + i].first == UINT_E_MAX);
      merge_order[merge_order_idx + i] = std::make_pair(std::get<1>(merge_seq[i]), std::get<2>(merge_seq[i]));
    });
    merge_order_idx += merge_seq.size();

    std::cout << "Finished unite merge." << std::endl;
    um.next("Unite merge time");
  }

  clustered_graph(Graph& G, Weights& weights) : G(G), weights(weights) {
    n = G.n;
    last_cluster_id = n;
    num_merges_performed = 0;
    clusters = sequence<clustered_vertex>(n);
    dendrogram = sequence<std::pair<uintE, W>>(2*n - 1, std::make_pair(UINT_E_MAX, W()));

    neighbor_map::reserve(G.m);

    auto edges = sequence<edge>::uninitialized(G.m);
    auto offsets = sequence<size_t>::from_function(G.n, [&] (size_t i) {
      return G.get_vertex(i).out_degree();
    });
    pbbslib::scan_inplace(make_slice(offsets));
    parallel_for(0, n, [&] (size_t i) {
      size_t off = offsets[i];
      auto map_f = [&] (const uintE& u, const uintE& v, const auto& wgh, size_t j) {
        edges[off + j] = {v, wgh};
      };
      G.get_vertex(i).out_neighbors().map_with_index(map_f);
    });

    parallel_for(0, n, [&] (size_t i) {
      auto orig = G.get_vertex(i);
      size_t off = offsets[i];
      clusters[i] = clustered_vertex(i, orig, weights, edges.begin() + off);
    }, 1);
    std::cout << "Built all vertices" << std::endl;

    merge_order = sequence<std::pair<uintE, float>>(n, std::make_pair(UINT_E_MAX, float()));
    merge_order_idx = 0;
  }

  sequence<std::pair<uintE, double>> get_dendrogram() {
    auto dendrogram = sequence<std::pair<uintE, double>>(2*n-1, std::make_pair(UINT_E_MAX, double()));
    std::cout << "merge order idx = " << merge_order_idx << " n = " << n << std::endl;
    std::cout << "dendrogram length = " << dendrogram.size() << std::endl;

    auto mapping = sequence<uintE>(n, UINT_E_MAX);

    size_t cluster_id = n;  // Next new cluster id is n.
    for (size_t i=0; i<merge_order_idx; i++) {
      auto [u, wgh] = merge_order[i];
      assert(u != UINT_E_MAX);
      uintE merge_target = clusters[u].current_id;

      //std::cout << "Merging " << u << " to " << merge_target << " with weight " << wgh << std::endl;

      // Get a new id.
      uintE new_id = cluster_id;
      cluster_id++;

      assert(u < n);
      if (mapping[u] != UINT_E_MAX) {
        //std::cout << "remapped u to: " << mapping[u] << std::endl;
        u = mapping[u];
      }

      assert(merge_target < n);
      uintE old_merge_target = merge_target;
      if (mapping[merge_target] != UINT_E_MAX) {  // Actually a different (new) cluster.
        //std::cout << "remapped merge_target to: " << mapping[merge_target] << std::endl;
        merge_target = mapping[merge_target];
      }
      mapping[old_merge_target] = new_id;  // Subsequent merge to merge_target uses new_id.
      // std::cout << "new cluster id: " << new_id << std::endl;

      dendrogram[u] = std::make_pair(new_id, wgh);
      dendrogram[merge_target] = std::make_pair(new_id, wgh);
    }


    auto all_vertices = parlay::delayed_seq<uintE>(n, [&] (size_t i) {
      return (clusters[i].current_id == i) ? i : UINT_E_MAX; });
    auto bad = parlay::filter(all_vertices, [&] (uintE u) -> bool { return u != UINT_E_MAX; });

    std::queue<uintE> bad_queue;
    for (size_t i=0; i<bad.size(); i++) {
      bad_queue.push(bad[i]);
    }

    while (bad_queue.size() > 1) {
      uintE fst = bad_queue.front();
      bad_queue.pop();
      uintE snd = bad_queue.front();
      bad_queue.pop();

      if (fst < n && mapping[fst] != UINT_E_MAX) {  // An original vertex.
        fst = mapping[fst];
      }
      if (snd < n && mapping[snd] != UINT_E_MAX) {
        snd = mapping[snd];
      }

      uintE new_id = cluster_id;  // increments next_id
      cluster_id++;
      dendrogram[fst] = {new_id, double(0)};
      dendrogram[snd] = {new_id, double(0)};

      debug(std::cout << "Merged components for: " << fst << " " << snd << " dend_size = " << dendrogram.size() << std::endl);

      bad_queue.push(new_id);
    }

    // Check dendrogram.

    for (size_t i=0; i<(2*n - 2); i++) {
      std::cout << "Checking i = " << i << std::endl;
      cluster_id = i;
      double wgh = std::numeric_limits<double>::max();
      while (true) {
        auto parent = dendrogram[cluster_id].first;
        auto merge_wgh = dendrogram[cluster_id].second;
        std::cout << "id = " << cluster_id << " parent = " << parent << " wgh = " << merge_wgh << std::endl;
        assert(wgh >= merge_wgh);
        wgh = merge_wgh;
        if (cluster_id == parent || parent == UINT_E_MAX) break;
        cluster_id = parent;
      }
      std::cout << "i = " << i << " is good." << std::endl;
    }

    return dendrogram;
  }

};


}  // namespace clustering
}  // namespace gbbs
