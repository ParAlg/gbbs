#pragma once

#include "gbbs/bucket.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"
#include "gbbs/helpers/sparse_table.h"

#include "induced_hybrid.h"
#include "intersect.h"

namespace gbbs {

// Wrapper for a hash function
struct hashtup {
  inline size_t operator()(const uintE& a) const { return parlay::hash64_2(a); }
};

// For triangle peeling
template <class Graph, class Graph2, class F, class H>
void vtx_intersect(Graph& G, Graph2& DG, F f, H ignore_f, uintE vg, uintE vdg) {
  // Intersect the neighbors of vg (vertex in G) with the neighbros of vdg
  // (vertex in DG)
  size_t v_deg = G.get_vertex(vg).out_degree();
  size_t i_deg = DG.get_vertex(vdg).out_degree();
  auto i_iter = DG.get_vertex(vdg).out_neighbors().get_iter();
  auto v_iter = G.get_vertex(vg).out_neighbors().get_iter();
  size_t i_iter_idx = 0;
  size_t v_iter_idx = 0;
  // Count of number of neighbors in the intersection
  size_t j = 0;
  while (i_iter_idx < i_deg && v_iter_idx < v_deg) {
    if (std::get<0>(i_iter.cur()) == std::get<0>(v_iter.cur())) {
      auto nbhr = std::get<0>(i_iter.cur());
      if (ignore_f(vg, nbhr) && ignore_f(vdg, nbhr)) {
        // If the neighbor is valid (ignore_f), one triangle is added to the
        // neighbor
        f(nbhr, 1);
        j++;
      }
      i_iter_idx++;
      v_iter_idx++;
      if (i_iter.has_next()) i_iter.next();
      if (v_iter.has_next()) v_iter.next();
    } else if (std::get<0>(i_iter.cur()) < std::get<0>(v_iter.cur())) {
      i_iter_idx++;
      if (i_iter.has_next()) i_iter.next();
    } else {
      v_iter_idx++;
      if (v_iter.has_next()) v_iter.next();
    }
  }
  // j triangles are added to the neighbor vdg
  f(vdg, j);
}

template <class Graph, class Graph2, class F, class I>
inline size_t triUpdate_serial(Graph& G, Graph2& DG, F get_active,
                               size_t active_size, char* still_active,
                               sequence<uintE>& rank,
                               sequence<size_t>& per_processor_counts,
                               bool do_update_changed, I update_changed,
                               sequence<uintE>& update_idxs) {
  using W = typename Graph::weight_type;

  // Mark every vertex in the active set
  parallel_for(0, active_size,
               [&](size_t j) { still_active[get_active(j)] = 1; }, 2048);

  size_t num_updates = 0;

  // Function that dictates which edges to consider in first level of recursion
  auto ignore_f = [&](const uintE& u, const uintE& v) {
    auto status_u = still_active[u];
    auto status_v = still_active[v];
    if (status_u == 2 || status_v == 2) return false;  // deleted edge
    if (status_v == 0) return true;  // non-deleted, non-active edge
    return rank[u] < rank[v];        // orient edges if in active set
  };

  // Collate triangle counts by processor
  auto update_d = [&](uintE vtx, size_t count) {
    size_t ct = per_processor_counts[vtx];
    per_processor_counts[vtx] += count;
    if (ct == 0) {
      update_idxs[num_updates] = vtx;
      num_updates++;
    }
  };

  // Triangle count updates
  for (size_t i = 0; i < active_size; i++) {
    auto j = get_active(i);
    // For each neighbor v of active vertex j = u, intersect N(v) with N(j)
    auto map_label_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (ignore_f(u, v)) vtx_intersect(G, DG, update_d, ignore_f, u, v);
    };
    G.get_vertex(j).out_neighbors().map(map_label_f, false);
  };

  // Perform update_changed on each vertex with changed triangle counts
  if (do_update_changed) {
    for (size_t i = 0; i < num_updates; i++) {
      update_changed(per_processor_counts, i, update_idxs[i]);
    }
  }

  // Mark every vertex in the active set as deleted
  parallel_for(0, active_size,
               [&](size_t j) { still_active[get_active(j)] = 2; }, 2048);

  return num_updates;
}

template <class Graph, class Graph2, class F, class H, class I>
inline size_t triUpdate(Graph& G, Graph2& DG, F get_active, size_t active_size,
                        size_t granularity, char* still_active,
                        sequence<uintE>& rank,
                        sequence<size_t>& per_processor_counts, H update,
                        bool do_update_changed, I update_changed) {
  using W = typename Graph::weight_type;
  auto n = G.n;

  // Mark every vertex in the active set
  parallel_for(0, active_size,
               [&](size_t j) { still_active[get_active(j)] = 1; }, 2048);

  // Sum of all out-degrees in active set
  size_t active_deg = 0;
  auto degree_map = parlay::delayed_seq<size_t>(active_size, [&](size_t i) {
    return G.get_vertex(get_active(i)).out_degree();
  });
  active_deg += parlay::reduce(degree_map);

  // Hash table to contain triangle count updates
  size_t edge_table_size = (size_t)(active_deg < n ? active_deg : n);
  auto edge_table = gbbs::sparse_table<uintE, bool, hashtup>(
      edge_table_size, std::make_tuple(UINT_E_MAX, false), hashtup());

  // Function that dictates which edges to consider in first level of recursion
  auto ignore_f = [&](const uintE& u, const uintE& v) {
    auto status_u = still_active[u];
    auto status_v = still_active[v];
    if (status_u == 2 || status_v == 2) return false;  // deleted edge
    if (status_v == 0) return true;  // non-deleted, non-active edge
    return rank[u] < rank[v];        // orient edges if in active set
  };

  // Collate triangle counts by processor
  auto update_d = [&](uintE vtx, size_t count) {
    size_t worker = worker_id();
    size_t ct = per_processor_counts[worker * n + vtx];
    per_processor_counts[worker * n + vtx] += count;
    if (ct == 0) edge_table.insert(std::make_tuple(vtx, true));
  };

  // Triangle count updates
  parallel_for(0, active_size,
               [&](size_t i) {
                 auto j = get_active(i);
                 // For each neighbor v of active vertex j = u, intersect N(v)
                 // with N(j)
                 auto map_label_f = [&](const uintE& u, const uintE& v,
                                        const W& wgh) {
                   if (ignore_f(u, v))
                     vtx_intersect(G, DG, update_d, ignore_f, u, v);
                 };
                 G.get_vertex(j).out_neighbors().map(map_label_f, false);
               },
               granularity, false);

  // Extract all vertices with changed triangle counts
  auto changed_vtxs = edge_table.entries();

  // Aggregate the updated counts across all worker's local arrays, as specified
  // by update
  parallel_for(0, changed_vtxs.size(),
               [&](size_t i) {
                 size_t nthreads = num_workers();
                 uintE v = std::get<0>(changed_vtxs[i]);
                 for (size_t j = 0; j < nthreads; j++) {
                   update(per_processor_counts, j, v);
                 }
               },
               128);

  // Perform update_changed on each vertex with changed triangle counts
  if (do_update_changed) {
    parallel_for(0, changed_vtxs.size(), [&](size_t i) {
      update_changed(per_processor_counts, i, std::get<0>(changed_vtxs[i]));
    });
  }

  // Mark every vertex in the active set as deleted
  parallel_for(0, active_size,
               [&](size_t j) { still_active[get_active(j)] = 2; }, 2048);

  return changed_vtxs.size();
}

template <class Graph, class Graph2, class F, class H, class I>
inline size_t cliqueUpdate(Graph& G, Graph2& DG, size_t k, size_t max_deg,
                           bool label, F get_active, size_t active_size,
                           size_t granularity, char* still_active,
                           sequence<uintE>& rank,
                           sequence<size_t>& per_processor_counts, H update,
                           bool do_update_changed, I update_changed) {
  const size_t n = G.n;

  // Set up space for clique counting
  auto init_induced = [&](HybridSpace_lw* induced) {
    induced->alloc(max_deg, k, n, label, true);
  };
  auto finish_induced = [&](HybridSpace_lw* induced) {
    if (induced != nullptr) {
      delete induced;
    }
  };

  // Sum of all out-degrees in active set
  size_t active_deg = 0;
  auto degree_map = parlay::delayed_seq<size_t>(active_size, [&](size_t i) {
    return G.get_vertex(get_active(i)).out_degree();
  });
  active_deg += parlay::reduce(degree_map);

  // Mark every vertex in the active set
  parallel_for(0, active_size,
               [&](size_t j) { still_active[get_active(j)] = 1; }, 2048);

  // Hash table to contain clique count updates
  size_t edge_table_size = (size_t)(active_deg < n ? active_deg : n);
  auto edge_table = gbbs::sparse_table<uintE, bool, hashtup>(
      edge_table_size, std::make_tuple(UINT_E_MAX, false), hashtup());

  // Function that dictates which edges to consider in first level of recursion
  auto ignore_f = [&](const uintE& u, const uintE& v) {
    auto status_u = still_active[u];
    auto status_v = still_active[v];
    if (status_u == 2 || status_v == 2) return false;  // deleted edge
    if (status_v == 0) return true;  // non-deleted, non-active edge
    return rank[u] < rank[v];        // orient edges if in active set
  };

  // Collate clique counts by processor
  auto update_d = [&](uintE vtx, size_t count) {
    size_t worker = worker_id();
    size_t ct = per_processor_counts[worker * n + vtx];
    per_processor_counts[worker * n + vtx] += count;
    if (ct == 0) edge_table.insert(std::make_tuple(vtx, true));
  };

  // Clique count updates
  parallel_for_alloc<HybridSpace_lw>(
      init_induced, finish_induced, 0, active_size,
      [&](size_t i, HybridSpace_lw* induced) {
        auto vert = get_active(i);
        if (G.get_vertex(vert).out_degree() != 0) {
          induced->setup(G, DG, k, vert, ignore_f);
          induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced,
                                                     update_d);
        }
      },
      granularity, false);

  // Extract all vertices with changed clique counts
  auto changed_vtxs = edge_table.entries();

  // Aggregate the updated counts across all worker's local arrays, as specified
  // by update
  parallel_for(0, changed_vtxs.size(),
               [&](size_t i) {
                 size_t nthreads = num_workers();
                 uintE v = std::get<0>(changed_vtxs[i]);
                 for (size_t j = 0; j < nthreads; j++) {
                   update(per_processor_counts, j, v);
                 }
               },
               128);

  // Perform update_changed on each vertex with changed clique counts
  if (do_update_changed) {
    parallel_for(0, changed_vtxs.size(), [&](size_t i) {
      update_changed(per_processor_counts, i, std::get<0>(changed_vtxs[i]));
    });
  }

  // Mark every vertex in the active set as deleted
  parallel_for(0, active_size,
               [&](size_t j) { still_active[get_active(j)] = 2; }, 2048);

  return changed_vtxs.size();
}

template <class Graph, class Graph2, class F, class I>
inline size_t cliqueUpdate_serial(Graph& G, Graph2& DG, size_t k,
                                  size_t max_deg, bool label, F get_active,
                                  size_t active_size, char* still_active,
                                  sequence<uintE>& rank,
                                  sequence<size_t>& per_processor_counts,
                                  bool do_update_changed, I update_changed,
                                  sequence<uintE>& update_idxs) {
  const size_t n = G.n;

  // Mark every vertex in the active set
  parallel_for(0, active_size,
               [&](size_t j) { still_active[get_active(j)] = 1; }, 2048);

  size_t num_updates = 0;

  // Set up space for clique counting
  HybridSpace_lw* induced = new HybridSpace_lw();
  induced->alloc(max_deg, k, n, label, true);

  // Function that dictates which edges to consider in first level of recursion
  auto ignore_f = [&](const uintE& u, const uintE& v) {
    auto status_u = still_active[u];
    auto status_v = still_active[v];
    if (status_u == 2 || status_v == 2) return false;  // deleted edge
    if (status_v == 0) return true;  // non-deleted, non-active edge
    return rank[u] < rank[v];        // orient edges if in active set
  };

  // Collate clique counts by processor
  auto update_d = [&](uintE vtx, size_t count) {
    size_t ct = per_processor_counts[vtx];
    per_processor_counts[vtx] += count;
    if (ct == 0) {
      update_idxs[num_updates] = vtx;
      num_updates++;
    }
  };

  // Clique count updates
  for (size_t i = 0; i < active_size; i++) {
    auto vtx = get_active(i);
    if (G.get_vertex(vtx).out_degree() != 0) {
      induced->setup(G, DG, k, vtx, ignore_f);
      induced_hybrid::KCliqueDir_fast_hybrid_rec(G, 1, k, induced, update_d);
    }
  }

  // Clean up space for clique counting
  if (induced != nullptr) {
    delete induced;
  }

  // Perform update_changed on each vertex with changed clique counts
  if (do_update_changed) {
    for (size_t i = 0; i < num_updates; i++) {
      update_changed(per_processor_counts, i, update_idxs[i]);
    }
  }

  // Mark every vertex in the active set as deleted
  parallel_for(0, active_size,
               [&](size_t j) { still_active[get_active(j)] = 2; }, 2048);

  return num_updates;
}

template <typename bucket_t, class Graph, class Graph2>
sequence<bucket_t> Peel(Graph& G, Graph2& DG, size_t k, size_t* cliques,
                        bool label, sequence<uintE>& rank,
                        size_t num_buckets = 16) {
  timer t2;
  t2.start();
  const size_t n = G.n;
  auto D = sequence<bucket_t>::from_function(
      n, [&](size_t i) { return cliques[i]; });
  auto D_filter = sequence<std::tuple<uintE, bucket_t>>::uninitialized(n);

  auto b = make_vertex_custom_buckets<bucket_t>(n, D, increasing, num_buckets);

  auto per_processor_counts =
      sequence<size_t>(n * num_workers(), static_cast<size_t>(0));

  char* still_active = (char*)calloc(n, sizeof(char));
  size_t max_deg =
      induced_hybrid::get_max_deg(G);  // could instead do max_deg of active?
  auto update_idxs = sequence<uintE>(max_deg);

  size_t rounds = 0;
  size_t finished = 0;
  bucket_t cur_bkt = 0;
  bucket_t max_bkt = 0;
  double max_density = 0;
  bool use_max_density = false;
  // Peel each bucket
  auto update_clique = [&](sequence<size_t>& ppc, size_t j, uintE v) {
    if (j == 0) return;
    ppc[v] += ppc[j * n + v];
    ppc[j * n + v] = 0;
  };
  while (finished != n) {
    // for max_density
    if (use_max_density) {
      auto degree_f = [&](size_t i) { return cliques[i]; };
      auto degree_seq = parlay::delayed_seq<size_t>(n, degree_f);
      auto edges_remaining = parlay::reduce(degree_seq);
      auto vtxs_remaining = n - finished;
      double current_density =
          ((double)edges_remaining) / ((double)vtxs_remaining);
      if (current_density > max_density) max_density = current_density;
    }

    // Retrieve next bucket
    auto bkt = b.next_bucket();
    auto active = vertexSubset(n, std::move(bkt.identifiers));
    cur_bkt = bkt.id;

    finished += active.size();
    max_bkt = std::max(cur_bkt, max_bkt);

    auto get_active = [&](size_t j) { return active.vtx(j); };

    size_t granularity = (cur_bkt * active.size() < 10000) ? 1024 : 1;

    size_t filter_size = 0;

    if (active.size() > 1) {
      auto update_changed = [&](sequence<size_t>& ppc, size_t i, uintE v) {
        /* Update the clique count for v, and zero out first worker's count */
        cliques[v] -= ppc[v];
        ppc[v] = 0;
        bucket_t deg = D[v];
        if (deg > cur_bkt) {
          bucket_t new_deg = std::max((bucket_t)cliques[v], (bucket_t)cur_bkt);
          D[v] = new_deg;
          // store (v, bkt) in an array now, pass it to apply_f below instead of
          // what's there right now -- maybe just store it in D_filter?
          D_filter[i] = std::make_tuple(v, b.get_bucket(deg, new_deg));
        } else
          D_filter[i] = std::make_tuple(UINT_E_MAX, 0);
      };
      if (k == 2)
        filter_size = triUpdate(G, DG, get_active, active.size(), granularity,
                                still_active, rank, per_processor_counts,
                                update_clique, true, update_changed);
      else
        filter_size =
            cliqueUpdate(G, DG, k, max_deg, label, get_active, active.size(),
                         granularity, still_active, rank, per_processor_counts,
                         update_clique, true, update_changed);
    } else {
      size_t filter_size_serial = 0;
      auto update_changed_serial = [&](sequence<size_t>& ppc, size_t i,
                                       uintE v) {
        /* Update the clique count for v, and zero out first worker's count */
        cliques[v] -= ppc[v];
        ppc[v] = 0;
        bucket_t deg = D[v];
        if (deg > cur_bkt) {
          bucket_t new_deg = std::max((bucket_t)cliques[v], (bucket_t)cur_bkt);
          D[v] = new_deg;
          D_filter[filter_size_serial] =
              std::make_tuple(v, b.get_bucket(deg, new_deg));
          filter_size_serial++;
        }
      };
      if (k == 2)
        triUpdate_serial(G, DG, get_active, active.size(), still_active, rank,
                         per_processor_counts, true, update_changed_serial,
                         update_idxs);
      else
        cliqueUpdate_serial(G, DG, k, max_deg, label, get_active, active.size(),
                            still_active, rank, per_processor_counts, true,
                            update_changed_serial, update_idxs);
      filter_size = filter_size_serial;
    }

    auto apply_f = [&](size_t i) -> std::optional<std::tuple<uintE, bucket_t>> {
      uintE v = std::get<0>(D_filter[i]);
      bucket_t bucket = std::get<1>(D_filter[i]);
      if (v != UINT_E_MAX && still_active[v] != 2) return wrap(v, bucket);
      return std::nullopt;
    };

    b.update_buckets(apply_f, filter_size);

    parallel_for(0, active.size(),
                 [&](size_t j) { cliques[active.vtx(j)] = 0; }, 2048);

    rounds++;
  }

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  std::cout.precision(17);
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "clique core: " << max_bkt << std::endl;
  if (use_max_density) std::cout << "max density: " << max_density << std::endl;

  free(still_active);

  return D;
}

template <class Graph, class Graph2>
double ApproxPeel(Graph& G, Graph2& DG, size_t k, size_t* cliques,
                  size_t num_cliques, bool label, sequence<uintE>& rank,
                  double eps) {
  std::cout << "eps: " << eps << "\n";
  timer t2;
  t2.start();
  const size_t n = G.n;
  auto D =
      sequence<size_t>::from_function(n, [&](size_t i) { return cliques[i]; });
  auto vertices_remaining =
      parlay::delayed_seq<uintE>(n, [&](size_t i) { return i; });

  size_t round = 1;
  sequence<uintE> last_arr;
  size_t remaining_offset = 0;
  size_t num_vertices_remaining = n;
  double density_multiplier = (k + 1) * (1. + eps);

  double max_density = 0.0;
  char* still_active = (char*)calloc(n, sizeof(char));
  size_t max_deg = induced_hybrid::get_max_deg(G);
  auto per_processor_counts =
      sequence<size_t>(n * num_workers(), static_cast<size_t>(0));

  auto update_clique = [&](sequence<size_t>& ppc, size_t j, uintE v) {
    D[v] -= ppc[j * n + v];
    ppc[j * n + v] = 0;
  };
  auto nop = [&](sequence<size_t>& ppc, size_t i, uintE v) { return; };

  // First round
  {
    size_t edges_remaining = num_cliques;
    // Update density
    double current_density = ((double)edges_remaining) / ((double)(n));
    double target_density = (density_multiplier * ((double)edges_remaining)) /
                            ((double)vertices_remaining.size());
    auto rho = target_density;
    if (current_density > max_density) max_density = current_density;

    auto keep_seq = parlay::delayed_seq<bool>(
        n, [&](size_t i) { return !(D[i] <= target_density); });

    sequence<uintE> this_arr;
    size_t num_removed;
    std::tie(this_arr, num_removed) =
        parlay::split_two(vertices_remaining, keep_seq);
    size_t active_size = num_removed;

    // remove this_arr vertices ************************************************
    size_t granularity = (rho * active_size < 10000) ? 1024 : 1;
    auto get_active = [&](size_t j) { return this_arr[j]; };
    if (k == 2)
      triUpdate(G, DG, get_active, active_size, granularity, still_active, rank,
                per_processor_counts, update_clique, false, nop);
    else
      cliqueUpdate(G, DG, k, max_deg, label, get_active, active_size,
                   granularity, still_active, rank, per_processor_counts,
                   update_clique, false, nop);
    parallel_for(0, active_size, [&](size_t j) { D[this_arr[j]] = 0; }, 2048);
    //***************

    round++;
    last_arr = std::move(this_arr);
    remaining_offset = num_removed;
    num_vertices_remaining -= num_removed;
  }

  while (num_vertices_remaining > 0) {
    uintE* start = last_arr.begin() + remaining_offset;
    uintE* end = start + num_vertices_remaining;
    auto vtxs_remaining = gbbs::make_slice(start, end);

    auto degree_f = [&](size_t i) {
      uintE v = vtxs_remaining[i];
      return static_cast<size_t>(D[v]);
    };
    auto degree_seq =
        parlay::delayed_seq<size_t>(vtxs_remaining.size(), degree_f);
    long edges_remaining = parlay::reduce(degree_seq);

    // Update density
    double current_density =
        ((double)edges_remaining) / ((double)vtxs_remaining.size());
    double target_density = (density_multiplier * ((double)edges_remaining)) /
                            ((double)vtxs_remaining.size());
    auto rho = target_density;
    if (current_density > max_density) max_density = current_density;

    auto keep_seq = parlay::delayed_seq<bool>(
        vtxs_remaining.size(),
        [&](size_t i) { return !(D[vtxs_remaining[i]] <= target_density); });

    sequence<uintE> this_arr;
    size_t num_removed;
    std::tie(this_arr, num_removed) =
        parlay::split_two(vtxs_remaining, keep_seq);

    num_vertices_remaining -= num_removed;
    if (num_vertices_remaining > 0) {
      size_t active_size = num_removed;

      // remove this_arr vertices
      // ************************************************
      size_t granularity = (rho * active_size < 10000) ? 1024 : 1;
      auto get_active = [&](size_t j) { return this_arr[j]; };
      if (k == 2)
        triUpdate(G, DG, get_active, active_size, granularity, still_active,
                  rank, per_processor_counts, update_clique, false, nop);
      else
        cliqueUpdate(G, DG, k, max_deg, label, get_active, active_size,
                     granularity, still_active, rank, per_processor_counts,
                     update_clique, false, nop);
      parallel_for(0, active_size, [&](size_t j) { D[this_arr[j]] = 0; }, 2048);
      //***************
    }

    round++;
    last_arr = std::move(this_arr);
    remaining_offset = num_removed;
  }

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;
  std::cout << "rho: " << round << std::endl;
  std::cout << "### Density of (2(1+\eps))-Densest Subgraph is: " << max_density
            << std::endl;

  free(still_active);

  return max_density;
}

}  // namespace gbbs
