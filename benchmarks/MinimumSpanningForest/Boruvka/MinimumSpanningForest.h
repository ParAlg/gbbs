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

#include <cassert>
#include "gbbs/gbbs.h"

#include "benchmarks/SpanningForest/SDB14/SpanningForest.h"

namespace gbbs {
namespace MinimumSpanningForest_boruvka {

constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;

// the number of filtering steps to run
constexpr size_t n_filter_steps = 5;

// Returns edge ids of edges in the mst in the last argument (uintE* mst) and
// returns the number of edges written out as n_in_mst.
template <class W, class M, class P, class D>
inline size_t Boruvka(edge_array<W>& E, uintE*& vtxs, uintE*& next_vtxs,
                      M& min_edges, P& parents, D& exhausted, size_t& n,
                      uintE* mst) {
  using vtxid_wgh_pair = std::pair<uintE, W>;

  using Edge = std::tuple<uintE, uintE, W>;
  size_t m = E.size();
  auto& edges = E.E;
  auto less = [](const vtxid_wgh_pair& a, const vtxid_wgh_pair& b) {
    // returns true if (weight is <) or (weight = and index is <)
    return (a.second < b.second) || (a.second == b.second && a.first < b.first);
  };

  uintE* edge_ids = gbbs::new_array_no_init<uintE>(m);
  parallel_for(0, m, kDefaultGranularity, [&](size_t i) { edge_ids[i] = i; });
  uintE* next_edge_ids = gbbs::new_array_no_init<uintE>(m);

  auto new_mst_edges = sequence<uintE>(n, UINT_E_MAX);
  auto is_root = sequence<bool>(n);

  // Stores edge indices that join the MinimumSpanningForest.
  size_t n_in_mst = 0;
  size_t round = 0;

  while (n > 1 && m > 0) {
    debug(std::cout << "Boruvka round: " << round << " n: " << n << " m: " << m
                    << "\n";);

    timer init_t;
    init_t.start();
    parallel_for(0, n, [&](size_t i) {
      uintE v = vtxs[i];
      // stores an (index, weight) pair
      min_edges[v] =
          std::make_pair(UINT_E_MAX, std::numeric_limits<int32_t>::max());
    });
    init_t.stop();
    debug(init_t.next("init time"););

    // 1. write_min to select the minimum edge out of each component.
    timer min_t;
    min_t.start();
    parallel_for(0, m, kDefaultGranularity, [&](size_t i) {
      uintE e_id = edge_ids[i];
      const Edge& e = edges[e_id];
      vtxid_wgh_pair cas_e(e_id, std::get<2>(e));
      gbbs::write_min(min_edges + std::get<0>(e), cas_e, less);
      gbbs::write_min(min_edges + std::get<1>(e), cas_e, less);
    });
    min_t.stop();
    debug(min_t.next("write min time"););

    // 2. test whether vertices found an edge incident to them
    timer mark_t;
    mark_t.start();
    parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
      uintE v = vtxs[i];
      const auto& e = min_edges[v];
      if (e.first == UINT_E_MAX) {
        // no more edges incident to v in this batch.
        exhausted[v] = true;
        is_root[i] = false;
        new_mst_edges[i] = UINT_E_MAX;
        debug(uintE pv = parents[v]; assert(pv == v););
      } else {
        uintE ind = e.first;
        const auto& edge = edges[ind];
        uintE u = std::get<0>(edge) ^ std::get<1>(edge) ^ v;
        // pick the lower endpoint as the root.
        if (u > v && ind == min_edges[u].first) {
          parents[v] = v;
          is_root[i] = true;
          new_mst_edges[i] = UINT_E_MAX;
        } else {
          // v is satellite: hook onto u.
          parents[v] = u;
          is_root[i] = false;
          new_mst_edges[i] = ind;
        }
      }
    });
    mark_t.stop();
    debug(mark_t.next("mark time"););

    // 3. filter out the new MinimumSpanningForest edges.
    timer filter_t;
    filter_t.start();
    n_in_mst += parlay::filterf(new_mst_edges.begin(), mst + n_in_mst, n,
                                [](uintE v) { return v != UINT_E_MAX; });
    debug(std::cout << "      " << n_in_mst << " edges added to mst."
                    << "\n";);
    filter_t.stop();
    debug(filter_t.next("filter time"););

    // 4. pointer jump to find component centers.
    timer jump_t;
    jump_t.start();
    parallel_for(0, n, [&](size_t i) {
      uintE v = vtxs[i];
      size_t ctr = 0;
      while (parents[v] != parents[parents[v]]) {
        parents[v] = parents[parents[v]];
        ctr++;
      }
    });
    jump_t.stop();
    debug(jump_t.next("jump time"););

    // 5. compact the vertices (pack out the roots)
    timer compact_t;
    compact_t.start();
    auto vtxs_im = gbbs::make_slice<uintE>(vtxs, n);

    n = parlay::pack_out(vtxs_im, is_root, gbbs::make_slice(next_vtxs, m));
    std::swap(vtxs, next_vtxs);
    compact_t.stop();
    debug(compact_t.next("compact time"););
    debug(std::cout << "      " << n << " vertices remain."
                    << "\n";);

    // 6. relabel the edges with the new roots.
    timer relab_t;
    relab_t.start();
    parallel_for(0, m, kDefaultGranularity, [&](size_t i) {
      size_t e_id = edge_ids[i];
      Edge& e = edges[e_id];
      uintE u = std::get<0>(e);
      uintE v = std::get<1>(e);
      uintE pu = parents[std::get<0>(e)];
      uintE pv = parents[std::get<1>(e)];
      if (u != pu || v != pv) {
        W w = std::get<2>(e);
        edges[e_id] = std::make_tuple(pu, pv, w);
      }
      if (pu == pv) {
        edge_ids[i] |= TOP_BIT;
      }
    });
    relab_t.stop();
    debug(relab_t.next("relabel time"););

    // 7. filter (or ignore) self-edges.
    auto self_loop_f = [&](size_t i) { return !(edge_ids[i] & TOP_BIT); };
    auto self_loop_im = parlay::delayed_seq<bool>(n, self_loop_f);
    auto edge_ids_im = gbbs::make_slice(edge_ids, m);
    m = parlay::pack_out(edge_ids_im, self_loop_im,
                         gbbs::make_slice(next_edge_ids, m));

    debug(std::cout << "filter, m is now " << m << " n is now " << n << "\n";);
    std::swap(edge_ids, next_edge_ids);
    round++;
  }

  std::cout << "Boruvka finished: total edges added to MinimumSpanningForest = "
            << n_in_mst << "\n";
  gbbs::free_array(edge_ids, m);
  gbbs::free_array(next_edge_ids, m);
  return n_in_mst;
  //  auto mst_im = sequence<uintE>(mst, n_in_mst); // allocated
  //  return mst_im;
}

constexpr size_t sample_size = 2000;
inline size_t hash_to_range(size_t hsh, size_t range) { return hsh & range; }

inline size_t key_for_pair(uint32_t k1, uintE k2, parlay::random rnd) {
  size_t key = (static_cast<size_t>(k1) << 32) + static_cast<size_t>(k2);
  return rnd.ith_rand(key);
}

template <template <class W> class vertex, class W>
inline edge_array<W> get_all_edges(symmetric_graph<vertex, W>& G) {
  auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return true;
  };
  return filter_all_edges(G, pred);
}

template <template <class W> class vertex, class W>
inline edge_array<W> get_top_k(symmetric_graph<vertex, W>& G, size_t k,
                               parlay::random r, bool first_round = false) {
  if (k == static_cast<size_t>(G.m)) {
    return get_all_edges(G);
  }
  using edge = std::tuple<uintE, uintE, W>;
  timer st;
  st.start();

  size_t n = G.n;
  size_t m = G.m;

  auto vertex_offs = sequence<long>(G.n);
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    vertex_offs[i] = G.get_vertex(i).out_degree();
  });
  parlay::scan_inclusive_inplace(make_slice(vertex_offs));

  auto sample_edges = sequence<edge>(sample_size);
  auto lte = [&](const size_t& left, const size_t& right) {
    return left <= right;
  };

  parallel_for(0, sample_size, kDefaultGranularity, [&](size_t i) {
    size_t sample_edge = r.ith_rand(i) % m;
    uintE vtx = parlay::binary_search(vertex_offs, sample_edge, lte);
    size_t ith = vertex_offs[vtx] - sample_edge - 1;
    uintE ngh;
    W wgh;
    std::tie(ngh, wgh) =
        G.get_vertex(vtx).out_neighbors().get_ith_neighbor(ith);
    sample_edges[i] = std::make_tuple(vtx, ngh, wgh);
  });

  auto cmp_by_wgh = [](const edge& left, const edge& right) {
    return std::get<2>(left) < std::get<2>(right);
  };
  parlay::sample_sort_inplace(make_slice(sample_edges), cmp_by_wgh);

  // 2. find approximate splitter.
  size_t ind = ((double)(k * sample_edges.size())) / G.m;
  auto splitter = sample_edges[ind];
  W split_weight = std::get<2>(splitter);

  size_t first_ind = 0;
  size_t last_ind = 0;
  size_t ssize = sample_edges.size();
  parallel_for(0, ssize, kDefaultGranularity, [&](size_t i) {
    if (std::get<2>(sample_edges[i]) == split_weight) {
      if (i == 0 || (std::get<2>(sample_edges[i - 1]) != split_weight)) {
        first_ind = i;
      }
      if (i == (ssize - 1) ||
          (std::get<2>(sample_edges[i + 1]) != split_weight)) {
        last_ind = i;
      }
    }
  });

  size_t weight_size = (last_ind - first_ind + 1);
  double split_wgh_fraction = ((1.0 * (last_ind - first_ind + 1)) / ssize);
  debug(std::cout << "split wgh is: " << split_weight << "\n";
        std::cout << "fraction of sample composed by split_wgh = "
                  << split_wgh_fraction << "\n";);
  st.stop();
  debug(st.next("startup time"););

  // 3. filter edges based on splitter
  if (split_wgh_fraction < 0.2) {
    auto filter_pred = [&](const uintE& src, const uintE& ngh,
                           const W& wgh) -> int {
      if (src > ngh) return 1;  // filter out
      if (wgh <= split_weight) {
        return 2;  // return in array
      }
      return 0;
    };
    return filter_edges(G, filter_pred);
  } else {
    // Special pred for performing extra hashing on edges with weight ==
    // split_wgh.
    double frac_to_take = (1.0 * (ind - first_ind)) / weight_size;
    std::cout << "frac of split_weight to take: " << frac_to_take << "\n";
    size_t range = (1L << parlay::log2_up(G.m)) - 1;
    size_t threshold = frac_to_take * range;
    // account for filtering directed edges
    threshold *= (first_round ? 2.0 : 1.0);
    auto filter_split_wgh_pred = [&](const uintE& src, const uintE& ngh,
                                     const W& wgh) -> int {
      if (src > ngh) return 1;
      if (wgh < split_weight) {
        return 2;  // return in array
      } else if (wgh == split_weight) {
        if (hash_to_range(key_for_pair(src, ngh, r), range) < threshold)
          return 2;
      }
      return 0;
    };
    return filter_edges(G, filter_split_wgh_pred);
  }
}

template <template <class W> class vertex, class W,
          typename std::enable_if<!std::is_same<W, gbbs::empty>::value,
                                  int>::type = 0>
inline sequence<std::tuple<uintE, uintE, W>> MinimumSpanningForest(
    symmetric_graph<vertex, W>& GA, bool largemem = false) {
  using edge = std::tuple<uintE, uintE, W>;
  using vtxid_wgh_pair = std::pair<uintE, W>;

  size_t n = GA.n;
  std::cout << "n = " << n << "\n";
  auto r = parlay::random();

  auto exhausted =
      sequence<bool>::from_function(n, [](size_t i) { return false; });
  auto parents = sequence<uintE>::from_function(n, [](size_t i) { return i; });
  auto mst_edges = parlay::sequence<edge>();

  auto min_edges = gbbs::new_array_no_init<vtxid_wgh_pair>(n);

  size_t n_active = n;
  uintE* vtxs = gbbs::new_array_no_init<uintE>(n_active);
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { vtxs[i] = i; });
  uintE* next_vtxs = gbbs::new_array_no_init<uintE>(n_active);

  //  for (size_t i=0; i<n; i++) {
  //    auto v = GA.get_vertex(i);
  //    auto map_edges = [&] (const uintE& u, const uintE& v, const W& wgh) {
  //      std::cout << u << " " << v << " " << wgh << std::endl;
  //    };
  //    v.out_neighbors().map(map_edges, false);
  //  }

  size_t round = 0;
  while (GA.m > 0) {
    timer round_t;
    round_t.start();
    std::cout << "\n";
    std::cout << "round = " << round << " n_active = " << n_active
              << " GA.m = " << GA.m
              << " MinimumSpanningForest size = " << mst_edges.size() << "\n";

    // find a prefix of lowest-weight edges.
    size_t split_idx = std::min((3 * n) / 2, (size_t)GA.m);
    if (largemem) {
      split_idx = (round == 0) ? std::min((9 * n) / 8, (size_t)GA.m)
                               : std::min(n / 2, (size_t)GA.m);
    }

    timer get_t;
    get_t.start();
    auto E = (round < n_filter_steps) ? get_top_k(GA, split_idx, r, round == 0)
                                      : get_all_edges(GA);
    get_t.stop();
    debug(get_t.next("get time"););
    size_t n_edges = E.size();
    std::cout << "Prefix size = " << split_idx << " #edges = " << n_edges
              << " G.m is now = " << GA.m << "\n";

    // relabel edges
    auto& edges = E.E;
    auto edges_save = E.E;  // copy
    if (round > 0) {
      parallel_for(0, n_edges, kDefaultGranularity, [&](size_t i) {
        edge& e = edges[i];
        uintE u = std::get<0>(e);
        uintE v = std::get<1>(e);
        std::get<0>(e) = parents[u];
        std::get<1>(e) = parents[v];
      });
    }

    // run Boruvka on the prefix and add new edges to mst_edges
    timer bt;
    bt.start();
    uintE* mst = gbbs::new_array_no_init<uintE>(n);
    size_t n_in_mst = Boruvka(E, vtxs, next_vtxs, min_edges, parents, exhausted,
                              n_active, mst);
    auto edge_ids = gbbs::make_slice(mst, n_in_mst);
    bt.stop();
    debug(bt.next("boruvka time"););

    auto edges_to_add = parlay::delayed_seq<edge>(edge_ids.size(), [&](size_t i) { return edges_save[edge_ids[i]]; });
    mst_edges.append(edges_to_add);
    edges_save.clear();
    gbbs::free_array(mst, n);

    // reactivate vertices and reset exhausted
    timer pack_t;
    pack_t.start();

    auto vtx_range = gbbs::make_slice(vtxs + n_active, vtxs + n);
    n_active += parlay::pack_index_out(make_slice(exhausted), vtx_range);
    pack_t.stop();
    debug(pack_t.next("reactivation pack"););

    parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
      if (exhausted[i]) exhausted[i] = false;
    });

    // pointer jump: vertices that were made inactive could have had their
    // parents change.
    parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
      size_t ctr = 0;
      while (parents[i] != parents[parents[i]]) {
        parents[i] = parents[parents[i]];
        ctr++;
      }
    });

    // pack out all edges in the graph that are shortcut by the added edges
    auto filter_pred = [&](const uintE& src, const uintE& ngh,
                           const W& wgh) -> int {
      auto c_src = parents[src];
      auto c_ngh = parents[ngh];
      return c_src == c_ngh;
    };
    debug(std::cout << "Filtering G, m = " << GA.m << "\n";);
    timer filter_t;
    filter_t.start();
    filter_edges(GA, filter_pred);
    filter_t.stop();
    filter_t.next("filter time");
    std::cout << "After filter, m is now " << GA.m << "\n";
    round++;
    round_t.stop();
    round_t.next("round time");

    r = r.next();
  }
  std::cout << "#edges in output mst: " << mst_edges.size() << "\n";
  auto wgh_imap_f = [&](size_t i) { return std::get<2>(mst_edges[i]); };
  auto wgh_imap = parlay::delayed_seq<size_t>(mst_edges.size(), wgh_imap_f);
  std::cout << "total weight = " << parlay::reduce(wgh_imap) << "\n";

  gbbs::free_array(min_edges, n);
  return mst_edges;
}

template <
    template <class W> class vertex, class W,
    typename std::enable_if<std::is_same<W, gbbs::empty>::value, int>::type = 0>
inline sequence<std::pair<uintE, uintE>> MinimumSpanningForest(
    symmetric_graph<vertex, W>& GA) {
  return workefficient_sf::SpanningForest(GA);
}

}  // namespace MinimumSpanningForest_boruvka
}  // namespace gbbs
