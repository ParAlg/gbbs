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

#include "ligra.h"
#include "union_find.h"

#include "lib/binary_search.h"
#include "lib/dyn_arr.h"
#include "lib/random.h"
#include "lib/sample_sort.h"
#include "lib/speculative_for.h"

namespace MST_boruvka {

constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;

// the number of filtering steps to run
constexpr size_t n_filter_steps = 5;

struct cas_type {
  uintE index;
  int32_t weight;
  cas_type() {
    index = UINT_E_MAX;
    weight = std::numeric_limits<int32_t>::max();
  }
  cas_type(uintE _index, int32_t _w) : index(_index), weight(_w) {}
};

template <class W, class M, class P, class D>
inline array_imap<uintE> Boruvka(edge_array<W>& E, uintE*& vtxs,
                                 uintE*& next_vtxs, M& min_edges, P& parents,
                                 D& exhausted, size_t& n) {
  using ct = cas_type;
  using edge = std::tuple<uintE, uintE, W>;
  size_t m = E.non_zeros;
  auto edges = E.E;
  auto less = [](const ct& a, const ct& b) {
    return (a.weight < b.weight) || (a.weight == b.weight && a.index < b.index);
  };

  uintE* edge_ids = newA(uintE, m);
  parallel_for_bc(i, 0, m, (m > pbbs::kSequentialForThreshold),
                  { edge_ids[i] = i; });
  uintE* next_edge_ids = nullptr;

  auto new_mst_edges = array_imap<uintE>(n, UINT_E_MAX);
  auto is_root = array_imap<bool>(n);

  uintE* mst = newA(uintE, n);
  size_t n_in_mst = 0;
  size_t round = 0;

  while (n > 1 && m > 0) {
    cout << "Boruvka round: " << round << " n: " << n << " m: " << m << endl;

    timer init_t;
    init_t.start();
    parallel_for_bc(i, 0, n, (n > 2000), {
      uintE v = vtxs[i];
      min_edges[v] = ct();
    });
    init_t.stop();  // init_t.reportTotal("init time");

    // 1. writeMin to select the minimum edge out of each component.
    timer min_t;
    min_t.start();
    parallel_for_bc(i, 0, m, (m > pbbs::kSequentialForThreshold), {
      uintE e_id = edge_ids[i];
      const edge& e = edges[e_id];
      ct cas_e(e_id, std::get<2>(e));
      writeMin(min_edges + std::get<0>(e), cas_e, less);
      writeMin(min_edges + std::get<1>(e), cas_e, less);
    });
    min_t.stop();  // min_t.reportTotal("write min time");

    // 2. test whether vertices found an edge incident to them
    timer mark_t;
    mark_t.start();
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      uintE v = vtxs[i];
      const auto& e = min_edges[v];
      if (e.index == UINT_E_MAX) {
        // no more edges incident to v in this batch.
        exhausted[v] = true;
        is_root[i] = false;
        new_mst_edges[i] = UINT_E_MAX;
        assert(parents[v] == v);
      } else {
        uintE ind = e.index;
        const auto& edge = edges[ind];
        uintE u = std::get<0>(edge) ^ std::get<1>(edge) ^ v;
        // pick the higher endpoint as the root.
        if (u < v && ind == min_edges[u].index) {
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
    mark_t.stop();  // mark_t.reportTotal("mark time");

    // 3. filter out the new MST edges.
    timer filter_t;
    filter_t.start();
    n_in_mst += pbbs::filterf(new_mst_edges.start(), mst + n_in_mst, n,
                              [](uintE v) { return v != UINT_E_MAX; });
    cout << "      " << n_in_mst << " edges added to mst." << endl;
    filter_t.stop();  // filter_t.reportTotal("filter time");

    // 4. pointer jump to find component centers.
    timer jump_t;
    jump_t.start();
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      uintE v = vtxs[i];
      size_t ctr = 0;
      while (parents[v] != parents[parents[v]]) {
        parents[v] = parents[parents[v]];
        ctr++;
      }
    });
    jump_t.stop();  // jump_t.reportTotal("jump time");

    // 5. compact the vertices (pack out the roots)
    timer compact_t;
    compact_t.start();
    auto vtxs_im = make_array_imap(vtxs, n);
    auto V_im = pbbs::pack(vtxs_im, is_root, pbbs::no_flag, next_vtxs);
    assert(!V_im.allocated);
    n = V_im.size();
    std::swap(vtxs, next_vtxs);
    compact_t.stop();  // compact_t.reportTotal("compact time");
    cout << "      " << n << " vertices remain." << endl;

    // 6. relabel the edges with the new roots.
    timer relab_t;
    relab_t.start();
    parallel_for_bc(i, 0, m, (m > pbbs::kSequentialForThreshold), {
      size_t e_id = edge_ids[i];
      edge& e = edges[e_id];
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
    relab_t.stop();  // relab_t.reportTotal("relabel time");

    // 7. filter (or ignore) self-edges.
    auto self_loop_im = make_in_imap<bool>(
        n, [&](size_t i) { return !(edge_ids[i] & TOP_BIT); });
    auto edge_ids_im = make_array_imap(edge_ids, m);
    if (round == 0) {
      auto A = pbbs::pack(edge_ids_im, self_loop_im);
      m = A.size();
      next_edge_ids = A.get_array();
    } else {
      auto A =
          pbbs::pack(edge_ids_im, self_loop_im, pbbs::no_flag, next_edge_ids);
      m = A.size();
      assert(!A.allocated);
    }
    cout << "filter, m is now " << m << " n is now " << n << endl;
    std::swap(edge_ids, next_edge_ids);
    round++;
  }

  cout << "Boruvka finished: total edges added to MST = " << n_in_mst << endl;
  auto mst_im = make_array_imap(mst, n_in_mst);
  mst_im.allocated = true;
  return mst_im;
}

constexpr size_t sample_size = 2000;
inline size_t hash_to_range(size_t hsh, size_t range) { return hsh & range; }

inline size_t key_for_pair(uint32_t k1, uintE k2, pbbs::random rnd) {
  size_t key = (static_cast<size_t>(k1) << 32) + static_cast<size_t>(k2);
  return rnd.ith_rand(key);
}

template <template <class W> class vertex, class W>
inline edge_array<W> get_all_edges(graph<vertex<W>>& G) {
  auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return true;
  };
  return filter_all_edges(G, pred);
}

template <template <class W> class vertex, class W>
inline edge_array<W> get_top_k(graph<vertex<W>>& G, size_t k, pbbs::random r,
                               bool first_round = false) {
  if (k == G.m) {
    return get_all_edges(G);
  }
  using edge = std::tuple<uintE, uintE, W>;

  size_t n = G.n;
  size_t m = G.m;

  auto vertex_offs = array_imap<long>(G.n);
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                  { vertex_offs[i] = G.V[i].getOutDegree(); });
  pbbs::scan_add(vertex_offs, vertex_offs, pbbs::fl_scan_inclusive);

  auto sample_edges = array_imap<edge>(sample_size);
  auto lte = [&](const size_t& l, const size_t& r) { return l <= r; };

  parallel_for_bc(
      i, 0, sample_size, (sample_size > pbbs::kSequentialForThreshold), {
        size_t sample_edge = r.ith_rand(i) % m;
        uintE vtx = pbbs::binary_search(vertex_offs, sample_edge, lte);
        size_t ith = vertex_offs[vtx] - sample_edge - 1;
        uintE ngh;
        W wgh;
        std::tie(ngh, wgh) = G.V[vtx].get_ith_out_neighbor(vtx, ith);
        sample_edges[i] = std::make_tuple(vtx, ngh, wgh);
      });

  auto cmp_by_wgh = [](const edge& l, const edge& r) {
    return std::get<2>(l) < std::get<2>(r);
  };
  pbbs::sample_sort(sample_edges.start(), sample_edges.size(), cmp_by_wgh);

  // 2. find approximate splitter.
  size_t ind = ((double)(k * sample_edges.size())) / G.m;
  auto splitter = sample_edges[ind];
  W split_weight = std::get<2>(splitter);

  size_t first_ind = 0;
  size_t last_ind = 0;
  size_t ssize = sample_edges.size();
  parallel_for_bc(i, 0, ssize, (ssize > pbbs::kSequentialForThreshold), {
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
  cout << "split wgh is: " << split_weight << endl;
  cout << "fraction of sample composed by split_wgh = " << split_wgh_fraction
       << endl;

  // 3. filter edges based on splitter
  if (split_wgh_fraction < 0.2) {
    auto filter_pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
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
    cout << "frac of split_weight to take: " << frac_to_take << endl;
    size_t range = (1L << pbbs::log2_up(G.m)) - 1;
    size_t threshold = frac_to_take * range;
    // account for filtering directed edges
    threshold *= (first_round ? 2.0 : 1.0);
    auto filter_split_wgh_pred = [&](const uintE& src, const uintE& ngh,
                                     const W& wgh) {
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
          typename std::enable_if<!std::is_same<W, pbbs::empty>::value,
                                  int>::type = 0>
inline void MST(graph<vertex<W>>& GA, bool largemem = false) {
  using w_vertex = vertex<W>;
  using edge = std::tuple<uintE, uintE, W>;
  using ct = cas_type;

  size_t n = GA.n;
  cout << "n = " << n << endl;
  auto r = pbbs::default_random;

  auto exhausted = array_imap<bool>(n, [](size_t i) { return false; });
  auto parents = array_imap<uintE>(n, [](size_t i) { return i; });
  auto mst_edges = dyn_arr<edge>(n);

  auto min_edges = newA(ct, n);

  size_t n_active = n;
  uintE* vtxs = newA(uintE, n_active);
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                  { vtxs[i] = i; });
  uintE* next_vtxs = newA(uintE, n_active);

  size_t round = 0;
  while (GA.m > 0) {
    timer round_t;
    round_t.start();
    cout << endl;
    cout << "round = " << round << " n_active = " << n_active
         << " GA.m = " << GA.m << " MST size = " << mst_edges.size << endl;

    // find a prefix of lowest-weight edges.
    size_t split_idx = min((3 * n) / 2, (size_t)GA.m);
    if (largemem) {
      split_idx = (round == 0) ? min((9 * n) / 8, (size_t)GA.m)
                               : min(n / 2, (size_t)GA.m);
    }

    timer get_t;
    get_t.start();
    auto E = (round < n_filter_steps) ? get_top_k(GA, split_idx, r, round == 0)
                                      : get_all_edges(GA);
    get_t.stop();
    get_t.reportTotal("get time");
    size_t n_edges = E.non_zeros;
    cout << "Prefix size = " << split_idx << " #edges = " << n_edges
         << " G.m is now = " << GA.m << endl;

    // relabel edges
    auto edges = E.E;
    parallel_for_bc(i, 0, n_edges, (n_edges > pbbs::kSequentialForThreshold), {
      edge& e = edges[i];
      uintE u = std::get<0>(e);
      uintE v = std::get<1>(e);
      std::get<0>(e) = parents[u];
      std::get<1>(e) = parents[v];
    });

    // run Boruvka on the prefix and add new edges to mst_edges
    timer bt;
    bt.start();
    auto edge_ids =
        Boruvka(E, vtxs, next_vtxs, min_edges, parents, exhausted, n_active);
    bt.stop();
    bt.reportTotal("boruvka time");
    mst_edges.copyInF([&](size_t i) { return E.E[edge_ids[i]]; },
                      edge_ids.size());

    // reactivate vertices and reset exhausted
    auto vtx_im = make_array_imap(vtxs, n);
    timer pack_t;
    pack_t.start();
    n_active += ligra_utils::seq::packIndex(vtxs + n_active, exhausted.start(),
                                            (uintE)n);
    pack_t.stop();  // pack_t.reportTotal("reactivation pack");

    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      if (exhausted[i]) exhausted[i] = false;
    });

    // pointer jump: vertices that were made inactive could have had their
    // parents change.
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      size_t ctr = 0;
      while (parents[i] != parents[parents[i]]) {
        parents[i] = parents[parents[i]];
        ctr++;
      }
    });

    // pack out all edges in the graph that are shortcut by the added edges
    auto filter_pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      auto c_src = parents[src];
      auto c_ngh = parents[ngh];
      return c_src == c_ngh;
    };
    cout << "Filtering G, m = " << GA.m << endl;
    timer filter_t;
    filter_t.start();
    filter_edges(GA, filter_pred);
    filter_t.stop();
    filter_t.reportTotal("filter time");
    cout << "After filter, m is now " << GA.m << endl;
    round++;
    round_t.stop();
    round_t.reportTotal("round time");

    r = r.next();
  }
  cout << "#edges in output mst: " << mst_edges.size << endl;
  auto wgh_imap = make_in_imap<size_t>(
      mst_edges.size, [&](size_t i) { return std::get<2>(mst_edges.A[i]); });
  cout << "total weight = " << pbbs::reduce_add(wgh_imap) << endl;

  mst_edges.del();
  free(min_edges);
}

template <
    template <class W> class vertex, class W,
    typename std::enable_if<std::is_same<W, pbbs::empty>::value, int>::type = 0>
inline uint32_t* MST(graph<vertex<W>>& GA, bool largeem = false) {
  cout << "Unimplemented for unweighted graphs" << endl;
  exit(0);
}
}  // namespace MST_boruvka

namespace MST_spec_for {
constexpr size_t sample_size = 10000;
inline size_t hash_to_range(size_t hsh, size_t range) { return hsh & range; }

inline size_t key_for_pair(uint32_t k1, uintE k2, pbbs::random rnd) {
  size_t key = (static_cast<size_t>(k1) << 32) + static_cast<size_t>(k2);
  return rnd.ith_rand(key);
}

template <template <class W> class vertex, class W, class UF>
inline edge_array<W> get_remaining(graph<vertex<W>>& G, size_t k, UF& uf,
                                   pbbs::random r) {
  auto filter_pred = [&](const uint32_t& src, const uintE& ngh, const W& wgh) {
    if (src < ngh) {
      return 2;  // return in array
    } else {
      return 1;  // filter, don't return in array
    }
    return 0;
  };
  return filter_edges(G, filter_pred);
}

template <template <class W> class vertex, class W, class UF>
inline void pack_shortcut_edges(graph<vertex<W>>& G, UF& uf) {
  auto filter_pred = [&](const uint32_t& src, const uintE& ngh,
                         const W& wgh) -> bool {
    if (src > ngh) {
      return true;
    }
    auto c_src = uf.find(src);
    auto c_ngh = uf.find(ngh);
    return c_src == c_ngh;
  };
  cout << "Calling filter, G.m = " << G.m << endl;
  filter_edges(G, filter_pred);
  cout << "G.m is now " << G.m << endl;
}

template <template <class W> class vertex, class W, class UF>
inline edge_array<W> get_top_k(graph<vertex<W>>& G, size_t k, UF& uf,
                               pbbs::random r, bool first_round = false) {
  if (k == G.m) {
    return get_remaining(G, k, uf, r);
  }

  size_t m = (first_round) ? (G.m / 2) : G.m;
  size_t range = (1L << pbbs::log2_up(G.m)) - 1;
  size_t scaled_size = sample_size * (((double)range) / m);

  // 1. Sample sample_size many edges and sort them.
  auto pred = [&](const uint32_t& src, const uintE& ngh, const W& wgh) -> bool {
    if (src < ngh) {
      size_t hash_val = hash_to_range(key_for_pair(src, ngh, r), range);
      return hash_val < scaled_size;
    }
    return 0;
  };
  auto sampled_e = sample_edges(G, pred);
  if (sampled_e.non_zeros == 0) {
    cout << "non_zeros = 0" << endl;
    exit(0);
    return get_remaining(G, k, uf, r);
  }
  auto cmp_by_wgh = [](const std::tuple<uint32_t, uintE, intE>& l,
                       const std::tuple<uintE, uintE, intE>& r) {
    return std::get<2>(l) < std::get<2>(r);
  };
  pbbs::sample_sort(sampled_e.E, sampled_e.non_zeros, cmp_by_wgh);

  // 2. Get approximate splitter.
  size_t ind = ((double)(k * sampled_e.non_zeros)) / G.m;
  auto splitter = sampled_e.E[ind];
  int32_t split_weight = std::get<2>(splitter);
  sampled_e.del();
  cout << "split wgh is: " << split_weight << endl;

  // 3. Filter edges based on splitter
  auto filter_pred = [&](const uint32_t& src, const uintE& ngh, const W& wgh) {
    if (wgh <= split_weight) {
      if (src < ngh) {
        return 2;  // return in array
      } else {
        return 1;  // filter, but don't return in array
      }
    }
    return 0;
  };
  return filter_edges(G, filter_pred);
}

template <template <class W> class vertex, class W,
          typename std::enable_if<!std::is_same<W, pbbs::empty>::value,
                                  int>::type = 0>
inline void MST(graph<vertex<W>>& GA) {
  using w_vertex = vertex<W>;
  using res = reservation<uintE>;
  using edge_t = std::tuple<uintE, uintE, W>;

  size_t n = GA.n;
  auto r = pbbs::default_random;
  auto uf = UnionFind(n);

  auto mst_edges = dyn_arr<edge_t>(n);

  size_t iter = 0;
  while (GA.m > 0) {
    cout << "iter = " << iter << " m = " << GA.m << endl;
    // 1. get weight of k'th smallest edge, and all edges smaller.
    size_t split_idx = min(n, (size_t)GA.m);
    timer get_t;
    get_t.start();
    auto edges = get_top_k(GA, split_idx, uf, r, (iter == 0));
    get_t.stop();
    get_t.reportTotal("get time");
    size_t n_edges = edges.non_zeros;
    auto cmp_by_wgh = [](const std::tuple<uint32_t, uintE, W>& l,
                         const std::tuple<uintE, uintE, W>& r) {
      return std::get<2>(l) < std::get<2>(r);
    };
    pbbs::sample_sort(edges.E, n_edges, cmp_by_wgh);
    cout << "Prefix size = " << split_idx << " #edges = " << n_edges
         << " G.m is now = " << GA.m << endl;

    // 2. initialize reservations, copy edge info, and run UF step.
    auto R = newA(res, n);
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold),
                    { R[i] = res(); });
    array_imap<bool> mstFlags =
        array_imap<bool>(n_edges, [](size_t i) { return 0; });

    auto UFStep = make_uf_step<uintE>(edges, R, mstFlags, uf);
    speculative_for<uintE>(UFStep, 0, n_edges, 8);

    UFStep.del();
    free(R);
    auto edge_im =
        make_in_imap<edge_t>(n_edges, [&](size_t i) { return edges.E[i]; });
    auto edges_ret = pbbs::pack(edge_im, mstFlags);
    cout << "added " << edges_ret.size() << endl;
    mst_edges.copyIn(edges_ret, edges_ret.size());
    edges.del();
    mstFlags.del();

    timer pack_t;
    pack_t.start();
    pack_shortcut_edges(GA, uf);
    pack_t.stop();
    pack_t.reportTotal("pack time");
    iter++;
  }
  cout << "n in mst: " << mst_edges.size << endl;
  auto wgh_imap = make_in_imap<size_t>(
      mst_edges.size, [&](size_t i) { return std::get<2>(mst_edges.A[i]); });
  cout << "wgh = " << pbbs::reduce_add(wgh_imap) << endl;

  mst_edges.del();
}

template <
    template <class W> class vertex, class W,
    typename std::enable_if<std::is_same<W, pbbs::empty>::value, int>::type = 0>
inline uint32_t* MST(graph<vertex<W>>& GA) {
  exit(0);
}
}  // namespace MST_spec_for
