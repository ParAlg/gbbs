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
#include "ligra.h"
#include "speculative_for.h"
#include "union_find.h"

#include "pbbslib/binary_search.h"
#include "pbbslib/dyn_arr.h"
#include "pbbslib/random.h"
#include "pbbslib/sample_sort.h"

namespace MST_boruvka {

constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;

// the number of filtering steps to run
constexpr size_t n_filter_steps = 5;

//struct alignas(8) cas_type {
//  uintE index;
//  int32_t weight;
//  cas_type() {
//    index = UINT_E_MAX;
//    weight = std::numeric_limits<int32_t>::max();
//  }
//  cas_type(uintE _index, int32_t _w) : index(_index), weight(_w) {}
//};


using cas_type = std::pair<uintE, int32_t>;

template <class W, class M, class P, class D>
inline sequence<uintE> Boruvka(edge_array<W>& E, uintE*& vtxs,
                                 uintE*& next_vtxs, M& min_edges, P& parents,
                                 D& exhausted, size_t& n) {
  using ct = cas_type;
  using edge = std::tuple<uintE, uintE, W>;
  size_t m = E.non_zeros;
  auto edges = E.E;
  auto less = [](const ct& a, const ct& b) {
    // returns true if (weight is <) or (weight = and index is <)
    return (a.second < b.second) || (a.second == b.second && a.first < b.first);
  };

  uintE* edge_ids = pbbslib::new_array_no_init<uintE>(m);
  par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) { edge_ids[i] = i; });
  uintE* next_edge_ids = nullptr;

  auto new_mst_edges = sequence<uintE>(n, UINT_E_MAX);
  auto is_root = sequence<bool>(n);

  // Stores edge indices that join the MST.
  uintE* mst = pbbslib::new_array_no_init<uintE>(n);
  size_t n_in_mst = 0;
  size_t round = 0;

  while (n > 1 && m > 0) {
    std::cout << "Boruvka round: " << round << " n: " << n << " m: " << m
              << "\n";

    timer init_t;
    init_t.start();
    par_for(0, n, [&] (size_t i) {
      uintE v = vtxs[i];
      // stores an (index, weight) pair
      min_edges[v] = std::make_pair(UINT_E_MAX, std::numeric_limits<int32_t>::max());
    });
    init_t.stop();  // init_t.reportTotal("init time");

    // 1. write_min to select the minimum edge out of each component.
    timer min_t;
    min_t.start();
    par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      uintE e_id = edge_ids[i];
      const edge& e = edges[e_id];
      ct cas_e(e_id, std::get<2>(e));
      pbbslib::write_min(min_edges + std::get<0>(e), cas_e, less);
      pbbslib::write_min(min_edges + std::get<1>(e), cas_e, less);
    });
    min_t.stop();  // min_t.reportTotal("write min time");

    // 2. test whether vertices found an edge incident to them
    timer mark_t;
    mark_t.start();
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      uintE v = vtxs[i];
      const auto& e = min_edges[v];
      if (e.first == UINT_E_MAX) {
        // no more edges incident to v in this batch.
        exhausted[v] = true;
        is_root[i] = false;
        new_mst_edges[i] = UINT_E_MAX;
        debug(uintE pv = parents[v];
        assert(pv == v););
      } else {
        uintE ind = e.first;
        const auto& edge = edges[ind];
        uintE u = std::get<0>(edge) ^ std::get<1>(edge) ^ v;
        // pick the higher endpoint as the root.
        if (u < v && ind == min_edges[u].first) {
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
    n_in_mst += pbbslib::filterf(new_mst_edges.begin(), mst + n_in_mst, n,
                              [](uintE v) { return v != UINT_E_MAX; });
    std::cout << "      " << n_in_mst << " edges added to mst."
              << "\n";
    filter_t.stop();  // filter_t.reportTotal("filter time");

    // 4. pointer jump to find component centers.
    timer jump_t;
    jump_t.start();
    par_for(0, n, [&] (size_t i) {
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
    auto vtxs_im = pbbslib::make_sequence<uintE>(vtxs, n);

    n = pbbslib::pack_out(vtxs_im, is_root, pbbslib::make_sequence(next_vtxs, m));
    std::swap(vtxs, next_vtxs);
    compact_t.stop();  // compact_t.reportTotal("compact time");
    std::cout << "      " << n << " vertices remain."
              << "\n";

    // 6. relabel the edges with the new roots.
    timer relab_t;
    relab_t.start();
    par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
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
    auto self_loop_f = [&](size_t i) { return !(edge_ids[i] & TOP_BIT); };
    auto self_loop_im = pbbslib::make_sequence<bool>(n, self_loop_f);
    auto edge_ids_im = pbbslib::make_sequence<uintE>(edge_ids, m);
    if (round == 0) {
      auto A = pbbslib::pack(edge_ids_im, self_loop_im);
      m = A.size();
      next_edge_ids = A.to_array();
    } else {
      m = pbbslib::pack_out(edge_ids_im, self_loop_im, pbbslib::make_sequence(next_edge_ids, m));
    }
    std::cout << "filter, m is now " << m << " n is now " << n << "\n";
    std::swap(edge_ids, next_edge_ids);
    round++;
  }

  // TODO check about freeing next_edge_ids and edge_ids
  std::cout << "Boruvka finished: total edges added to MST = " << n_in_mst
            << "\n";
  pbbslib::free_array(edge_ids);
  pbbslib::free_array(next_edge_ids);
  auto mst_im = sequence<uintE>(mst, n_in_mst); // allocated
  return mst_im;
}

constexpr size_t sample_size = 2000;
inline size_t hash_to_range(size_t hsh, size_t range) { return hsh & range; }

inline size_t key_for_pair(uint32_t k1, uintE k2, pbbslib::random rnd) {
  size_t key = (static_cast<size_t>(k1) << 32) + static_cast<size_t>(k2);
  return rnd.ith_rand(key);
}

template <class G>
inline edge_array<W> get_all_edges(G& G) {
  using W = typename G::weight_type;
  auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return true;
  };
  return filter_all_edges(G, pred);
}

template <class G>
inline edge_array<W> get_top_k(G& G, size_t k, pbbslib::random r,
                               bool first_round = false) {
  using W = typename G::weight_type;
  if (k == static_cast<size_t>(G.m)) {
    return get_all_edges(G);
  }
  using edge = std::tuple<uintE, uintE, W>;
  timer st; st.start();

  size_t n = G.n;
  size_t m = G.m;

  auto vertex_offs = sequence<long>(G.n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { vertex_offs[i] = G.V[i].getOutDegree(); });
  pbbslib::scan_add_inplace(vertex_offs.slice(), pbbslib::fl_scan_inclusive);

  auto sample_edges = sequence<edge>(sample_size);
  auto lte = [&](const size_t& l, const size_t& r) { return l <= r; };

  par_for(0, sample_size, pbbslib::kSequentialForThreshold, [&] (size_t i) {
        size_t sample_edge = r.ith_rand(i) % m;
        uintE vtx = pbbslib::binary_search(vertex_offs, sample_edge, lte);
        size_t ith = vertex_offs[vtx] - sample_edge - 1;
        uintE ngh;
        W wgh;
        std::tie(ngh, wgh) = G.V[vtx].get_ith_out_neighbor(vtx, ith);
        sample_edges[i] = std::make_tuple(vtx, ngh, wgh);
      });

  auto cmp_by_wgh = [](const edge& l, const edge& r) {
    return std::get<2>(l) < std::get<2>(r);
  };
  pbbslib::sample_sort_inplace(sample_edges.slice(), cmp_by_wgh);

  // 2. find approximate splitter.
  size_t ind = ((double)(k * sample_edges.size())) / G.m;
  auto splitter = sample_edges[ind];
  W split_weight = std::get<2>(splitter);

  size_t first_ind = 0;
  size_t last_ind = 0;
  size_t ssize = sample_edges.size();
  par_for(0, ssize, pbbslib::kSequentialForThreshold, [&] (size_t i) {
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
  std::cout << "split wgh is: " << split_weight << "\n";
  std::cout << "fraction of sample composed by split_wgh = "
            << split_wgh_fraction << "\n";
  st.stop(); st.reportTotal("startup time");

  // 3. filter edges based on splitter
  if (split_wgh_fraction < 0.2) {
    auto filter_pred = [&](const uintE& src, const uintE& ngh, const W& wgh) -> int {
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
    size_t range = (1L << pbbslib::log2_up(G.m)) - 1;
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

template <class G,
          typename std::enable_if<!std::is_same<typename G::weight_type, pbbslib::empty>::value,
                                  int>::type = 0>
inline void MST(G& GA, bool largemem = false) {
  using W = typename G::weight_type;
  using edge = std::tuple<uintE, uintE, W>;
  using ct = cas_type;

  size_t n = GA.n;
  std::cout << "n = " << n << "\n";
  auto r = pbbslib::random();

  auto exhausted = sequence<bool>(n, [](size_t i) { return false; });
  auto parents = sequence<uintE>(n, [](size_t i) { return i; });
  auto mst_edges = pbbslib::dyn_arr<edge>(n);

  auto min_edges = pbbslib::new_array_no_init<ct>(n);

  size_t n_active = n;
  uintE* vtxs = pbbslib::new_array_no_init<uintE>(n_active);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { vtxs[i] = i; });
  uintE* next_vtxs = pbbslib::new_array_no_init<uintE>(n_active);

  size_t round = 0;
  while (GA.m > 0) {
    timer round_t;
    round_t.start();
    std::cout << "\n";
    std::cout << "round = " << round << " n_active = " << n_active
              << " GA.m = " << GA.m << " MST size = " << mst_edges.size << "\n";

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
    get_t.reportTotal("get time");
    size_t n_edges = E.non_zeros;
    std::cout << "Prefix size = " << split_idx << " #edges = " << n_edges
              << " G.m is now = " << GA.m << "\n";

    // relabel edges
    auto edges = E.E;
    if (round > 0) {
      par_for(0, n_edges, pbbslib::kSequentialForThreshold, [&] (size_t i) {
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
    auto edge_ids =
        Boruvka(E, vtxs, next_vtxs, min_edges, parents, exhausted, n_active);
    bt.stop();
    bt.reportTotal("boruvka time");
    mst_edges.copyInF([&](size_t i) { return E.E[edge_ids[i]]; },
                      edge_ids.size());

    // reactivate vertices and reset exhausted
    timer pack_t;
    pack_t.start();
    // TODO: remove dependency on ligra_utils


    auto vtx_range = pbbs::make_range(vtxs+n_active, vtxs+n);
    n_active += pbbslib::pack_index_out(exhausted.slice(), vtx_range);
//    n_active += ligra_utils::seq::packIndex(vtxs + n_active, exhausted.begin(),
//                                            (uintE)n);
    pack_t.stop();  // pack_t.reportTotal("reactivation pack");

    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      if (exhausted[i]) exhausted[i] = false;
    });


    // pointer jump: vertices that were made inactive could have had their
    // parents change.
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      size_t ctr = 0;
      while (parents[i] != parents[parents[i]]) {
        parents[i] = parents[parents[i]];
        ctr++;
      }
    });

    // pack out all edges in the graph that are shortcut by the added edges
    auto filter_pred = [&](const uintE& src, const uintE& ngh, const W& wgh) -> int {
      auto c_src = parents[src];
      auto c_ngh = parents[ngh];
      return c_src == c_ngh;
    };
    std::cout << "Filtering G, m = " << GA.m << "\n";
    timer filter_t;
    filter_t.start();
    filter_edges(GA, filter_pred);
    filter_t.stop();
    filter_t.reportTotal("filter time");
    std::cout << "After filter, m is now " << GA.m << "\n";
    round++;
    round_t.stop();
    round_t.reportTotal("round time");

    r = r.next();
  }
  std::cout << "#edges in output mst: " << mst_edges.size << "\n";
  auto wgh_imap_f = [&](size_t i) { return std::get<2>(mst_edges.A[i]); };
  auto wgh_imap = pbbslib::make_sequence<size_t>(
      mst_edges.size, wgh_imap_f);
  std::cout << "total weight = " << pbbslib::reduce_add(wgh_imap) << "\n";

  mst_edges.clear();
  pbbslib::free_array(min_edges);
}

template <class G,
    typename std::enable_if<std::is_same<typename G::weight_type, pbbslib::empty>::value, int>::type = 0>
inline uint32_t* MST(G& GA, bool largeem = false) {
  std::cout << "Unimplemented for unweighted graphs"
            << "\n";
  exit(0);
}
}  // namespace MST_boruvka

//namespace MST_spec_for {
//constexpr size_t sample_size = 10000;
//inline size_t hash_to_range(size_t hsh, size_t range) { return hsh & range; }
//
//inline size_t key_for_pair(uint32_t k1, uintE k2, pbbslib::random rnd) {
//  size_t key = (static_cast<size_t>(k1) << 32) + static_cast<size_t>(k2);
//  return rnd.ith_rand(key);
//}
//
//template <class G, class UF>
//inline edge_array<W> get_remaining(G& G, size_t k, UF& uf,
//                                   pbbslib::random r) {
//  using W = typename G::weight_type;
//  auto filter_pred = [&](const uint32_t& src, const uintE& ngh, const W& wgh) {
//    if (src < ngh) {
//      return 2;  // return in array
//    } else {
//      return 1;  // filter, don't return in array
//    }
//    return 0;
//  };
//  return filter_edges(G, filter_pred);
//}
//
//template <template <class W> class vertex, class W, class UF>
//inline void pack_shortcut_edges(G& G, UF& uf) {
//  using W = typename G::weight_type;
//  auto filter_pred = [&](const uint32_t& src, const uintE& ngh,
//                         const W& wgh) -> int {
//    if (src > ngh) {
//      return true;
//    }
//    auto c_src = uf.find(src);
//    auto c_ngh = uf.find(ngh);
//    return c_src == c_ngh;
//  };
//  std::cout << "Calling filter, G.m = " << G.m << "\n";
//  filter_edges(G, filter_pred, no_output);
//  std::cout << "G.m is now " << G.m << "\n";
//}
//
//template <class G, class UF>
//inline edge_array<W> get_top_k(G& G, size_t k, UF& uf,
//                               pbbslib::random r, bool first_round = false) {
//  using W = typename G::weight_type;
//  if (k == static_cast<size_t>(G.m)) {
//    return get_remaining(G, k, uf, r);
//  }
//
//  size_t m = (first_round) ? (G.m / 2) : G.m;
//  size_t range = (1L << pbbslib::log2_up(G.m)) - 1;
//  size_t scaled_size = sample_size * (((double)range) / m);
//
//  // 1. Sample sample_size many edges and sort them.
//  auto pred = [&](const uint32_t& src, const uintE& ngh, const W& wgh) -> bool {
//    if (src < ngh) {
//      size_t hash_val = hash_to_range(key_for_pair(src, ngh, r), range);
//      return hash_val < scaled_size;
//    }
//    return 0;
//  };
//  auto sampled_e = sample_edges(G, pred);
//  if (sampled_e.non_zeros == 0) {
//    std::cout << "non_zeros = 0"
//              << "\n";
//    exit(0);
//    return get_remaining(G, k, uf, r);
//  }
//  auto cmp_by_wgh = [](const std::tuple<uint32_t, uintE, intE>& l,
//                       const std::tuple<uintE, uintE, intE>& r) {
//    return std::get<2>(l) < std::get<2>(r);
//  };
//  pbbslib::sample_sort(pbbslib::make_sequence(sampled_e.E, sampled_e.non_zeros), cmp_by_wgh);
//
//  // 2. Get approximate splitter.
//  size_t ind = ((double)(k * sampled_e.non_zeros)) / G.m;
//  auto splitter = sampled_e.E[ind];
//  int32_t split_weight = std::get<2>(splitter);
//  sampled_e.del();
//  std::cout << "split wgh is: " << split_weight << "\n";
//
//  // 3. Filter edges based on splitter
//  auto filter_pred = [&](const uint32_t& src, const uintE& ngh, const W& wgh) {
//    if (wgh <= split_weight) {
//      if (src < ngh) {
//        return 2;  // return in array
//      } else {
//        return 1;  // filter, but don't return in array
//      }
//    }
//    return 0;
//  };
//  return filter_edges(G, filter_pred);
//}
//
//template <class G
//          typename std::enable_if<!std::is_same<W, pbbslib::empty>::value,
//                                  int>::type = 0>
//inline void MST(G& GA) {
//  using W = typename G::weight_type;
//  using res = reservation<uintE>;
//  using edge_t = std::tuple<uintE, uintE, W>;
//
//  size_t n = GA.n;
//  auto r = pbbslib::random();
//  auto uf = UnionFind(n);
//
//  auto mst_edges = pbbslib::dyn_arr<edge_t>(n);
//
//  size_t iter = 0;
//  while (GA.m > 0) {
//    std::cout << "iter = " << iter << " m = " << GA.m << "\n";
//    // 1. get weight of k'th smallest edge, and all edges smaller.
//    size_t split_idx = std::min(n, (size_t)GA.m);
//    timer get_t;
//    get_t.start();
//    auto edges = get_top_k(GA, split_idx, uf, r, (iter == 0));
//    get_t.stop();
//    get_t.reportTotal("get time");
//    size_t n_edges = edges.non_zeros;
//    auto cmp_by_wgh = [](const std::tuple<uint32_t, uintE, W>& l,
//                         const std::tuple<uintE, uintE, W>& r) {
//      return std::get<2>(l) < std::get<2>(r);
//    };
//    pbbslib::sample_sort(pbbslib::make_sequence(edges.E, n_edges), cmp_by_wgh);
//    std::cout << "Prefix size = " << split_idx << " #edges = " << n_edges
//              << " G.m is now = " << GA.m << "\n";
//
//    // 2. initialize reservations, copy edge info, and run UF step.
//    auto R = pbbslib::new_array_no_init<res>(n);
//    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
//                    { R[i] = res(); });
//    sequence<bool> mstFlags =
//        sequence<bool>(n_edges, [](size_t i) { return 0; });
//
//    auto UFStep = make_uf_step<uintE>(edges, R, mstFlags, uf);
//    speculative_for<uintE>(UFStep, 0, n_edges, 8);
//
//    UFStep.clear();
//    pbbslib::free_array(R);
//    auto edge_imap_f = [&](size_t i) { return edges.E[i]; };
//    auto edge_im =
//        pbbslib::make_sequence<edge_t>(n_edges, edge_imap_f);
//    auto edges_ret = pbbslib::pack(edge_im, mstFlags);
//    std::cout << "added " << edges_ret.size() << "\n";
//    mst_edges.copyIn(edges_ret, edges_ret.size());
//    edges.del();
//    mstFlags.clear();
//
//    timer pack_t;
//    pack_t.start();
//    pack_shortcut_edges(GA, uf);
//    pack_t.stop();
//    pack_t.reportTotal("pack time");
//    iter++;
//  }
//  std::cout << "n in mst: " << mst_edges.size << "\n";
//  auto wgh_imap_f = [&](size_t i) { return std::get<2>(mst_edges.A[i]); };
//  auto wgh_imap = pbbslib::make_sequence<size_t>(
//      mst_edges.size, wgh_imap_f);
//  std::cout << "wgh = " << pbbslib::reduce_add(wgh_imap) << "\n";
//
//  mst_edges.del();
//}
//
//template <class G
//    typename std::enable_if<std::is_same<typename G::weight_type, pbbslib::empty>::value, int>::type = 0>
//inline uint32_t* MST(G& GA) {
//  exit(0);
//}
//}  // namespace MST_spec_for
