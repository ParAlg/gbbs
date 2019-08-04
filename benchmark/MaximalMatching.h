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

#include "bridge.h"
#include "ligra.h"
#include "speculative_for.h"

#include "pbbslib/dyn_arr.h"
#include "pbbslib/sparse_table.h"

namespace mm {

  constexpr size_t n_filter_steps = 5;

inline size_t hash_to_range(size_t hsh, size_t range) { return hsh & range; }
inline size_t key_for_pair(uintE k1, uintE k2, pbbslib::random rnd) {
  size_t l = std::min(k1, k2);
  size_t r = std::max(k1, k2);
  size_t key = (l << 32) + r;
  return rnd.ith_rand(key);
}

  struct matchStep {
    using edge = std::tuple<uintE, uintE, pbbs::empty>;
    edge* E;
    uintE* R;
    uintE* M;
    bool* matched;
    matchStep(edge* _E, uintE* _R, uintE* M, bool* m) : E(_E), R(_R), M(M), matched(m) {}

    bool reserve(uintE i) {
      uintE u = std::get<0>(E[i]);
      uintE v = std::get<1>(E[i]);
      if (matched[u] || matched[v] || (u == v)) return 0;
      reserveLoc<uintE>(&R[u], i);
      reserveLoc<uintE>(&R[v], i);
      return 1;
    }

    bool commit(uintE i) {
      uintE u = std::get<0>(E[i]);
      uintE v = std::get<1>(E[i]);
      if (R[v] == i) {
        R[v] = UINT_E_MAX; // reset
        if (R[u] == i) {
          matched[u] = matched[v] = 1;
          // mark edge in both directions
          M[u] = v;
          M[v] = u;
          return 1;
        }
      } else if (R[u] == i)
        R[u] = UINT_E_MAX; // reset
      return 0;
    }
  };

  template <class G>
  inline edge_array<pbbs::empty> get_all_edges(G& GA, bool* matched,
                                 pbbslib::random r) {
    cout << "fetching all edges!" << endl;
    using W = typename G::weight_type;
    using edge = std::tuple<uintE, uintE, pbbs::empty>;
    size_t n = GA.n;

    auto allvtxs = pbbslib::make_sequence<uintE>(n, [&] (size_t i) { return i; });
    auto is_live = [&] (uintE v) { return !matched[v]; };
    auto live_vtxs = pbbs::filter(allvtxs, is_live);

    assert(live_vtxs.size() > 0);

    auto vtx_degs = sequence<size_t>(live_vtxs.size());
    parallel_for(0, live_vtxs.size(), [&] (size_t i) {
      vtx_degs[i] = GA.getOutDegree(live_vtxs[i]);
    }, 1024);

    size_t m = pbbslib::scan_add_inplace(vtx_degs.slice());

    auto eout = sequence<edge>(m);
    parallel_for(0, live_vtxs.size(), [&] (size_t i) {
      uintE vtx_id = live_vtxs[i];
      auto vtx = GA.get_vertex(vtx_id);
      size_t off = vtx_degs[i];
      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        return true;
      };
      auto write_f = [&] (const uintE& ngh, const size_t& offset, const bool& val) {
        eout[offset] = std::make_tuple(vtx_id, ngh, pbbs::empty());
      };
      vtx.copyOutNgh(vtx_id, off, map_f, write_f);
    }, 1);

    edge_array<pbbs::empty> e_arr;
    e_arr.num_rows = GA.n;
    e_arr.num_cols = GA.n;
    e_arr.non_zeros = eout.size();
    e_arr.E = eout.to_array();
    return e_arr;
  }

  template <class G>
  inline edge_array<pbbs::empty> get_edges(G& GA, size_t k, bool* matched,
                                 pbbslib::random r) {
    using W = typename G::weight_type;
    using edge = std::tuple<uintE, uintE, pbbs::empty>;
    size_t n = GA.n;

    size_t edges_remaining = GA.m;

    using K = std::tuple<uintE, uintE>;
    using V = pbbs::empty;
    auto key_hash = [&] (const K& k) -> size_t {
      return key_for_pair(std::get<0>(k), std::get<1>(k), pbbs::random());
    };
    auto empty = std::make_tuple(std::make_tuple(UINT_E_MAX, UINT_E_MAX), pbbs::empty());
    auto tab = make_sparse_table(2*k, empty, key_hash);

    size_t range = (1L << pbbslib::log2_up(edges_remaining)) - 1;
    // edges_remaining edges going into [0, pbbslib::log2_up(edges_remaining))
    // we want ~k edges, meaning we scale k by the same factor

    double scale = double(range) / edges_remaining; // 2 >= scale > 1
    assert(k > 1);
    size_t scaled_k = (k * scale);

    double frac_to_take = double(scaled_k) / double(range);
    size_t threshold = frac_to_take * range;

    auto map_f = [&](const uintE& src, const uintE& ngh, const W& w) {
      if (hash_to_range(key_for_pair(src, ngh, r), range) < threshold) {
        auto kv = std::make_tuple(std::make_tuple(src, ngh), pbbs::empty());
        tab.insert(kv);
      }
    };
    GA.map_edges(map_f);

    auto ee = tab.entries();

    edge_array<pbbs::empty> e_arr;
    e_arr.num_rows = GA.n;
    e_arr.num_cols = GA.n;
    e_arr.non_zeros = ee.size();
    e_arr.E = (std::tuple<uintE, uintE, pbbs::empty>*)ee.to_array();
    return e_arr;
  }

};  // namespace mm

// Runs a constant number of filters on the whole graph, and runs the
// prefix-based algorithm on them. Finishes off the rest of the graph with the
// prefix-based algorithm.
template <class G>
inline sequence<std::tuple<uintE, uintE, pbbs::empty>> MaximalMatching(G& GA) {
  using W = typename G::weight_type;
  using edge = std::tuple<uintE, uintE, pbbs::empty>;

  timer ft; ft.start();
  auto filter_pred = [&] (const uintE& s, const uintE& d, const W& wgh) {
    return s < d;
  };
  auto PG = filter_graph(GA, filter_pred);
  ft.stop(); ft.reportTotal("filter time");

  timer mt;
  mt.start();
  size_t n = GA.n;
  auto r = pbbslib::random();

  auto R = sequence<uintE>(n, [&](size_t i) { return UINT_E_MAX; });
  auto matching = sequence<uintE>(n, [&](size_t i) { return UINT_E_MAX; });
  auto matched = sequence<bool>(n, [&](size_t i) { return false; }); // bitvector

  size_t k = GA.n; // ((3 * GA.n) / 2);

  size_t round = 0;
  timer gete;
  timer eff;
  while (PG.m > 0) {
    gete.start();
    timer get_t; get_t.start();
    auto e_arr = (PG.m <= k) ? mm::get_all_edges(PG, matched.begin(), r) : mm::get_edges(PG, k, matched.begin(), r);
    get_t.stop(); get_t.reportTotal("get time");

    auto eim_f = [&](size_t i) { return e_arr.E[i]; };
    auto eim = pbbslib::make_sequence<edge>(e_arr.non_zeros, eim_f);
    gete.stop();

    std::cout << "Got: " << e_arr.non_zeros << " edges "
              << " PG.m is now: " << PG.m << "\n";
    mm::matchStep mStep(e_arr.E, R.begin(), matching.begin(), matched.begin());
    eff.start();
    eff_for<size_t>(mStep, 0, e_arr.non_zeros, 50, 0, PG.n);
    eff.stop();
    e_arr.del();

    auto removed_seq = pbbs::delayed_seq<size_t>(n, [&] (size_t i) {
      if (matching[i] != UINT_E_MAX && PG.getOutDegree(i) > 0) {
        size_t deg =  PG.getOutDegree(i);
        PG.get_vertex(i).clear_vertex();
        return deg;
      }
      return static_cast<size_t>(0);
    });
    size_t deg_removed = pbbslib::reduce_add(removed_seq);
    PG.m -= deg_removed;
    std::cout << "removed: " << deg_removed << " many edges" << "\n";

    auto pack_pred = [&] (const uintE& u, const uintE& v, const W& wgh) {
      assert(!matched[u]);
      return !matched[v]; // only keep edges to unmatched v's
    };
    filter_graph(PG, pack_pred);

    // TODO: delete eid, other data returned this round
    round++;
    r = r.next();
  }

  auto matching_seq = pbbs::delayed_seq<std::tuple<uintE, uintE, pbbs::empty>>(n, [&] (size_t i) {
    return std::make_tuple(i, matching[i], pbbs::empty());
  });
  auto matching_filter = [&] (const std::tuple<uintE, uintE, pbbs::empty>& edge) {
    uintE u = std::get<0>(edge);
    uintE v = std::get<1>(edge);
    return (v != UINT_E_MAX) && (u < v);
  };

  auto ret = pbbs::filter(matching_seq, matching_filter);
  mt.stop();
  eff.reportTotal("eff for time");
  gete.reportTotal("get edges time");
  mt.reportTotal("Matching time");
  PG.del();
  std::cout << "matching size = " << ret.size() << "\n";
  return std::move(ret);
}

template <class G, class Seq>
inline void verify_matching(G& GA, Seq& matching) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  auto ok = sequence<bool>(n, [](size_t i) { return 1; });
  auto matched = sequence<uintE>(n, [](size_t i) { return 0; });

  // Check that this is a valid matching
  par_for(0, matching.size(), [&] (size_t i) {
                    const auto& edge = matching[i];
                    uintE u = std::get<0>(edge);
                    uintE v = std::get<1>(edge);
                    pbbslib::write_add(&matched[u], 1);
                    pbbslib::write_add(&matched[v], 1);
                  });

  bool valid = true;
  par_for(0, n, [&] (size_t i) {
    if (matched[i] > 1) valid = false;
  });
  assert(valid == true);

  // Check maximality of the matching
  auto map2_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    if (!matched[src] && !matched[ngh]) {
      // could have added this edge, increasing the size of the matching
      ok[src] = 0;
      ok[ngh] = 0;
    }
  };
  par_for(0, n, 1, [&] (size_t i) { GA.get_vertex(i).mapOutNgh(i, map2_f); });

  auto ok_f = [&](size_t i) { return ok[i]; };
  auto ok_im = pbbslib::make_sequence<size_t>(n, ok_f);
  size_t n_ok = pbbslib::reduce_add(ok_im);
  if (n == n_ok) {
    std::cout << "Matching OK! matching size is: " << matching.size() << "\n";
  } else {
    std::cout << "Matching invalid---" << (n - n_ok)
              << " vertices saw bad neighborhoods."
              << "\n";
  }
}
