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

#include "lib/dyn_arr.h"
#include "lib/random.h"
#include "lib/random_shuffle.h"
#include "lib/sample_sort.h"
#include "lib/speculative_for.h"
#include "ligra.h"

using namespace std;

namespace mm {
constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;
constexpr size_t n_filter_steps = 5;

template <class W>
struct matchStep {
  using edge = tuple<uintE, uintE, W>;
  edge* E;
  uintE* R;
  bool* matched;
  matchStep(edge* _E, uintE* _R, bool* m) : E(_E), R(_R), matched(m) {}

  bool reserve(uintE i) {
    uintE u = get<0>(E[i]);
    uintE v = get<1>(E[i]);
    if (matched[u] || matched[v] || (u == v)) return 0;
    reserveLoc<uintE>(R[u], i);
    reserveLoc<uintE>(R[v], i);
    return 1;
  }

  bool commit(uintE i) {
    uintE u = get<0>(E[i]);
    uintE v = get<1>(E[i]);
    if (R[v] == i) {
      R[v] = UINT_E_MAX;
      if (R[u] == i) {
        matched[u] = matched[v] = 1;
        R[u] = UINT_E_MAX;
        // mark edge
        E[i] = make_tuple(get<0>(E[i]) |= TOP_BIT, get<1>(E[i]), get<2>(E[i]));
        return 1;
      }
    } else if (R[u] == i)
      R[u] = UINT_E_MAX;
    return 0;
  }
};

size_t hash_to_range(size_t hsh, size_t range) { return hsh & range; }

size_t key_for_pair(uintE k1, uintE k2, pbbs::random rnd) {
  size_t l = std::min(k1, k2);
  size_t r = std::max(k1, k2);
  size_t key = (l << 32) + r;
  return rnd.ith_rand(key);
}

template <template <class W> class vertex, class W>
edge_array<W> get_all_edges(graph<vertex<W>>& G, bool* matched,
                            pbbs::random rnd) {
  using edge = tuple<uintE, uintE, W>;
  auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    return !(matched[src] || matched[ngh]) && (src < ngh);
  };
  auto E = filter_all_edges(G, pred);

  auto e_arr = E.E;
  timer perm_t;
  perm_t.start();
  auto perm = pbbs::random_permutation<uintT>(E.non_zeros);
  auto out = array_imap<edge>(E.non_zeros);
  parallel_for_bc(i, 0, E.non_zeros,
                  (E.non_zeros > pbbs::kSequentialForThreshold), {
                    out[i] = e_arr[perm[i]];  // gather or scatter?
                  });
  E.del();
  E.E = out.get_array();
  perm_t.stop();
  perm_t.reportTotal("permutation time");
  return E;
}

template <template <class W> class vertex, class W>
edge_array<W> get_edges(graph<vertex<W>>& G, size_t k, bool* matched,
                        pbbs::random r) {
  using edge = tuple<uintE, uintE, W>;
  size_t m = G.m / 2;  // assume sym
  bool finish = (m <= k);

  cout << "Threshold, using m = " << m << endl;
  size_t range = pbbs::log2_up(G.m);
  range = 1L << range;
  range -= 1;

  size_t threshold = k;
  auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    size_t hash_val = hash_to_range(key_for_pair(src, ngh, r), range);
    if ((src > ngh) || matched[src] || matched[ngh]) {
      return 1;  // pack out, not in edgearr
    } else if (hash_val < threshold || finish) {
      return 2;  // pack out, returned in edgearr
    }
    return 0;  // keep in graph, not in edgearr
  };
  auto E = filter_edges(G, pred);

  // permute the retrieved edges
  auto e_arr = E.E;
  timer perm_t;
  perm_t.start();
  auto perm = pbbs::random_permutation<uintT>(E.non_zeros);
  auto out = array_imap<edge>(E.non_zeros);
  parallel_for_bc(i, 0, E.non_zeros,
                  (E.non_zeros > pbbs::kSequentialForThreshold), {
                    out[i] = e_arr[perm[i]];  // gather or scatter?
                  });
  E.del();
  E.E = out.get_array();
  perm_t.stop();
  perm_t.reportTotal("permutation time");
  return E;
}

};  // namespace mm

// Runs a constant number of filters on the whole graph, and runs the
// prefix-based algorithm on them. Finishes off the rest of the graph with the
// prefix-based algorithm.
template <template <class W> class vertex, class W>
auto MaximalMatching(graph<vertex<W>>& G) {
  using edge = tuple<uintE, uintE, W>;

  timer mt;
  mt.start();
  size_t n = G.n;
  auto r = pbbs::default_random;

  auto R = array_imap<uintE>(n, [&](size_t i) { return UINT_E_MAX; });
  auto matched = array_imap<bool>(n, [&](size_t i) { return false; });

  size_t k = ((3 * G.n) / 2);
  auto matching = dyn_arr<edge>(n);

  size_t round = 0;
  timer gete;
  timer eff;
  while (G.m > 0) {
    gete.start();
    auto e_arr = (round < mm::n_filter_steps)
                     ? mm::get_edges(G, k, matched.start(), r)
                     : mm::get_all_edges(G, matched.start(), r);

    auto eim = make_in_imap<edge>(e_arr.non_zeros,
                                  [&](size_t i) { return e_arr.E[i]; });
    gete.stop();

    cout << "Got: " << e_arr.non_zeros << " edges "
         << " G.m is now: " << G.m << endl;
    mm::matchStep<W> mStep(e_arr.E, R.start(), matched.start());
    eff.start();
    eff_for<uintE>(mStep, 0, e_arr.non_zeros, 50, 0, G.n);
    eff.stop();

    auto e_added =
        pbbs::filter(eim, [](edge e) { return get<0>(e) & mm::TOP_BIT; });
    auto sizes = array_imap<size_t>(e_added.size());
    parallel_for_bc(i, 0, e_added.size(),
                    (e_added.size() > pbbs::kSequentialForThreshold), {
                      const auto& e = e_added[i];
                      uintE u = get<0>(e) & mm::VAL_MASK;
                      uintE v = get<1>(e) & mm::VAL_MASK;
                      uintE deg_u = G.V[u].getOutDegree();
                      uintE deg_v = G.V[v].getOutDegree();
                      G.V[u].setOutDegree(0);
                      G.V[v].setOutDegree(0);
                      sizes[i] = deg_u + deg_v;
                    });
    size_t total_size = pbbs::reduce_add(sizes);
    G.m -= total_size;
    cout << "removed: " << total_size << " many edges" << endl;

    matching.copyIn(e_added, e_added.size());

    round++;
    r = r.next();
  }
  cout << "matching size = " << matching.size << endl;
  auto ret = make_array_imap<edge>(matching.A, matching.size);
  ret.allocated = true;
  mt.stop();
  mt.reportTotal("Matching time");
  eff.reportTotal("eff for time");
  gete.reportTotal("get edges time");
  return std::move(ret);
}

template <template <class W> class vertex, class W, class Seq>
void verify_matching(graph<vertex<W>>& G, Seq& matching) {
  size_t n = G.n;
  auto ok = array_imap<bool>(n, [](size_t i) { return 1; });
  auto matched = array_imap<uintE>(n, [](size_t i) { return 0; });

  // Check that this is a valid matching
  parallel_for_bc(i, 0, matching.size(),
                  (matching.size() > pbbs::kSequentialForThreshold), {
                    const auto& edge = matching[i];
                    pbbs::write_add(&matched[get<0>(edge)], 1);
                    pbbs::write_add(&matched[get<1>(edge)], 1);
                  });

  bool valid = true;
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
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
  parallel_for_bc(i, 0, n, true, { G.V[i].mapOutNgh(i, map2_f); });

  auto ok_im = make_in_imap<size_t>(n, [&](size_t i) { return ok[i]; });
  size_t n_ok = pbbs::reduce_add(ok_im);
  if (n == n_ok) {
    cout << "Matching OK! matching size is: " << matching.size() << endl;
  } else {
    cout << "Matching invalid---" << (n - n_ok)
         << " vertices saw bad neighborhoods." << endl;
  }
}
