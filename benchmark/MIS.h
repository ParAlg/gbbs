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

#include "pbbslib/random_shuffle.h"
#include "pbbslib/sparse_table.h"
#include "pbbslib/speculative_for.h"
#include "ligra.h"

namespace MIS_rootset {

template <template <class W> class vertex, class W, class Fl>
inline void verify_mis(graph<vertex<W>>& GA, Fl& in_mis) {
  auto d = sequence<uintE>(GA.n, (uintE)0);
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    if (!d[ngh]) {
      d[ngh] = 1;
    }
  };
  par_for(0, GA.n, [&] (size_t i) {
    if (in_mis[i]) {
      GA.V[i].mapOutNgh(i, map_f);
    }
  });
  par_for(0, GA.n, [&] (size_t i) {
    if (in_mis[i]) {
      assert(!d[i]);
    }
  });
  auto mis_f = [&](size_t i) { return (size_t)in_mis[i]; };
  auto mis_int =
      make_sequence<size_t>(GA.n, mis_f);
  size_t mis_size = pbbslib::reduce_add(mis_int);
  if (pbbslib::reduce_add(d) != (GA.n - mis_size)) {
    std::cout << "MIS incorrect"
              << "\n";
    assert(false);
  }
  std::cout << "MIS Ok"
            << "\n";
}

template <template <class W> class vertex, class W, class VS, class P>
inline vertexSubset get_nghs(graph<vertex<W>>& GA, VS& vs, P p) {
  vs.toSparse();
  assert(!vs.isDense);
  auto deg_f =  [&](size_t i) { return GA.V[vs.vtx(i)].getOutDegree(); };
  auto deg_im = make_sequence<size_t>(
      vs.size(), deg_f);
  size_t sum_d = pbbslib::reduce_add(deg_im);

  if (sum_d > GA.m / 100) {  // dense forward case
    auto dense = sequence<bool>(GA.n, false);
    auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      if (p(ngh) && !dense[ngh]) {
        dense[ngh] = 1;
      }
    };
    par_for(0, vs.size(), [&] (size_t i) {
      uintE v = vs.vtx(i);
      GA.V[v].mapOutNgh(v, map_f);
    });
    return vertexSubset(GA.n, dense.get_array());
  } else {  // sparse --- iterate, and add nghs satisfying P to a hashtable
    std::cout << "sum_d = " << sum_d << std::endl;
    auto ht = make_sparse_table<uintE, pbbslib::empty>(
        sum_d, std::make_tuple(UINT_E_MAX, pbbslib::empty()),
        [&](const uintE& k) { return pbbslib::hash64(k); });
    vs.toSparse();
    par_for(0, vs.size(), [&] (size_t i) {
      auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
        if (p(ngh)) {
          ht.insert(std::make_tuple(ngh, pbbslib::empty()));
        }
      };
      uintE v = vs.vtx(i);
      GA.V[v].mapOutNgh(v, map_f);
    });
    auto nghs = ht.entries();
    ht.del();
    return vertexSubset(GA.n, nghs.size(), (uintE*)nghs.get_array());
  }
}

inline bool hash_lt(const uintE& src, const uintE& ngh) {
  uint32_t src_h = pbbslib::hash32(src);
  uint32_t ngh_h = pbbslib::hash32(ngh);
  return (src_h < ngh_h) || ((src_h == ngh_h) && src < ngh);
};

template <class W>
struct mis_f {
  intE* p;
  uintE* perm;
  mis_f(intE* _p, uintE* _perm) : p(_p), perm(_perm) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (perm[s] < perm[d]) {
      p[d]--;
      return p[d] == 0;
    }
    return false;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    if (perm[s] < perm[d]) {
      return (pbbslib::xadd(&p[d], -1) == 1);
    }
    return false;
  }
  inline bool cond(uintE d) { return (p[d] > 0); }
};

template <class W>
struct mis_f_2 {
  intE* p;
  mis_f_2(intE* _p) : p(_p) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (hash_lt(s, d)) {
      p[d]--;
      return p[d] == 0;
    }
    return false;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    if (hash_lt(s, d)) {
      return (writeAdd(&p[d], -1) == 0);
    }
    return false;
  }
  inline bool cond(uintE d) { return (p[d] > 0); }  // still live
};

template <template <class W> class vertex, class W>
inline sequence<bool> MIS(graph<vertex<W>>& GA) {
  timer init_t;
  init_t.start();
  size_t n = GA.n;

  // compute the priority DAG
  auto priorities = sequence<intE>(n);
  auto perm = pbbslib::random_permutation<uintE>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintE our_pri = perm[i];
    auto count_f = [&](uintE src, uintE ngh, const W& wgh) {
      uintE ngh_pri = perm[ngh];
      return ngh_pri < our_pri;
    };
    priorities[i] = GA.V[i].countOutNgh(i, count_f);
  });
  init_t.stop();
  init_t.reportTotal("init");

  // compute the initial rootset
  auto zero_f = [&](size_t i) { return priorities[i] == 0; };
  auto zero_map =
      make_sequence<bool>(n, zero_f);
  auto init = pbbslib::pack_index<uintE>(zero_map);
  auto roots = vertexSubset(n, init.size(), init.get_array());

  auto in_mis = sequence<bool>(n, false);
  size_t finished = 0;
  size_t rounds = 0;
  while (finished != n) {
    assert(roots.size() > 0);
    std::cout << "round = " << rounds << " size = " << roots.size()
              << " remaining = " << (n - finished) << "\n";

    // set the roots in the MIS
    vertexMap(roots, [&](uintE v) { in_mis[v] = true; });

    // compute neighbors of roots that are still live
    auto removed = get_nghs(
        GA, roots, [&](const uintE& ngh) { return priorities[ngh] > 0; });
    vertexMap(removed, [&](uintE v) { priorities[v] = 0; });

    // compute the new roots: neighbors of removed that have their priorities
    // set to 0 after eliminating all nodes in removed
    intE* pri = priorities.start();
    timer nr;
    nr.start();
    auto new_roots =
        edgeMap(GA, removed, mis_f<W>(pri, perm.start()), -1, sparse_blocked);
    nr.stop();
    nr.reportTotal("new roots time");

    // update finished with roots and removed. update roots.
    finished += roots.size();
    finished += removed.size();
    removed.del();
    roots.del();

    roots = new_roots;
    rounds++;
  }
  return in_mis;
}
}  // namespace MIS_rootset

namespace MIS_spec_for {
// For each vertex:
//   Flags = 0 indicates undecided
//   Flags = 1 indicates chosen
//   Flags = 2 indicates a neighbor is chosen
template <template <class W> class vertex, class W>
struct MISstep {
  char* FlagsNext;
  char* Flags;
  graph<vertex<W>>& G;

  MISstep(char* _PF, char* _F, graph<vertex<W>>& _G)
      : FlagsNext(_PF), Flags(_F), G(_G) {}

  bool reserve(intT i) {
    // decode neighbor
    FlagsNext[i] = 1;

    auto map_f = [&](const uintE& src, const uintE& ngh,
                     const W& wgh) -> std::tuple<int, int> {
      if (ngh < src) {
        auto fl = Flags[ngh];
        return std::make_tuple(fl == 1, fl == 0);
      }
      return std::make_tuple(0, 0);
    };
    auto red_f = [&](const std::tuple<int, int>& l,
                     const std::tuple<int, int>& r) {
      return std::make_tuple(std::get<0>(l) + std::get<0>(r),
                             std::get<1>(l) + std::get<1>(r));
    };
    using E = std::tuple<int, int>;
    auto id = std::make_tuple(0, 0);
    auto monoid = make_monoid(red_f, id);
    auto res = G.V[i].template reduceOutNgh<E>(i, map_f, monoid);
    if (std::get<0>(res) > 0) {
      FlagsNext[i] = 2;
    } else if (std::get<1>(res) > 0) {
      FlagsNext[i] = 0;
    }
    return 1;
  }

  bool commit(intT i) { return (Flags[i] = FlagsNext[i]) > 0; }
};

template <template <class W> class vertex, class W>
inline sequence<char> MIS(graph<vertex<W>>& GA) {
  size_t n = GA.n;
  auto Flags = sequence<char>(n, [&](size_t i) { return 0; });
  auto FlagsNext = sequence<char>(n);
  auto mis = MISstep<vertex, W>(FlagsNext.start(), Flags.start(), GA);
  eff_for<uintE>(mis, 0, n, 50);
  return Flags;
}
};  // namespace MIS_spec_for

template <template <class W> class vertex, class W, class Seq>
inline void verify_MIS(graph<vertex<W>>& GA, Seq& mis) {
  size_t n = GA.n;
  auto ok = sequence<bool>(n, [&](size_t i) { return 1; });
  par_for(0, n, [&] (size_t i) {
    auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      return mis[ngh];
    };
    size_t ct = GA.V[i].countOutNgh(i, pred);
    ok[i] = (mis[i]) ? (ct == 0) : (ct > 0);
  });
  auto ok_f = [&](size_t i) { return ok[i]; };
  auto ok_imap = make_sequence<size_t>(n, ok_f);
  size_t n_ok = pbbslib::reduce_add(ok_imap);
  if (n_ok == n) {
    std::cout << "valid MIS"
              << "\n";
  } else {
    std::cout << "invalid MIS, " << (n - n_ok)
              << " vertices saw bad neighborhoods"
              << "\n";
  }
}
