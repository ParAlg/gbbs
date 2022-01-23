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

#include "gbbs/gbbs.h"
#include "gbbs/helpers/speculative_for.h"

namespace gbbs {
namespace MaximalIndependentSet_rootset {

template <class Graph, class Fl>
inline void verify_mis(Graph& G, Fl& in_mis) {
  using W = typename Graph::weight_type;
  auto d = sequence<uintE>(G.n, (uintE)0);
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    if (!d[ngh]) {
      d[ngh] = 1;
    }
  };
  parallel_for(0, G.n, [&](size_t i) {
    if (in_mis[i]) {
      G.get_vertex(i).out_neighbors().map(map_f);
    }
  });
  parallel_for(0, G.n, [&](size_t i) {
    if (in_mis[i]) {
      assert(!d[i]);
    }
  });
  auto mis_f = [&](size_t i) { return (size_t)in_mis[i]; };
  auto mis_int = parlay::delayed_seq<size_t>(G.n, mis_f);
  size_t mis_size = parlay::reduce(mis_int);
  if (parlay::reduce(d) != (G.n - mis_size)) {
    std::cout << "MaximalIndependentSet incorrect"
              << "\n";
    assert(false);
  }
  std::cout << "MaximalIndependentSet Ok"
            << "\n";
}

template <class P, class W>
struct GetNghs {
  P& p;
  GetNghs(P& p) : p(p) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (p[d] > 0) {
      p[d] = 0;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    auto p_d = p[d];
    if (p_d > 0 && gbbs::atomic_compare_and_swap(&p[d], p_d, 0)) {
      return true;
    }
    return false;
  }
  inline bool cond(uintE d) { return (p[d] > 0); }
};

template <class Graph, class VS, class P>
inline vertexSubset get_nghs(Graph& G, VS& vs, P& p) {
  using W = typename Graph::weight_type;
  return neighbor_map(G, vs, GetNghs<P, W>(p));
}

inline bool hash_lt(const uintE& src, const uintE& ngh) {
  uint32_t src_h = parlay::hash32(src);
  uint32_t ngh_h = parlay::hash32(ngh);
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
      return (gbbs::fetch_and_add(&p[d], -1) == 1);
    }
    return false;
  }
  inline bool cond(uintE d) { return (p[d] > 0); }
};

template <class Graph>
inline sequence<bool> MaximalIndependentSet(Graph& G) {
  using W = typename Graph::weight_type;
  timer init_t;
  init_t.start();
  size_t n = G.n;

  // compute the priority DAG
  auto priorities = sequence<intE>(n);  // why intE?
  auto perm = parlay::random_permutation<uintE>(n);
  parallel_for(0, n, 1, [&](size_t i) {
    uintE our_pri = perm[i];
    auto count_f = [&](uintE src, uintE ngh, const W& wgh) {
      uintE ngh_pri = perm[ngh];
      return ngh_pri < our_pri;
    };
    priorities[i] = G.get_vertex(i).out_neighbors().count(count_f);
  });
  init_t.stop();
  debug(init_t.next("init"););

  // compute the initial rootset
  auto zero_f = [&](size_t i) { return priorities[i] == 0; };
  auto zero_map = parlay::delayed_seq<bool>(n, zero_f);
  auto init = parlay::pack_index<uintE>(zero_map);
  auto roots = vertexSubset(n, std::move(init));

  auto in_mis = sequence<bool>(n, false);
  size_t finished = 0;
  size_t rounds = 0;
  while (finished != n) {
    assert(roots.size() > 0);
    std::cout << "## round = " << rounds << " size = " << roots.size()
              << " remaining = " << (n - finished) << "\n";

    // set the roots in the MaximalIndependentSet
    vertexMap(roots, [&](uintE v) { in_mis[v] = true; });

    // compute neighbors of roots that are still live using nghMap
    auto removed = get_nghs(G, roots, priorities);
    vertexMap(removed, [&](uintE v) { priorities[v] = 0; });
    std::cout << "## removed: " << removed.size() << " many vertices"
              << std::endl;

    // compute the new roots: neighbors of removed that have their priorities
    // set to 0 after eliminating all nodes in removed
    intE* pri = priorities.begin();
    timer nr;
    nr.start();
    auto new_roots =
        edgeMap(G, removed, mis_f<W>(pri, perm.begin()), -1, sparse_blocked);
    nr.stop();
    nr.next("## new roots time");

    // update finished with roots and removed. update roots.
    finished += roots.size();
    finished += removed.size();

    roots = std::move(new_roots);
    rounds++;
  }
  return in_mis;
}
}  // namespace MaximalIndependentSet_rootset

namespace MaximalIndependentSet_spec_for {
// For each vertex:
//   Flags = 0 indicates undecided
//   Flags = 1 indicates chosen
//   Flags = 2 indicates a neighbor is chosen
template <class Graph>
struct MaximalIndependentSetstep {
  using W = typename Graph::weight_type;
  char* FlagsNext;
  char* Flags;
  Graph& G;

  MaximalIndependentSetstep(char* _PF, char* _F, Graph& _G)
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
    auto id = std::make_tuple(0, 0);
    auto monoid = parlay::make_monoid(red_f, id);
    auto res = G.get_vertex(i).out_neighbors().reduce(map_f, monoid);
    if (std::get<0>(res) > 0) {
      FlagsNext[i] = 2;
    } else if (std::get<1>(res) > 0) {
      FlagsNext[i] = 0;
    }
    return 1;
  }

  bool commit(intT i) { return (Flags[i] = FlagsNext[i]) > 0; }
};

template <class Graph>
inline sequence<char> MaximalIndependentSet(Graph& G) {
  size_t n = G.n;
  auto Flags = sequence<char>::from_function(n, [&](size_t i) { return 0; });
  auto FlagsNext = sequence<char>(n);
  auto mis =
      MaximalIndependentSetstep<Graph>(FlagsNext.begin(), Flags.begin(), G);
  eff_for<uintE>(mis, 0, n, 50);
  return Flags;
}
};  // namespace MaximalIndependentSet_spec_for

template <class Graph, class Seq>
inline void verify_MaximalIndependentSet(Graph& G, Seq& mis) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  auto ok = sequence<bool>::from_function(n, [&](size_t i) { return 1; });
  parallel_for(0, n, [&](size_t i) {
    auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      return mis[ngh];
    };
    size_t ct = G.get_vertex(i).out_neighbors().count(pred);
    ok[i] = (mis[i]) ? (ct == 0) : (ct > 0);
  });
  auto ok_f = [&](size_t i) { return ok[i]; };
  auto ok_imap = parlay::delayed_seq<size_t>(n, ok_f);
  size_t n_ok = parlay::reduce(ok_imap);
  if (n_ok == n) {
    std::cout << "valid MaximalIndependentSet"
              << "\n";
  } else {
    std::cout << "invalid MaximalIndependentSet, " << (n - n_ok)
              << " vertices saw bad neighborhoods"
              << "\n";
  }
}
}  // namespace gbbs
