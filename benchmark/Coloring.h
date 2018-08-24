// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

#pragma once

#include "bucket.h"
#include "ligra.h"

#include "lib/index_map.h"
#include "lib/random_shuffle.h"
#include "lib/sparse_table.h"

namespace coloring {
template <template <typename W> class vertex, class W, class Seq>
uintE color(graph<vertex<W> >& GA, uintE v, Seq& colors) {
  uintE deg = GA.V[v].getOutDegree();
  if (deg > 0) {
    bool* bits;
    bool s_bits[1000];
    if (deg > 1000)
      bits = newA(bool, deg);
    else
      bits = (bool*)s_bits;

    granular_for(i, 0, deg, (deg > 2000), { bits[i] = 0; });
    auto map_f = wrap_f<W>([&](uintE src, uintE ngh) {
      uintE color = colors[ngh];
      if (color < deg) {
        bits[color] = 1;
      }
    });
    GA.V[v].mapOutNgh(v, map_f);
    auto im = make_in_imap<uintE>(
        deg, [&](size_t i) { return (bits[i] == 0) ? (uintE)i : UINT_E_MAX; });
    auto min_f = [](uintE l, uintE r) { return std::min(l, r); };
    uintE color = pbbs::reduce(im, min_f);
    if (deg > 1000) {
      free(bits);
    }
    return (color == UINT_E_MAX) ? (deg + 1) : color;
  }
  return 0;
}
}  // namespace coloring

template <class W>
struct coloring_f {
  intE* p;
  coloring_f(intE* _p) : p(_p) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (s == d) {
      cout << "error" << endl;
      exit(-1);
    }
    p[d]--;
    return p[d] == 0;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    if (s == d) {
      cout << "error" << endl;
      exit(-1);
    }
    return (pbbs::xadd(&p[d], -1) == 1);
  }
  inline bool cond(uintE d) { return (p[d] > 0); }
};

template <template <typename W> class vertex, class W>
array_imap<uintE> Coloring(graph<vertex<W> >& GA, bool lf=false) {
  timer color_t;
  color_t.start();
  timer initt;
  initt.start();
  const size_t n = GA.n;
  const size_t m = GA.m;

  // For each vertex count the number of out-neighbors with log-degree >= us
  auto priorities = array_imap<intE>(n);
  auto colors = array_imap<uintE>(n, [](size_t i) { return UINT_E_MAX; });

  if (lf) {
    cout << "Running LF" << endl;
    // LF heuristic
    auto P = pbbs::random_permutation<uintE>(n);
    parallel_for(size_t i=0; i<n; i++) {
      uintE our_deg = GA.V[i].getOutDegree();
      uintE i_p = P[i];
      auto count_f = wrap_f<W>([&] (uintE src, uintE ngh) {
        uintE ngh_deg = GA.V[ngh].getOutDegree();
        return (ngh_deg > our_deg) || ((ngh_deg == our_deg) && P[ngh] < i_p);
      });
      priorities[i] = GA.V[i].countOutNgh(i, count_f);
    }
  } else {
    cout << "Running LLF" << endl;
    // LLF heuristic
    auto P = pbbs::random_permutation<uintE>(n);
    parallel_for(size_t i = 0; i < n; i++) {
      uintE our_deg = pbbs::log2_up(GA.V[i].getOutDegree());
      uintE i_p = P[i];
      // breaks ties using P
      auto count_f = wrap_f<W>([&](uintE src, uintE ngh) {
        uintE ngh_deg = pbbs::log2_up(GA.V[ngh].getOutDegree());
        return (ngh_deg > our_deg) || ((ngh_deg == our_deg) && P[ngh] < i_p);
      });
      priorities[i] = GA.V[i].countOutNgh(i, count_f);
    }
  }

  auto zero_map =
      make_in_imap<bool>(n, [&](size_t i) { return priorities[i] == 0; });
  auto init = pbbs::pack_index<uintE>(zero_map);
  auto roots = vertexSubset(n, init.size(), init.get_array());
  initt.reportTotal("init time");

  size_t finished = 0, rounds = 0;
  while (finished != n) {
    assert(roots.size() > 0);
    finished += roots.size();

    // color the rootset
    vertexMap(roots,
              [&](uintE v) { colors[v] = coloring::color(GA, v, colors); });

    // compute the new rootset
    auto new_roots = edgeMap(GA, roots, coloring_f<W>(priorities.start()), -1,
                             sparse_blocked);
    roots.del();
    roots = new_roots;
    rounds++;
  }
  cout << "Total rounds = " << rounds << endl;
  return std::move(colors);
}

template <template <typename W> class vertex, class W, class Seq>
void verify_coloring(graph<vertex<W>>& G, Seq& colors) {
  size_t n = G.n;
  auto ok = array_imap<bool>(n);
  parallel_for(size_t i=0; i<n; i++) {
    uintE src_color = colors[i];
    auto pred = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
      uintE ngh_color = colors[ngh];
      return src_color == ngh_color;
    };
    size_t ct = G.V[i].countOutNgh(i, pred);
    ok[i] = (ct > 0);
  }
  auto im = make_in_imap<size_t>(n, [&] (size_t i) { return (size_t)ok[i]; });
  size_t ct = pbbs::reduce_add(im);
  cout << "ct = " << ct << endl;
  if (ct > 0) {
    cout << "Invalid coloring" << endl;
  } else {
    cout << "Valid coloring" << endl;
  }
}
