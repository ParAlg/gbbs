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

#include "ligra/bridge.h"
#include "ligra/ligra.h"

#include <iostream>
#include <limits.h>
#include <vector>
#include <mutex>
#include <atomic>


namespace union_find {

/* ================================== CSR templates ================================== */

// returns a component labeling
template <class Find, class Unite, class G>
inline pbbs::sequence<uintE> UnionFindHookTemplate(G& GA, Unite& unite, Find& find) {
  using W = typename G::weight_type;
  size_t n = GA.n;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });
  auto hooks = pbbs::sequence<uintE>(n, [&] (size_t i) { return UINT_E_MAX; });

  timer ut; ut.start();
  parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (uintE u, uintE v, const W& wgh) {
      if (u < v) {
        unite(u, v, parents, hooks);
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f); // in parallel
  }, 1);
  ut.stop(); debug(ut.reportTotal("union time"));

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    parents[i] = find(i,parents);
  });
  ft.stop(); debug(ft.reportTotal("find time"););
  // filter UINT_E_MAX from hooks if spanning forest edges are desired.
  return parents;
}

// returns a component labeling
template <class Find, class Unite, class G>
inline pbbs::sequence<uintE> UnionFindTemplate(G& GA, Unite& unite, Find& find) {
  using W = typename G::weight_type;
  size_t n = GA.n;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });

  timer ut; ut.start();
  parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (uintE u, uintE v, const W& wgh) {
      if (u < v) {
        unite(u, v, parents);
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f); // in parallel
  }, 1);
  ut.stop(); debug(ut.reportTotal("union time"));

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    parents[i] = find(i,parents);
  });
  ft.stop(); debug(ft.reportTotal("find time"););
  // filter UINT_E_MAX from hooks if spanning forest edges are desired.
  return parents;
}

// returns a component labeling
template <class Find, class Unite, class G>
inline pbbs::sequence<uintE> UnionFindRemTemplate(G& GA, Unite& unite, Find& find) {
  using W = typename G::weight_type;
  size_t n = GA.n;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });
  auto locks = pbbs::sequence<std::mutex>(n);

  timer ut; ut.start();
  parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (uintE u, uintE v, const W& wgh) {
      if (u < v) {
        unite(u, v, parents, locks);
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f); // in parallel
  }, 1);
  ut.stop(); debug(ut.reportTotal("union time"));

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    parents[i] = find(i,parents);
  });
  ft.stop(); debug(ft.reportTotal("find time"););
  // filter UINT_E_MAX from hooks if spanning forest edges are desired.
  return parents;
}


/* ================================== COO templates ================================== */

// returns a component labeling
template <class Find, class Unite, class W>
inline pbbs::sequence<uintE> UnionFindHookTemplate_coo(edge_array<W>& G, Unite& unite, Find& find) {
  size_t n = G.num_rows;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });
  auto hooks = pbbs::sequence<uintE>(n, [&] (size_t i) { return UINT_E_MAX; });

  timer ut; ut.start();
  parallel_for(0, G.non_zeros, [&] (size_t i) {
    uintE u, v; W wgh;
    std::tie(u, v, wgh) = G.E[i];
    if (u < v) {
      unite(u, v, parents, hooks);
    }
  }, 512);
  ut.stop(); debug(ut.reportTotal("union time"));

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    parents[i] = find(i,parents);
  });
  ft.stop(); debug(ft.reportTotal("find time"););
  // filter UINT_E_MAX from hooks if spanning forest edges are desired.
  return parents;
}

template <class Find, class Unite, class W>
inline pbbs::sequence<uintE> UnionFindTemplate_coo(edge_array<W>& G, Unite& unite, Find& find) {
  size_t n = G.num_rows;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });

  timer ut; ut.start();
  parallel_for(0, G.non_zeros, [&] (size_t i) {
    uintE u, v; W wgh;
    std::tie(u, v, wgh) = G.E[i];
    if (u < v) {
      unite(u, v, parents);
    }
  }, 512);
  ut.stop(); debug(ut.reportTotal("union time"));

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    parents[i] = find(i,parents);
  });
  ft.stop(); debug(ft.reportTotal("find time"););
  return parents;
}

template <class Find, class Unite, class W>
inline pbbs::sequence<uintE> UnionFindRemTemplate_coo(edge_array<W>& G, Unite& unite, Find& find) {
  size_t n = G.num_rows;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });
  auto locks = pbbs::sequence<std::mutex>(n);

  timer ut; ut.start();
  parallel_for(0, G.non_zeros, [&] (size_t i) {
    uintE u, v; W wgh;
    std::tie(u, v, wgh) = G.E[i];
    if (u < v) {
      unite(u, v, parents, locks);
    }
  }, 512);
  ut.stop(); debug(ut.reportTotal("union time"));

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    parents[i] = find(i,parents);
  });
  ft.stop(); debug(ft.reportTotal("find time"););
  return parents;
}

}  // namespace union_find
