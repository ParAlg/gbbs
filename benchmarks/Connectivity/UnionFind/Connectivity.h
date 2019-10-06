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
#include "pbbslib/random.h"

#include <iostream>
#include <limits.h>
#include <vector>
#include <mutex>
#include <atomic>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <random>


namespace union_find {

/* ================================== CSR templates ================================== */

/* Used if the algorithm either requires a hooks array, or uses the hooks array
 * to emit a spanning forest */
template <class Find, class Unite, class G>
struct UnionFindHookTemplate {
  G& GA;
  Unite& unite;
  Find& find;
  UnionFindHookTemplate(G& GA, Unite& unite, Find& find) : GA(GA), unite(unite), find(find) {}

  pbbs::sequence<uintE> components() {
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
};

/* Used if the algorithm only requires a parents array */
template <class Find, class Unite, class G>
struct UnionFindTemplate {

  G& GA;
  Unite& unite;
  Find& find;
  UnionFindTemplate(G& GA, Unite& unite, Find& find) : GA(GA), unite(unite), find(find) {}

  pbbs::sequence<uintE> components() {
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
};

/* ****************************** Sampling ******************************/

// From gapbs/cc.c, Thanks to S. Beamer + M. Sutton for the well documented
// reference implementation of afforest.
template <class Seq>
typename Seq::value_type sample_frequent_element(Seq& S, uint64_t num_samples=1024) {
  using T = typename Seq::value_type;
  std::unordered_map<T, int> sample_counts(32);
  using kvp_type = typename std::unordered_map<T, int>::value_type;
  // Sample elements from 'S'
  std::mt19937 gen;
  std::uniform_int_distribution<T> distribution(0, S.size() - 1);
  for (T i = 0; i < num_samples; i++) {
    T n = distribution(gen);
    sample_counts[S[n]]++;
  }
  // Find most frequent element in samples (estimate of most frequent overall)
  auto most_frequent = std::max_element(
    sample_counts.begin(), sample_counts.end(),
    [](const kvp_type& a, const kvp_type& b) { return a.second < b.second; });

  float frac_of_graph = static_cast<float>(most_frequent->second) / num_samples;
  std::cout
    << "Skipping largest intermediate component (ID: " << most_frequent->first
    << ", approx. " << (frac_of_graph * 100)
    << "% of the graph)" << std::endl;
  return most_frequent->first;
}

// returns a component labeling
// Based on the implementation in gapbs/cc.c, Thanks to S. Beamer + M. Sutton
// for the well documented reference implementation of afforest.
template <class Find, class Unite, class G>
inline pbbs::sequence<uintE> UnionFindSampleTemplate(G& GA, Unite& unite,
    Find& find, uint32_t neighbor_rounds = 2) {
  using W = typename G::weight_type;
  size_t n = GA.n;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });

  pbbs::random rnd();

  for (uint32_t r=0; r<neighbor_rounds; r++) {
    parallel_for(0, n, [&] (size_t u) {
      auto u_rnd = rnd.fork(u);
      auto u_vtx = GA.get_vertex(u);
      if (u_vtx.getOutDegree() > 0) {
        size_t ngh_idx = u_rnd.rand() % u_vtx.getOutDegree();
        uintE ngh = u_vtx.get_ith_out_neighbor(u, ngh_idx);
        unite(u, ngh, parents);
      }
    }, 512);

    // compress nodes
    parallel_for(0, n, [&] (size_t u) {
      find(u, parents);
    }, 512);
  }

  uintE frequent_comp = sample_frequent_element(parents);

  timer ut; ut.start();
  parallel_for(0, n, [&] (size_t u) {
    // Only process edges for vertices not linked to the main component
    // note that this is safe only for undirected graphs. For directed graphs,
    // the in-edges must also be explored for all vertices.
    if (parents[u] != frequent_comp) {
      auto map_f = [&] (uintE u, uintE v, const W& wgh) {
        if (u < v) {
          unite(u, v, parents);
        }
      };
      GA.get_vertex(u).mapOutNgh(u, map_f); // in parallel
    }
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

}  // namespace union_find
