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

#include "benchmarks/Connectivity/common.h"
#include "gbbs/gbbs.h"
#include "union_find_rules.h"

#include <limits.h>
#include <algorithm>
#include <atomic>
#include <iostream>
#include <mutex>
#include <random>
#include <unordered_map>
#include <vector>
#include <vector>

namespace gbbs {
namespace union_find {

/* ================================== CSR templates
 * ================================== */

template <class Find, class Unite, class G>
struct UFAlgorithm {
  G& GA;
  Unite& unite;
  Find& find;
  UFAlgorithm(G& GA, Unite& unite, Find& find)
      : GA(GA), unite(unite), find(find) {}

  void initialize(sequence<parent>& P) {}

  template <SamplingOption sampling_option>
  void compute_components(sequence<parent>& parents,
                          uintE frequent_comp = UINT_E_MAX) {
    using W = typename G::weight_type;
    constexpr bool provides_frequent_comp = sampling_option != no_sampling;
    size_t n = GA.n;
    sequence<parent> clusters;
    uintE granularity;
    if
      constexpr(provides_frequent_comp) {
        clusters = parents;
        granularity = 512;
        std::cout << "# provides frequent comp" << std::endl;
      }
    else {
      granularity = 1;
    }
    std::cout << "# frequent_comp = " << frequent_comp
              << " gran = " << granularity << std::endl;

    timer ut;
    ut.start();
    parallel_for(0, n,
                 [&](size_t i) {
                   auto map_f = [&](uintE u, uintE v, const W& wgh) {
                     if
                       constexpr(provides_frequent_comp) {
                         unite(u, v, parents);
                       }
                     else {
                       if (v < u) {
                         unite(u, v, parents);
                       }
                     }
                   };
                   if
                     constexpr(provides_frequent_comp) {
                       if (clusters[i] != frequent_comp) {
                         GA.get_vertex(i).out_neighbors().map(map_f);
                       }
                     }
                   else {
                     GA.get_vertex(i).out_neighbors().map(map_f);
                   }
                 },
                 granularity);
    ut.stop();
    ut.next("union time");

    timer ft;
    ft.start();
    parallel_for(0, n, [&](size_t i) { parents[i] = find(i, parents); });
    ft.stop();
    debug(ft.next("find time"););
  }

  template <bool reorder_batch, class Seq>
  void process_batch(sequence<parent>& parents, Seq& updates) {
    if
      constexpr(reorder_batch == true) {
        auto ret = reorder_updates(updates);
        auto reordered_updates = ret.first;
        size_t update_end = ret.second;
        auto insertions = reordered_updates.cut(0, update_end);
        auto queries = reordered_updates.cut(update_end, updates.size());
        /* run updates */
        parallel_for(0, insertions.size(), [&](size_t i) {
          auto[u, v, optype] = insertions[i];
          unite(u, v, parents);
        });
        /* run queries */
        parallel_for(0, queries.size(), [&](size_t i) {
          auto[u, v, optype] = queries[i];
          u = find(u, parents); /* force */
          v = find(v, parents); /* force */
        });
      }
    else {
      /* run queries and updates together */
      parallel_for(0, updates.size(), [&](size_t i) {
        uintE u, v;
        UpdateType optype;
        std::tie(u, v, optype) = updates[i];
        if (optype == query_type) { /* query */
          find(u, parents);         /* force */
          find(v, parents);         /* force */
        } else {                    /* insert */
          unite(u, v, parents);
        }
      });
    }
  }
};

/* default union_find algorithm using find_compress and unite */
template <class Seq>
sequence<parent> find_compress_uf(size_t n, Seq& updates) {
  auto find = find_variants::find_compress;
  auto unite = unite_variants::Unite<decltype(find)>(find);
  auto parents =
      sequence<parent>::from_function(n, [&](size_t i) { return i; });
  parallel_for(0, updates.size(), [&](size_t i) {
    auto[u, v] = updates[i];
    unite(u, v, parents);
  });
  return parents;
}

}  // namespace union_find
}  // namespace gbbs
