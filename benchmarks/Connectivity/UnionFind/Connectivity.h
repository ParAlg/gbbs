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
#include "union_find_rules.h"

#include "benchmarks/LowDiameterDecomposition/MPX13/LowDiameterDecomposition.h"

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

template <class Find, class Unite, class G>
struct UFAlgorithm {
  G& GA;
  Unite& unite;
  Find& find;
  UFAlgorithm(G& GA, Unite& unite, Find& find) : GA(GA), unite(unite), find(find) {}

  template <bool provides_frequent_comp>
  void compute_components(pbbs::sequence<uintE>& parents, uintE frequent_comp = UINT_E_MAX) {
    using W = typename G::weight_type;
    size_t n = GA.n;

    timer ut; ut.start();
    parallel_for(0, n, [&] (size_t i) {
      auto map_f = [&] (uintE u, uintE v, const W& wgh) {
        if (u < v) {
          unite(u, v, parents);
        }
      };
      if constexpr (provides_frequent_comp) {
        if (parents[i] == frequent_comp) {
          GA.get_vertex(i).mapOutNgh(i, map_f);
        }
      } else {
        GA.get_vertex(i).mapOutNgh(i, map_f);
      }
    }, 1);
    ut.stop(); debug(ut.reportTotal("union time"));

    timer ft; ft.start();
    parallel_for(0, n, [&] (size_t i) {
      parents[i] = find(i,parents);
    });
    ft.stop(); debug(ft.reportTotal("find time"););
  }
};

}  // namespace union_find
