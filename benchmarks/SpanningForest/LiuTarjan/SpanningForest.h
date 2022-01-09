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
#include "gbbs/bridge.h"
#include "gbbs/gbbs.h"
#include "liu_tarjan_rules.h"

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
namespace lt {

template <LiuTarjanConnectOption connect_option>
auto get_connect_function() {
  if
    constexpr(connect_option == parent_connect) {
      return primitives::parent_connect;
    }
  else if
    constexpr(connect_option == extended_connect) {
      return primitives::extended_connect;
    }
  else {
    abort();
  }
}

template <LiuTarjanUpdateOption update_option>
auto get_update_function() {
  /* only root update is supported */
  if
    constexpr(update_option == root_update) { return primitives::root_update; }
  else {
    abort();
  }
}

template <LiuTarjanShortcutOption shortcut_option>
auto get_shortcut_function() {
  if
    constexpr(shortcut_option == shortcut) { return primitives::shortcut; }
  else if
    constexpr(shortcut_option == full_shortcut) {
      return primitives::root_shortcut;
    }
  else {
    abort();
  }
}

template <class Connect, LiuTarjanConnectOption connect_option, class Update,
          LiuTarjanUpdateOption update_option, class Shortcut,
          LiuTarjanShortcutOption shortcut_option, class Graph>
struct LiuTarjanAlgorithm {
  Graph& GA;
  size_t n;
  Connect& connect;
  Update& update;
  Shortcut& shortcut;
  using message = lt::message;
  sequence<message> messages;
  LiuTarjanAlgorithm(Graph& GA, size_t n, Connect& connect, Update& update,
                     Shortcut& shortcut)
      : GA(GA), n(n), connect(connect), update(update), shortcut(shortcut) {
    static_assert(update_option ==
                  root_update); /* only works for root_update algorithms */
  }

  void initialize(sequence<parent>& P, sequence<edge>& Edges) {
    messages = sequence<message>::uninitialized(P.size());
    parallel_for(0, n, [&](size_t i) { messages[i] = message(i, empty_edge); });
  }

  template <SamplingOption sampling_option>
  void compute_spanning_forest(sequence<parent>& P, sequence<edge>& Edges,
                               parent frequent_comp = UINT_E_MAX) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    auto parents_changed = true;
    while (parents_changed) {
      parents_changed = false;

      // connect
      parallel_for(0, n, [&](size_t i) {
        auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
          bool updated = connect(u, v, P, messages);
          if (updated && !parents_changed) {
            parents_changed = true;
          }
        };
        if
          constexpr(sampling_option != no_sampling) {
            if (P[i] != frequent_comp) {
              GA.get_vertex(i).out_neighbors().map(map_f);
            }
          }
        else {
          GA.get_vertex(i).out_neighbors().map(map_f);
        }
      });

      // Update
      parallel_for(0, n, [&](size_t u) {
        update(u, P, messages, Edges); /* root_update */
      });

      // Shortcut
      parallel_for(0, n, [&](size_t u) { shortcut(u, P); });
    }
  }
};

}  // namespace lt
}  // namespace gbbs
