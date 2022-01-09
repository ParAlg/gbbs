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
    constexpr(connect_option == simple_connect) { return primitives::connect; }
  else if
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
  if
    constexpr(update_option == simple_update) {
      return primitives::simple_update;
    }
  else if
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

template <LiuTarjanAlterOption alter_option>
auto get_alter_function() {
  return primitives::alter;
}

template <class Connect, LiuTarjanConnectOption connect_option, class Update,
          LiuTarjanUpdateOption update_option, class Shortcut,
          LiuTarjanShortcutOption shortcut_option, class Alter,
          LiuTarjanAlterOption alter_option, class Graph>
struct LiuTarjanAlgorithm {
  using W = typename Graph::weight_type;

  Graph& GA;
  size_t n;
  Connect& connect;
  Update& update;
  Shortcut& shortcut;
  Alter& alter;
  sequence<uintE> messages;
  sequence<bool> flags;
  LiuTarjanAlgorithm(Graph& GA, size_t n, Connect& connect, Update& update,
                     Shortcut& shortcut, Alter& alter)
      : GA(GA),
        n(n),
        connect(connect),
        update(update),
        shortcut(shortcut),
        alter(alter) {}

  void initialize(sequence<parent>& P) {
    messages = sequence<uintE>(P.size());
    parallel_for(0, n, [&](size_t i) { messages[i] = i; });
    flags = sequence<bool>(P.size(), false);
  }

  template <SamplingOption sampling_option>
  void compute_components(sequence<parent>& P,
                          parent frequent_comp = UINT_E_MAX) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    auto parents_changed = true;
    while (parents_changed) {
      parents_changed = false;

      uintE granularity;
      constexpr bool provides_frequent_comp = sampling_option != no_sampling;
      if
        constexpr(provides_frequent_comp) {
          granularity = 512;
          std::cout << "# provides frequent comp" << std::endl;
        }
      else {
        granularity = 1;
      }
      // connect
      parallel_for(0, n,
                   [&](size_t i) {
                     auto map_f = [&](const uintE& u, const uintE& v,
                                      const W& wgh) {
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
                   },
                   granularity);

      // Update
      parallel_for(0, n, [&](size_t u) { update(u, P, messages); });

      // Shortcut
      parallel_for(0, n, [&](size_t u) { shortcut(u, P); });

      // Note that alter is not applied since we do not process edges
      // in COO, and performing edge modification/filtering on CSR
      // would be prohibitively costly.
    }
  }

  template <bool reorder_batch, class Seq>
  void process_batch(sequence<parent>& parents, Seq& updates) {
    /* Must reorder updates since queries are not linearizable o.w. */
    auto ret = reorder_updates(updates);
    auto reordered_updates = ret.first;
    size_t update_end = ret.second;
    auto insertions = reordered_updates.cut(0, update_end);
    auto queries = reordered_updates.cut(update_end, updates.size());

    using edge = std::pair<uintE, uintE>;
    auto inserts =
        sequence<edge>::from_function(insertions.size(), [&](size_t i) {
          auto[u, v, typ] = insertions[i];
          (void)typ;
          return std::make_pair(u, v);
        });

    auto parents_changed = true;
    auto& P = parents;
    size_t round = 0;
    while (parents_changed) {
      round++;
      parents_changed = false;
      std::cout << "# round = " << round << std::endl;

      // Parent-Connect
      timer pc;
      pc.start();
      parallel_for(0, inserts.size(), [&](size_t i) {
        auto[u, v] = inserts[i];
        bool updated = connect(u, v, P, messages);
        if (updated && !parents_changed) {
          parents_changed = true;
        }
      });
      pc.stop();
      pc.next("# pc time");

      // Update local neighborhoods
      timer ut;
      ut.start();
      parallel_for(0, inserts.size(), [&](size_t i) {
        auto[u, v] = inserts[i];
        auto p_u = P[u];
        auto p_v = P[v];
        if (flags[u] == false &&
            gbbs::atomic_compare_and_swap(&flags[u], false, true)) {
          update(u, P, messages);
        }
        if (flags[v] == false &&
            gbbs::atomic_compare_and_swap(&flags[v], false, true)) {
          update(v, P, messages);
        }
        if (flags[p_u] == false &&
            gbbs::atomic_compare_and_swap(&flags[p_u], false, true)) {
          update(p_u, P, messages);
        }
        if (flags[p_v] == false &&
            gbbs::atomic_compare_and_swap(&flags[p_v], false, true)) {
          update(p_v, P, messages);
        }
      });
      parallel_for(0, inserts.size(), [&](size_t i) {
        auto[u, v] = inserts[i];
        if (flags[u]) {
          flags[u] = false;
        }
        if (flags[P[u]]) {
          flags[P[u]] = false;
        }
        if (flags[v]) {
          flags[v] = false;
        }
        if (flags[P[v]]) {
          flags[P[v]] = false;
        }
      });
      ut.stop();
      ut.next("# update time");

      // Shortcut
      timer sc;
      sc.start();
      parallel_for(0, inserts.size(), [&](size_t i) {
        auto[u, v] = inserts[i];
        if (flags[u] == false &&
            gbbs::atomic_compare_and_swap(&flags[u], false, true)) {
          shortcut(u, P);
          messages[u] = P[u];
        }
        if (flags[v] == false &&
            gbbs::atomic_compare_and_swap(&flags[v], false, true)) {
          shortcut(v, P);
          messages[v] = P[v];
        }
      });
      sc.stop();
      sc.next("# shortcut time");

      parallel_for(0, inserts.size(), [&](size_t i) {
        auto[u, v] = inserts[i];
        if (flags[u]) {
          flags[u] = false;
        }
        if (flags[v]) {
          flags[v] = false;
        }
      });

      timer at;
      at.start();
      if
        constexpr(alter_option != no_alter) {
          // Alter
          constexpr uintE null_comp = UINT_E_MAX - 1;
          constexpr edge nullary_edge = std::make_pair(null_comp, null_comp);
          parallel_for(0, inserts.size(), [&](size_t i) {
            auto[u, v] = inserts[i];
            uintE p_u, p_v;
            p_u = P[u];
            p_v = P[v];
            if (p_u == p_v) {
              inserts[i] = std::make_pair(null_comp, null_comp);
            } else {
              inserts[i] = std::make_pair(p_u, p_v);
            }
          });

          auto new_inserts = parlay::filter(
              inserts, [&](const edge& e) { return e != nullary_edge; });
          inserts = new_inserts;
        }
      at.stop();
      at.next("# alter time");
    }

    // Process queries
    parallel_for(0, queries.size(), [&](size_t i) {
      auto[u, v, utype] = updates[i];
      (void)utype;

      while (P[u] != P[P[u]]) {
        P[u] = P[P[u]];
      } /* found p_u */

      while (P[v] != P[P[v]]) {
        P[v] = P[P[v]];
      } /* found p_v */
    });
  }
};

template <class Connect, LiuTarjanConnectOption connect_option, class Update,
          LiuTarjanUpdateOption update_option, class Shortcut,
          LiuTarjanShortcutOption shortcut_option, class Alter,
          LiuTarjanAlterOption alter_option, class W>
struct LiuTarjanAlgorithmCOO {
  using edge = std::tuple<uintE, uintE, W>;
  sequence<edge> graph;
  size_t n;

  Connect& connect;
  Update& update;
  Shortcut& shortcut;
  Alter& alter;

  sequence<uintE> messages;

  LiuTarjanAlgorithmCOO(sequence<edge>&& graph, size_t n, Connect& connect,
                        Update& update, Shortcut& shortcut, Alter& alter)
      : graph(std::move(graph)),
        n(n),
        connect(connect),
        update(update),
        shortcut(shortcut),
        alter(alter) {}

  void initialize(sequence<parent>& P) {
    messages = sequence<uintE>(P.size());
    parallel_for(0, n, [&](size_t i) { messages[i] = i; });
  }

  void my_alter(sequence<parent>& P) {
    parallel_for(0, graph.size(), [&](size_t i) {
      edge& e = graph[i];
      e = std::make_tuple(
          (std::get<0>(e) == largest_comp) ? largest_comp : P[std::get<0>(e)],
          (std::get<1>(e) == largest_comp) ? largest_comp : P[std::get<1>(e)],
          W());
    });
  }

  bool my_connect(sequence<parent>& P) {
    bool parents_changed = false;
    parallel_for(0, graph.size(), [&](size_t i) {
      const edge& e = graph[i];
      if (std::get<0>(e) != std::get<1>(e)) {
        bool updated = connect(std::get<0>(e), std::get<1>(e), P, messages);
        if (updated && !parents_changed) {
          parents_changed = true;
        }
      }
    });
    return parents_changed;
  }

  void my_update(sequence<parent>& P) {
    parallel_for(0, n, [&](size_t u) { update(u, P, messages); });
  }

  void my_shortcut(sequence<parent>& P) {
    parallel_for(0, n, [&](size_t u) { shortcut(u, P); });
  }

  template <SamplingOption sampling_option>
  void compute_components(sequence<parent>& P,
                          parent frequent_comp = UINT_E_MAX) {
    auto parents_changed = true;
    static_assert(alter_option != no_alter);
    my_alter(P);
    while (parents_changed) {
      parents_changed = my_connect(P);
      my_update(P);
      my_shortcut(P);
      my_alter(P);
    }
  }
};

template <class Graph>
struct StergiouAlgorithm {
  Graph& GA;
  size_t n;
  StergiouAlgorithm(Graph& GA, size_t n) : GA(GA), n(n) {}

  void initialize(sequence<parent>& P) {}

  template <SamplingOption sampling_option>
  void compute_components(sequence<parent>& parents,
                          parent frequent_comp = UINT_E_MAX) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    sequence<parent> previous_parents(n);

    auto parents_changed = true;
    size_t rounds = 0;
    while (parents_changed) {
      parents_changed = false;
      if (rounds > 0) {
        parallel_for(0, n, [&](size_t i) { previous_parents[i] = parents[i]; });
      }

      parallel_for(0, n, [&](size_t i) {
        auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
          parent parent_u = previous_parents[u];
          parent parent_v = previous_parents[v];
          bool updated = false;
          if (parents[v] > parent_u) {
            updated |=
                gbbs::write_min(&parents[v], parent_u, std::less<parent>());
          }
          if (parents[u] > parent_v) {
            updated |=
                gbbs::write_min(&parents[u], parent_v, std::less<parent>());
          }
          if (updated && !parents_changed) {
            parents_changed = true;
          }
        };
        if
          constexpr(sampling_option != no_sampling) {
            if (parents[i] != frequent_comp) {
              GA.get_vertex(i).out_neighbors().map(map_f);
            }
          }
        else {
          GA.get_vertex(i).out_neighbors().map(map_f);
        }
      });

      // Shortcut
      parallel_for(0, n,
                   [&](size_t u) { lt::primitives::shortcut(u, parents); });
    }
  }
};

}  // namespace lt
}  // namespace gbbs
