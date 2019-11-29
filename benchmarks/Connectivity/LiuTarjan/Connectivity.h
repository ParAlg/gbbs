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
#include "liu_tarjan_rules.h"
#include "benchmarks/Connectivity/common.h"

#include <iostream>
#include <limits.h>
#include <vector>
#include <mutex>
#include <atomic>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <random>


namespace lt {

template <LiuTarjanConnectOption connect_option>
auto get_connect_function() {
  if constexpr (connect_option == simple_connect) {
    return primitives::connect;
  } else if constexpr (connect_option == parent_connect) {
    return primitives::parent_connect;
  } else if constexpr (connect_option == extended_connect) {
    return primitives::extended_connect;
  } else {
    abort();
  }
}

template <LiuTarjanUpdateOption update_option>
auto get_update_function() {
  if constexpr (update_option == simple_update) {
    return primitives::simple_update;
  } else if constexpr (update_option == root_update) {
    return primitives::root_update;
  } else {
    abort();
  }
}

template <LiuTarjanShortcutOption shortcut_option>
auto get_shortcut_function() {
  if constexpr (shortcut_option == shortcut) {
    return primitives::shortcut;
  } else if constexpr (shortcut_option == full_shortcut) {
    return primitives::root_shortcut;
  } else {
    abort();
  }
}

template <LiuTarjanAlterOption alter_option>
auto get_alter_function() {
  return primitives::alter;
}


template <class Connect,
          LiuTarjanConnectOption connect_option,
          class Update,
          LiuTarjanUpdateOption update_option,
          class Shortcut,
          LiuTarjanShortcutOption shortcut_option,
          class Graph>
struct LiuTarjanAlgorithm {
  Graph& GA;
  size_t n;
  Connect& connect;
  Update& update;
  Shortcut& shortcut;
  LiuTarjanAlgorithm(Graph& GA, size_t n, Connect& connect, Update& update, Shortcut& shortcut) :
    GA(GA), n(n), connect(connect), update(update), shortcut(shortcut) {}

  void initialize(pbbs::sequence<parent>& P) {}

  template <SamplingOption sampling_option>
  void compute_components(pbbs::sequence<parent>& P, uintE frequent_comp = UINT_E_MAX) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    auto parents_changed = true;
    while (parents_changed) {
      parents_changed = false;

      // Parent-Connect
      parallel_for(0, n, [&] (size_t i) {
        auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
          bool updated = connect(u, v, P);
          if (updated && !parents_changed) {
            parents_changed = true;
          }
        };
        if constexpr (sampling_option != no_sampling) {
          if (P[i] != frequent_comp) {
            GA.get_vertex(i).mapOutNgh(i, map_f);
          }
        } else {
          GA.get_vertex(i).mapOutNgh(i, map_f);
        }
      });

      // Can skip this step for a regular update
      if constexpr (update_option != simple_update) {
        // Update
        parallel_for(0, n, [&] (size_t u) {
          update(u, P);
        });
      }

      // Shortcut
      parallel_for(0, n, [&] (size_t u) {
        shortcut(u, P);
      });
    }
  }

  template <class Seq>
  void process_batch(pbbs::sequence<parent>& parents, Seq& batch, size_t insert_to_query) {

    auto parents_changed = true;
    auto& P = parents;
    size_t round = 0;
    while (parents_changed) {
      parents_changed = false;

      // Parent-Connect

      parallel_for(0, batch.size(), [&] (size_t i) {
        uintE u, v;
        std::tie(u,v) = batch[i];
        if (i % insert_to_query != 0) { /* is an update */
          bool updated = connect(u, v, P);
          if (updated && !parents_changed) {
            parents_changed = true;
          }
        }
      });

      // Can skip this step for a regular update
      if constexpr (update_option != simple_update) {
        // Update
        parallel_for(0, n, [&] (size_t u) {
          update(u, P);
        });
      }

      /* TODO: separate queries and updates here */
      // Shortcut
      parallel_for(0, batch.size(), [&] (size_t i) {
        uintE u, v;
        std::tie(u,v) = batch[i];
        if (i % insert_to_query != 0) { /* is an update */
          shortcut(u, P);
        } else { /* answer the query, if it is the first round*/
          if (round == 0) {
            while (P[u] != P[P[u]]) {
              P[u] = P[P[u]];
            }
            /* found p_u */

            while (P[v] != P[P[v]]) {
              P[v] = P[P[v]];
            }
            /* found p_v */
          }
        }
      });
      round++;
    }

    // Process queries, no need to skip.
    parallel_for(0, batch.size(), [&] (size_t i) {
      uintE u, v;
      std::tie(u,v) = batch[i];
      if (i % insert_to_query != 0) { /* is an update */
        shortcut(u, P);
      } else { /* answer the query, if it is the first round*/
        if (round == 0) {
          while (P[u] != P[P[u]]) {
            P[u] = P[P[u]];
          }
          /* found p_u */

          while (P[v] != P[P[v]]) {
            P[v] = P[P[v]];
          }
          /* found p_v */
        }
      }
    });

  }

};

template <class Graph>
struct StergiouAlgorithm {
  Graph& GA;
  size_t n;
  StergiouAlgorithm(Graph& GA, size_t n) :
    GA(GA), n(n) { }

  void initialize(pbbs::sequence<parent>& P) { }

  template <SamplingOption sampling_option>
  void compute_components(pbbs::sequence<parent>& parents, uintE frequent_comp = UINT_E_MAX) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    pbbs::sequence<parent> previous_parents(n);

    auto parents_changed = true;
    size_t rounds = 0;
    while (parents_changed) {
      parents_changed = false;
      if (rounds > 0) {
        parallel_for(0, n, [&] (size_t i) {
          previous_parents[i] = parents[i];
        });
      }

      parallel_for(0, n, [&] (size_t i) {
        auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
          uintE parent_u = previous_parents[u];
          uintE parent_v = previous_parents[v];
          bool updated = false;
          if (parents[v] > parent_u) {
            updated |= pbbs::write_min(&parents[v], parent_u, std::less<uintE>());
          }
          if (parents[u] > parent_v) {
            updated |= pbbs::write_min(&parents[u], parent_v, std::less<uintE>());
          }
          if (updated && !parents_changed) {
            parents_changed = true;
          }
        };
        if constexpr (sampling_option != no_sampling) {
          if (parents[i] != frequent_comp) {
            GA.get_vertex(i).mapOutNgh(i, map_f);
          }
        } else {
          GA.get_vertex(i).mapOutNgh(i, map_f);
        }
      });

      // Shortcut
      parallel_for(0, n, [&] (size_t u) {
          lt::primitives::shortcut(u, parents);
      });
    }
  }

  /* processes a batch of edge insertions and edge queries */
  template <class Seq>
  void process_batch(pbbs::sequence<parent>& parents, Seq& batch, size_t insert_to_query) {

//    auto parents_changed = true;
//    auto& P = parents;
//    size_t round = 0;
//    while (parents_changed) {
//      parents_changed = false;
//
//
//      parallel_for(0, batch.size(), [&] (size_t i) {
//        auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
//          uintE u, v;
//          std::tie(u,v) = batch[i];
//
//          uintE parent_u = previous_parent[u];
//          uintE parent_v = previous_parent[v];
//          bool updated = false;
//          if (parents[v] > parent_u) {
//            updated |= pbbs::write_min(&parents[v], parent_u);
//          }
//          if (parent[u] > parent_v) {
//            updated |= pbbs::write_min(&parents[u], parent_v);
//          }
//          if (updated && !parents_changed) {
//            parents_changed = true;
//          }
//        };
//        if constexpr (sampling_option != no_sampling) {
//          if (parents[i] != frequent_comp) {
//            GA.get_vertex(i).mapOutNgh(i, map_f);
//          }
//        } else {
//          GA.get_vertex(i).mapOutNgh(i, map_f);
//        }
//      });
//
//      // Shortcut
//      parallel_for(0, batch.size(), [&] (size_t i) {
//        uintE u, v;
//        std::tie(u,v) = batch[i];
//        if (i % insert_to_query != 0) { /* is an update */
//          shortcut(u, P);
//        } else { /* answer the query, if it is the first round*/
//          if (round == 0) {
//            while (P[u] != P[P[u]]) {
//              P[u] = P[P[u]];
//            }
//            /* found p_u */
//
//            while (P[v] != P[P[v]]) {
//              P[v] = P[P[v]];
//            }
//            /* found p_v */
//          }
//        }
//      });
//
//      round++;
//    }
//
//    // Process queries, no need to skip.
//    parallel_for(0, batch.size(), [&] (size_t i) {
//      uintE u, v;
//      std::tie(u,v) = batch[i];
//      if (i % insert_to_query != 0) { /* is an update */
//        shortcut(u, P);
//      } else { /* answer the query, if it is the first round*/
//        if (round == 0) {
//          while (P[u] != P[P[u]]) {
//            P[u] = P[P[u]];
//          }
//          /* found p_u */
//
//          while (P[v] != P[P[v]]) {
//            P[v] = P[P[v]];
//          }
//          /* found p_v */
//        }
//      }
//    });
//
  }

};



}  // namespace lt
