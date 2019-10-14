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
#include "benchmarks/Connectivity/Common/common.h"

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

enum LiuTarjanConnectOption {
  simple_connect, parent_connect, extended_connect
};
enum LiuTarjanUpdateOption {
  simple_update, root_update
};
enum LiuTarjanShortcutOption {
  shortcut, full_shortcut
};
enum LiuTarjanAlterOption {
  alter, no_alter
};

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
  Connect& connect;
  Update& update;
  Shortcut& shortcut;
  LiuTarjanAlgorithm(Graph& GA, Connect& connect, Update& update, Shortcut& shortcut) :
    GA(GA), connect(connect), update(update), shortcut(shortcut) {}

  template <bool provides_frequent_comp>
  void compute_components(pbbs::sequence<parent>& P, uintE frequent_comp = UINT_E_MAX) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    auto parents_changed = true;
    while (parents_changed) {
      parents_changed = false;

      // Parent-Connect
      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        bool updated = connect(u, v, P);
        if (updated && !parents_changed) {
          parents_changed = true;
        }
      };
      GA.map_edges(map_f);

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
};

template <class Connect,
          LiuTarjanConnectOption connect_option,
          class Update,
          LiuTarjanUpdateOption update_option,
          class Shortcut,
          LiuTarjanShortcutOption shortcut_option,
          class Alter,
          LiuTarjanAlterOption alter_option,
          class Graph>
struct LiuTarjanAlgorithmAlter {
  Graph& GA;
  Connect& connect;
  Update& update;
  Shortcut& shortcut;
  Alter& alter;
  LiuTarjanAlgorithmAlter(Graph& GA, Connect& connect, Update& update, Shortcut& shortcut, Alter& alter) :
    GA(GA), connect(connect), update(update), shortcut(shortcut), alter(alter) {}

  template <bool provides_frequent_comp>
  void compute_components(pbbs::sequence<parent>& P, uintE frequent_comp = UINT_E_MAX) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    auto parents_changed = true;
    while (parents_changed) {
      parents_changed = false;

      // Parent-Connect
      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        if (u != v) {
          bool updated = connect(u, v, P);
          if (updated && !parents_changed) {
            parents_changed = true;
          }
        }
      };
      GA.map_edges(map_f);

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

      auto alter_f = [&] (const uintE& u, const uintE& v, const W& wgh) -> std::tuple<uintE, uintE, W> {
        if (u != v) {
          uintE au, av;
          std::tie(au, av) = alter(u, v, P);
          return std::make_tuple(au, av, wgh);
        }
        return std::make_tuple(u, v, wgh);
      };
      GA.alter_edges(alter_f);
    }
  }
};


}  // namespace lt
