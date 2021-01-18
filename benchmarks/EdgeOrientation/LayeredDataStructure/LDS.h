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

#include <unordered_set>
#include <vector>
#include <stack>

#include "gbbs/gbbs.h"

namespace gbbs {

// Some notation and data structure description, basically copied from
// Henzinger et al. (2020). (https://arxiv.org/pdf/2002.10142.pdf)
//
// notation:
// levels : log^2(n) many, notated as l. level of a vertex is l(v)
// groups : log(n) many. group of a particular level if g(l)
//
// invariants:
// (1) each v \in V has at most 5 * 2^{g(l(v))} neighbors at its own, or higher
//     levels. (not too many neighbors above)
// (2) each v \in V where l(v) > 1 has at least 2^{g(l(v) - 1)} neighbors in
//     level l(v) - 1. (enough neighbors in the level below)
// vertices that do not satisfy these invariants are said to be _dirty_.
//
// The Henzinger et al. algorithm works by moving vertices that are dirty to
// a level where both invariants are satisfied.
//
// data structures maintained:
// For each vertex v, and each level l' < l(v), maintain the following:
// - doubly-linked list neighbors(v, l') containing all neighbors of V in
//   level l'
//
// Furthermore, we have another doubly-linked list up_neighbors(v)
// containing all neighbors v in levels l' >= l(v).
//
// For each edge (v,u), we store a pointer to the position of u in neighbors
// or up_neighbors (and vice versa).
//
// In addition, we store the size of each list (neighbors and up_neighbors)
// which lets us test invariants (1) and (2) in O(1) work.
struct LDS {
  size_t n;  // number of vertices

  using level = std::unordered_set<uintE>;
  using down_neighbors = std::vector<level>;
  using up_neighbors = level;
  using edge_type = std::pair<uintE, uintE>;

  struct LDSVertex {
    uintE level;  // level of this vertex
    down_neighbors
        down;         // neighbors in levels < level, bucketed by their level.
    up_neighbors up;  // neighbors

    LDSVertex() : level(0) {}

    void insert_neighbor(uintE v, uintE l_v) {
      if (l_v < level) {
        assert(down.size() > l_v);
        down[l_v].insert(v);
      } else {
        up.insert(v);
      }
    }

    void remove_neighbor(uintE v, uintE l_v) {
      if (l_v < level) {
        down[l_v].erase(down[l_v].find(v));
      } else {
        up.erase(up.find(v));
      }
    }

    inline bool upper_invariant(const size_t levels_per_group) const {
      uintE group = level / levels_per_group;
//      std::cout << "up.size = " << up.size() << " thresh = " << static_cast<size_t>(5 * (1 << group)) << std::endl;
      return up.size() <= static_cast<size_t>(5 * (1 << group));
    }

    inline bool lower_invariant(const size_t levels_per_group) const {
      if (level == 0) return true;
      uintE lower_group = (level - 1) / levels_per_group;
      return down[level - 1].size() >= static_cast<size_t>(1 << lower_group);
    }

  };

  // number of inner-levels per group,  O(\log n) many.
  size_t levels_per_group;
  pbbs::sequence<LDSVertex> L;
  std::stack<uintE> Dirty;

  LDS(size_t n) : n(n) {
    levels_per_group = pbbs::log2_up(n);
    L = pbbs::sequence<LDSVertex>(n);
  }

  uintE get_level(uintE ngh) {
    return L[ngh].level;
  }

  // Moving u from level to level - 1.
  template <class Levels>
  void level_decrease(uintE u, Levels& L) {
    uintE level = L[u].level;
    auto& up = L[u].up;
    assert(level > 0);
    auto& prev_level = L[u].down[level - 1];

    for (const auto& ngh : prev_level) {
      up.insert(ngh);
    }
    L[u].down.pop_back();  // delete the last level in u's structure.

    for (const auto& ngh : up) {
      if (get_level(ngh) == level) {
        L[ngh].up.erase(L[ngh].up.find(u));
        L[ngh].down[level-1].insert(u);
      } else if (get_level(ngh) > level) {
        L[ngh].down[level].erase(L[ngh].down[level].find(u));
        L[ngh].down[level - 1].insert(u);

        if (get_level(ngh) == level + 1) {
          Dirty.push(ngh);
        }
      } else {
        // u is still "up" for this ngh, no need to update.
        assert(get_level(ngh) == (level - 1));
      }
    }
    L[u].level--;  // decrease level
  }

  // Moving u from level to level + 1.
  template <class Levels>
  void level_increase(uintE u, Levels& L) {
    uintE level = L[u].level;
    std::vector<uintE> same_level;
    auto& up = L[u].up;

    for (auto it = up.begin(); it != up.end();) {
      uintE ngh = *it;
      if (L[ngh].level == level) {
        same_level.emplace_back(ngh);
        it = up.erase(it);
        // u is still "up" for this ngh, no need to update.
      } else {
        it++;
        // Must update ngh's accounting of u.
        if (L[ngh].level > level + 1) {
          L[ngh].down[level].erase(L[ngh].down[level].find(u));
          L[ngh].down[level + 1].insert(u);
        } else {
          assert(L[ngh].level == level + 1);
          L[ngh].down[level].erase(L[ngh].down[level].find(u));
          L[ngh].up.insert(u);

          Dirty.push(ngh);
        }
      }
    }
    // We've now split L[u].up into stuff in the same level (before the
    // update) and stuff in levels >= level + 1. Insert same_level elms
    // into down.
    auto& down = L[u].down;
    down.emplace_back(std::unordered_set<uintE>());
    assert(down.size() == level + 1);  // [0, level)
    for (const auto& ngh : same_level) {
      down[level].insert(ngh);
    }
    L[u].level++;  // Increase level.
  }

  void fixup() {
    while (!Dirty.empty()) {
      uintE u = Dirty.top();
      Dirty.pop();
      if (!L[u].upper_invariant(levels_per_group)) {
        // Move u to level i+1.
        level_increase(u, L);
        Dirty.push(u);  // u might need to move up more levels.
      } else if (!L[u].lower_invariant(levels_per_group)) {
        level_decrease(u, L);
      }
    }
  }

  bool edge_exists(edge_type e) {
    auto[u, v] = e;
    auto l_u = L[u].level;
    auto l_v = L[v].level;
    if (l_u < l_v) {  // look in up(u)
      return (L[u].up.find(v) != L[u].up.end());
    } else {  // look in up(v)
      return (L[v].up.find(u) != L[v].up.end());
    }
  }

  bool insert_edge(edge_type e) {
    if (edge_exists(e)) return false;
    auto[u, v] = e;
    std::cout << "inserting edge: " << u << " " << v << std::endl;
    auto l_u = L[u].level;
    auto l_v = L[v].level;
    L[u].insert_neighbor(v, l_v);
    L[v].insert_neighbor(u, l_u);

    Dirty.push(u);
    Dirty.push(v);
    fixup();
    return true;
  }

  bool delete_edge(edge_type e) {
    assert(edge_exists(e));
    if (!edge_exists(e)) return false;
    auto[u, v] = e;
    auto l_u = L[u].level;
    auto l_v = L[v].level;
    L[u].remove_neighbor(v, l_v);
    L[v].remove_neighbor(u, l_u);

    Dirty.push(u); Dirty.push(v);
    fixup();
    return true;
  }

  inline uintE group_for_level(uintE level) const {
    return level / levels_per_group;
  }
};

template <class Graph>
inline void RunLDS(Graph& G) {
  using W = typename Graph::weight_type;
  auto layers = LDS(G.n);
  for (size_t i = 0; i < 1000; i++) {
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (u < v) {
        // std::cout << "inserting u = " << u << " v = " << v << std::endl;
        bool ret = layers.insert_edge({u, v});
        assert(ret);
      }
    };
    G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
  }

  std::cout << "Finished all insertions!" << std::endl;

  for (size_t i = 0; i < 1000; i++) {
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (u < v) {
        // std::cout << "deleting u = " << u << " v = " << v << std::endl;
        bool ret = layers.delete_edge({u, v});
        assert(ret);
      }
    };
    G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
  }

  std::cout << "Finished all deletions!" << std::endl;

  size_t sum_lev = 0;
  for (size_t i=0; i<G.n; i++) {
    sum_lev += layers.L[i].level;
  }
  std::cout << "sum_lev = " << sum_lev << std::endl;
}

}  // namespace gbbs
