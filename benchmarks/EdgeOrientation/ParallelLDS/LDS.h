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
#include <stack>

#include "gbbs/gbbs.h"

#include "sparse_set.h"

namespace gbbs {

struct LDS {

  static constexpr double kDelta = 9.0;
  static constexpr double kUpperConstant = 2 + (static_cast<double>(3) / kDelta);
  static constexpr double kEpsilon = 3.0;
  static constexpr double kOnePlusEpsilon = 1 + kEpsilon;

  static constexpr uintE kUpLevel = UINT_E_MAX;

  using levelset = gbbs::sparse_set<uintE>;
  using down_neighbors = parlay::sequence<levelset>;
  using up_neighbors = levelset;
  using edge_type = std::pair<uintE, uintE>;

  struct LDSVertex {
    uintE level;  // level of this vertex
    down_neighbors
        down;         // neighbors in levels < level, bucketed by their level.
    up_neighbors up;  // neighbors

    LDSVertex() : level(0) {}

    inline double group_degree(size_t group) const {
      return pow(kOnePlusEpsilon, group);
    }

    inline bool upper_invariant(const size_t levels_per_group) const {
      uintE group = level / levels_per_group;
      uintE up_degree = kUpperConstant * group_degree(group);
      // std::cout << "up_size = " << up.num_elms() << " up_degree = " << up_degree << std::endl;
      return up.num_elms() <= up_degree;
    }

    inline bool lower_invariant(const size_t levels_per_group) const {
      if (level == 0) return true;
      uintE lower_group = (level - 1) / levels_per_group;
      auto up_size = up.num_elms();
      auto prev_level_size = down[level - 1].num_elms();
      size_t our_group_degree = static_cast<size_t>(group_degree(lower_group));
      return (up_size + prev_level_size) >= our_group_degree;
    }

    inline bool is_dirty(const size_t levels_per_group) const {
      bool upper = upper_invariant(levels_per_group);
      bool lower = lower_invariant(levels_per_group);
      return !(upper && lower);
    }

  };

  size_t n;  // number of vertices
  size_t levels_per_group;  // number of inner-levels per group,  O(\log n) many.
  parlay::sequence<LDSVertex> L;

  LDS(size_t n) : n(n) {
    levels_per_group = ceil(log(n) / log(kOnePlusEpsilon));
    L = parlay::sequence<LDSVertex>(n);
  }

  uintE get_level(uintE ngh) {
    return L[ngh].level;
  }

  bool edge_exists(edge_type e) {
    auto[u, v] = e;
    auto l_u = L[u].level;
    auto l_v = L[v].level;
    if (l_u < l_v) {  // look in up(u)
      return L[u].up.contains(v);
    } else {  // look in up(v)
      return L[v].up.contains(u);
    }
  }


  // Input: sequence of (level, vertex_id) pairs.
  template <class Seq, class Buckets>
  void update_buckets(Seq&& possibly_dirty, Buckets& buckets) {
    using level_and_vtx = std::pair<uintE, uintE>;
    // // Remap by level, sending vertices not wanting to move to UINT_E_MAX.
    // parallel_for(0, possibly_dirty.size(), [&] (size_t i) {
    //   auto v = possibly_dirty[i];
    //   if (!L[v].is_dirty(levels_per_group)) possibly_dirty[i] = UINT_E_MAX;
    // });

    // Compute dirty vertices, which have either lower / upper threshold
    // violated.
    auto level_and_vtx_seq = parlay::delayed_seq<level_and_vtx>(possibly_dirty.size(), [&] (size_t i) {
      uintE v = possibly_dirty[i];
      uintE level = UINT_E_MAX;
      if (L[v].is_dirty(levels_per_group)) level = L[v].level;
      return std::make_pair(level, v);
    });
    auto dirty = parlay::filter(level_and_vtx_seq, [&] (const level_and_vtx& lv) {
      return lv.first != UINT_E_MAX;
    });
    std::cout << "Num dirty elms = " << dirty.size() << std::endl;
    if (dirty.size() == 0) return;

    // Sort by level.
    auto get_key = [&] (const level_and_vtx& elm) { return elm.first; };
    parlay::integer_sort_inplace(parlay::make_slice(dirty), get_key);

    // Get unique levels.
    auto bool_seq = parlay::delayed_seq<bool>(dirty.size() + 1, [&] (size_t i) {
      return (i == 0) || (i == dirty.size()) || (dirty[i].first != dirty[i-1].first);
    });
    auto level_starts = parlay::pack_index(bool_seq);
    assert(level_starts[level_starts.size() - 1] == dirty.size());
    assert(level_starts.size() >= 2);

    std::cout << "Number of affected levels is: " << (level_starts.size()-1) << std::endl;

    // using dirty_elts = sequence<uintE>;
    uintE max_current_level = dirty[level_starts[level_starts.size() - 1]].first + 1;
    std::cout << "max_current_level = " << max_current_level << std::endl;
    if (buckets.size() < max_current_level) {
      buckets.resize(max_current_level);
    }

    parallel_for(0, level_starts.size() - 1, [&] (size_t i) {
      uintE idx = level_starts[i];
      uintE level = dirty[idx].first;
      uintE num_in_level = level_starts[i+1] - idx;

      auto stuff_in_level = parlay::delayed_seq<uintE>(num_in_level, [&] (size_t j) {
        return dirty[idx + j].second;
      });
      buckets[level].append(stuff_in_level);
    });
  }

  // returns the total number of moved vertices
  template <class Buckets>
  size_t rebalance_insertions(Buckets&& buckets, size_t cur_bucket_id, size_t total_moved = 0) {
    if (cur_bucket_id >= buckets.size())
      return total_moved;

    auto& cur_bucket = buckets[cur_bucket_id];
    if (cur_bucket.size() == 0)
      return rebalance_insertions(buckets, cur_bucket_id+1, total_moved);

    // 1. figure out the target level for each vertex in cur_bucket.
    using level_and_vtx = std::pair<uintE, uintE>;

//    auto level_and_vtx_seq = parlay::delayed_seq<level_and_vtx>(possibly_dirty.size(), [&] (size_t i) {
//      uintE v = possibly_dirty[i];
//      uintE level = UINT_E_MAX;
//      if (L[v].is_dirty(levels_per_group)) level = L[v].level;
//      return std::make_pair(level, v);
//    });
//    auto dirty = parlay::filter(level_and_vtx_seq, [&] (const level_and_vtx& lv) {
//      return lv.first != UINT_E_MAX;
//    });



  }


  template <class Seq>
  void batch_insertion(const Seq& insertions_unfiltered) {
    // Remove edges that already exist from the input.
    auto insertions_filtered = parlay::filter(parlay::make_slice(insertions_unfiltered),
        [&] (const edge_type& e) { return !edge_exists(e); });

    // Duplicate the edges in both directions and sort.
    auto insertions_dup = sequence<edge_type>::uninitialized(2*insertions_filtered.size());
    parallel_for(0, insertions_filtered.size(), [&] (size_t i) {
      auto [u, v] = insertions_filtered[i];
      insertions_dup[2*i] = {u, v};
      insertions_dup[2*i + 1] = {v, u};
    });
    auto compare_tup = [&] (const edge_type& l, const edge_type& r) { return l < r; };
    parlay::sort_inplace(parlay::make_slice(insertions_dup), compare_tup);

    // Remove duplicate edges to get the insertions.
    auto not_dup_seq = parlay::delayed_seq<bool>(insertions_dup.size(), [&] (size_t i) {
      auto [u, v] = insertions_dup[i];
      bool not_self_loop = u != v;
      return not_self_loop && ((i == 0) || (insertions_dup[i] != insertions_dup[i-1]));
    });
    auto insertions = parlay::pack(parlay::make_slice(insertions_dup), not_dup_seq);


    // Compute the starts of each (modified) vertex's new edges.
    auto bool_seq = parlay::delayed_seq<bool>(insertions.size() + 1, [&] (size_t i) {
      return (i == 0) || (i == insertions.size()) || (std::get<0>(insertions[i-1]) != std::get<0>(insertions[i]));
    });
    auto starts = parlay::pack_index(bool_seq);

    // Save the vertex ids and starts (hypersparse CSR format). The next step
    // will overwrite the edge pairs to store (neighbor, current_level).
    // (saving is not necessary if we modify + sort in a single parallel loop)
    auto affected = sequence<uintE>::from_function(starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE vtx_id = std::get<0>(insertions[idx]);
      return vtx_id;
    });

    // Map over the vertices being modified. Overwrite their incident edges with
    // (level_id, neighbor_id) pairs, sort by level_id, resize each level to the
    // correct size, and then insert in parallel.
    parallel_for(0, starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE vtx = std::get<0>(insertions[idx]);
      uintE our_level = L[vtx].level;
      uintE incoming_degree = starts[i+1] - starts[i];
      auto neighbors = parlay::make_slice(insertions.begin() + idx,
          insertions.begin() + idx + incoming_degree);

      // std::cout << "Processing vtx = " << vtx << " idx = " << idx << " incoming_degree = " << incoming_degree << std::endl;

      // Map the incident edges to (level, neighbor_id).
      parallel_for(0, incoming_degree, [&] (size_t off) {
        auto [u, v] = neighbors[off];
        assert(vtx == u);
        uintE neighbor_level = L[v].level;
        if (neighbor_level >= our_level) { neighbor_level = kUpLevel; }
        neighbors[off] = {neighbor_level, v};
      });

      // Sort neighbors by level.
      parlay::sort_inplace(neighbors);

      // Resize level containers. Can do this in parallel in many ways (e.g.,
      // pack), but a simple way that avoids memory allocation is to map over
      // everyone, and have the first index for each level search for the level.
      // If we use a linear search the algorithm is work eff. and has depth
      // proportional to the max incoming size of a level. We can also improve
      // to log(n) depth by using a doubling search.
      parallel_for(0, neighbors.size(), [&] (size_t i) {
        if ((i == 0) || neighbors[i].first != neighbors[i-1].first) {  // start of a new level
          uintE level_id = neighbors[i].first;
          size_t j = i;
          for (; j<neighbors.size(); j++) {
            if (neighbors[j].first != level_id) break;
          }
          uintE num_in_level = j - i;

          if (level_id != kUpLevel) {
            L[vtx].down[level_id].resize(num_in_level);
          } else {
            L[vtx].up.resize(num_in_level);
          }
        }
      });

      // std::cout << "Starting insertions into hash tables" << std::endl;
      // Insert neighbors into the correct level incident to us.
      parallel_for(0, neighbors.size(), [&] (size_t i) {
        auto [level_id, v] = neighbors[i];
        if (level_id != kUpLevel) {
          bool inserted = L[vtx].down[level_id].insert(v);
          assert(inserted);
        } else {
          bool inserted = L[vtx].up.insert(v);
          assert(inserted);
        }
      });
      // std::cout << "Finished insertions into hash tables" << std::endl;

    }, 10000000000);  // for testing

    // New edges are done being inserted. Update the level structure.
    // Interface: supply vertex seq -> process will settle everything.

    using dirty_elts = sequence<uintE>;
    sequence<dirty_elts> buckets;

    // Place the affected vertices into buckets based on their current level.
    update_buckets(std::move(affected), buckets);

    // Update the bucket structure.
    size_t total_moved = rebalance_insertions(std::move(buckets), 0);

    std::cout << "During insertions, " << total_moved << " vertices moved." << std::endl;
  }

  void check_invariants() {
    bool invs_ok = true;
    for (size_t i=0; i<n; i++) {
      bool upper_ok = L[i].upper_invariant(levels_per_group);
      bool lower_ok = L[i].lower_invariant(levels_per_group);
      assert(upper_ok);
      assert(lower_ok);
      invs_ok &= upper_ok;
      invs_ok &= lower_ok;
    }
    std::cout << "invs ok is: " << invs_ok << std::endl;
  }

  inline uintE group_for_level(uintE level) const {
    return level / levels_per_group;
  }
};

template <class Graph>
inline void RunLDS(Graph& G) {
  // using W = typename Graph::weight_type;
  size_t n = G.n;
  auto layers = LDS(n);

  auto insertions = sequence<LDS::edge_type>::uninitialized(100);
  for (size_t i=0; i<10; i++) {
    size_t off = 10*i;
    for (size_t j=0; j<10; j++) {
      insertions[off + j] = {i, i + j + 1};
      // std::cout << "inserting : " << i << " " << (i + j + 1) << std::endl;
    }
  }
  layers.batch_insertion(insertions);


//  for (size_t i = 0; i < n; i++) {
//    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
//      if (u < v) {
//        layers.insert_edge({u, v});
//      }
//    };
//    G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
//  }

  std::cout << "Finished all insertions!" << std::endl;
  layers.check_invariants();
  std::cout << "Finished check" << std::endl;

}


//  // Moving u from level to level - 1.
//  template <class Levels>
//  void level_decrease(uintE u, Levels& L) {
//    uintE level = L[u].level;
//    auto& up = L[u].up;
//    assert(level > 0);
//    auto& prev_level = L[u].down[level - 1];
//
//    for (const auto& ngh : prev_level) {
//      up.insert(ngh);
//    }
//    L[u].down.pop_back();  // delete the last level in u's structure.
//
//    for (const auto& ngh : up) {
//      if (get_level(ngh) == level) {
//        L[ngh].up.erase(L[ngh].up.find(u));
//        L[ngh].down[level-1].insert(u);
//
//      } else if (get_level(ngh) > level) {
//        L[ngh].down[level].erase(L[ngh].down[level].find(u));
//        L[ngh].down[level - 1].insert(u);
//
//        if (get_level(ngh) == level + 1) {
//          Dirty.push(ngh);
//        }
//      } else {
//        // u is still "up" for this ngh, no need to update.
//        assert(get_level(ngh) == (level - 1));
//      }
//    }
//    L[u].level--;  // decrease level
//  }
//
//  // Moving u from level to level + 1.
//  template <class Levels>
//  void level_increase(uintE u, Levels& L) {
//    uintE level = L[u].level;
//    std::vector<uintE> same_level;
//    auto& up = L[u].up;
//
//    for (auto it = up.begin(); it != up.end();) {
//      uintE ngh = *it;
//      if (L[ngh].level == level) {
//        same_level.emplace_back(ngh);
//        it = up.erase(it);
//        // u is still "up" for this ngh, no need to update.
//      } else {
//        it++;
//        // Must update ngh's accounting of u.
//        if (L[ngh].level > level + 1) {
//          L[ngh].down[level].erase(L[ngh].down[level].find(u));
//          L[ngh].down[level + 1].insert(u);
//        } else {
//          assert(L[ngh].level == level + 1);
//          L[ngh].down[level].erase(L[ngh].down[level].find(u));
//          L[ngh].up.insert(u);
//
//          Dirty.push(ngh);
//        }
//      }
//    }
//    // We've now split L[u].up into stuff in the same level (before the
//    // update) and stuff in levels >= level + 1. Insert same_level elms
//    // into down.
//    auto& down = L[u].down;
//    down.emplace_back(std::unordered_set<uintE>());
//    assert(down.size() == level + 1);  // [0, level)
//    for (const auto& ngh : same_level) {
//      down[level].insert(ngh);
//    }
//    L[u].level++;  // Increase level.
//  }
//
//  void fixup() {
//    while (!Dirty.empty()) {
//      uintE u = Dirty.top();
//      Dirty.pop();
//      if (!L[u].upper_invariant(levels_per_group)) {
//        // Move u to level i+1.
//        level_increase(u, L);
//        Dirty.push(u);  // u might need to move up more levels.
//        // std::cout << "(move up) pushing u = " << u << std::endl;
//      } else if (!L[u].lower_invariant(levels_per_group)) {
//        level_decrease(u, L);
//        Dirty.push(u);  // u might need to move down more levels.
//        // std::cout << "(move down) pushing u = " << u << std::endl;
//      }
//    }
//  }



}  // namespace gbbs
