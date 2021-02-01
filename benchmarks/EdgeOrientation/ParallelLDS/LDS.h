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
#include "gbbs/dynamic_graph_io.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

#include "sparse_set.h"

namespace gbbs {

struct LDS {

  double delta = 48.0;
  double UpperConstant = 2 + ((double) 3 / delta);
  double eps = 1.6;
  double OnePlusEps = 1 + eps;

  static constexpr uintE kUpLevel = UINT_E_MAX;
  // Used to indicate that this vertex is not moving.
  static constexpr uintE kNotMoving = UINT_E_MAX;

  using levelset = gbbs::sparse_set<uintE>;
  using down_neighbors = parlay::sequence<levelset>;
  using up_neighbors = levelset;
  using edge_type = std::pair<uintE, uintE>;

  struct LDSVertex {
    uintE level;         // The level of this vertex.
    uintE desire_level;  // The desire level of this vertex (set only when moving this vertex).
    down_neighbors down; // The neighbors in levels < level, bucketed by their level.
    up_neighbors up;     // The neighbors in levels >= level.

    LDSVertex() : level(0), desire_level(kNotMoving) {}

    // Used when Invariant 1 (upper invariant) is violated.
    template <class Levels>
    inline uintE get_desire_level_upwards(uintE vtx_id, Levels& L, const size_t levels_per_group,
            double upper_constant, double eps) const {
      assert(!upper_invariant(levels_per_group, upper_constant, eps));
      using LV = std::pair<uintE, uintE>;
      auto LV_seq = parlay::delayed_seq<LV>(up.size(), [&] (size_t i) {
        uintE v = up.table[i];
        uintE l_v = UINT_E_MAX;
        if (levelset::valid(v)) l_v = L[v].level;
        return std::make_pair(l_v, v);
      });

      auto LVs = parlay::filter(LV_seq, [&] (const LV& lv) {
        return lv.first != UINT_E_MAX;
      });

      parlay::sort_inplace(parlay::make_slice(LVs));

      uintE new_level = UINT_E_MAX;
      uintE prev_level = level;
      for (size_t i=0; i<LVs.size(); i++) {
        // Start of a new level among our up-neighbors
        if ((i == 0) || (LVs[i].first != LVs[i-1].first)) {
          uintE cur_level = LVs[i].first;
          if (cur_level == level) continue;  // we are surely not moving to the same level

          size_t degree_at_level = LVs.size() - i;
          bool done = false;
          // Consider all levels from prev_level to the cur_level.
          for (uintE j = prev_level+1; j <= cur_level; j++) {
            uintE group = j / levels_per_group;
            uintE up_degree = upper_constant*group_degree(group, eps);

            if (degree_at_level <= up_degree) {
              new_level = j;
              done = true;
              break;
            }
          }
          if (done) break;

          prev_level = cur_level;
        }
      }

      // Could new_level still be UINT_E_MAX? possibly, e.g., if all up
      // neighbors are on the same level, and there are too many on this level.
      // In this case, the desired level is one larger than the max level.
      if (new_level == UINT_E_MAX)
        new_level = prev_level + 1;


      return new_level;

//      uintE new_level = UINT_E_MAX;
//      uintE max_level = 0;
//      // todo: make search parallel
//      for (size_t i=0; i<LVs.size(); i++) {
//        if ((i == 0) || (LVs[i].first != LVs[i-1].first)) {
//          uintE level = LVs[i].first;
//          size_t degree_at_level = LVs.size() - i;
//
//          uintE group = level / levels_per_group;
//          uintE up_degree = kUpperConstant * group_degree(group);
//          max_level = level;
//          if (degree_at_level <= up_degree) {
//            new_level = level;
//            break;
//          }
//        }
//      }
//
//      // Could new_level still be UINT_E_MAX? possibly, e.g., if all up
//      // neighbors are on the same level, and there are too many on this level.
//      // In this case, the desired level is one larger than the max level.

//
//      if (new_level == UINT_E_MAX)
//        new_level = max_level + 1;
//
//      return new_level;
    }

    // Used when Invariant 2 (up* degree invariant) is violated.
    template <class Levels>
    inline uintE get_desire_level_downwards(uintE vtx_id, Levels& L, const size_t levels_per_group,
            double upper_constant, double eps) const {
        assert(!lower_invariant(levels_per_group, eps));
        assert(level > 0);

        // Currently sequential going level by level but we probably want to
        // optimize later.
        //
        // This is the current level of the node.
        uintE cur_level = level - 1;
        if (cur_level == 0) return cur_level;
        uintE num_up_neighbors = num_neighbors_higher_than_level((uintE) cur_level);
        uintE num_up_star_neighbors = num_up_neighbors;
        while (cur_level > 0) {
            num_up_neighbors = num_up_star_neighbors;
            num_up_star_neighbors += down[cur_level-1].num_elms();
            if (upstar_degree_satisfies_invariant(cur_level, levels_per_group, num_up_star_neighbors,
                        num_up_neighbors, upper_constant, eps)) {
                break;
            }
            cur_level -= 1;
        }

        return cur_level;
    }

    inline uintE num_up_neighbors() const {
      return up.num_elms();
    }

    inline uintE num_up_star_neighbors() const {
      if (level == 0) return up.num_elms();
      return up.num_elms() + down[level - 1].num_elms();
    }

    // Get the sum of the number of neighbors at or higher than start_level
    //
    // TODO: this can be optimized by making the summation parallel, currently
    // it is sequential.
    inline uintE num_neighbors_higher_than_level(uintE start_level) const {
        uintE num_flipped_neighbors = 0;
        while (start_level < level) {
            num_flipped_neighbors += down[start_level].num_elms();
            start_level++;
        }
        return num_flipped_neighbors + up.num_elms();
    }

    template <class OutputSeq>
    inline uintE emit_up_neighbors(OutputSeq output, uintE our_id) const {
      // TODO: Could filter using filter, but let's do this seq. for now...
      size_t off = 0;
      for (size_t i=0; i<up.table_seq.size(); i++) {
        uintE v = up.table[i];
        if (levelset::valid(v))
          output[off++] = std::make_pair(v, our_id);
      }
      return off;
    }

    template <class OutputSeq>
    inline uintE emit_neighbors_higher_than_level(OutputSeq output, uintE our_id, uintE start_level) const {
      // TODO: Could filter using filter, but let's do this seq. for now...
      size_t off = 0;
      for (size_t j = start_level; j < level; j++) {
        for (size_t i=0; i < down[j].table_seq.size(); i++) {
            uintE v = down[j].table[i];
            if (levelset::valid(v))
                output[off++] = std::make_pair(v, our_id);
        }
      }

      for (size_t i = 0; i < up.table_seq.size(); i++) {
        uintE v = up.table[i];
        if (levelset::valid(v))
            output[off++] = std::make_pair(v, our_id);
      }
      return off;
    }


    template <class Levels>
    inline void filter_up_neighbors(uintE vtx_id, Levels& L) {
      uintE removed = 0;
      auto all_up = up.entries();
      assert(desire_level != kNotMoving);

      auto resize_sizes_seq = parlay::sequence<uintE>(desire_level - level, (uintE)0);
      auto resize_sizes = resize_sizes_seq.begin();

      for (size_t i=0; i<up.table_seq.size(); i++) {
        uintE v = up.table[i];
        if (levelset::valid(v)) {
          if (L[v].level == level && L[v].desire_level == kNotMoving) {
            up.table[i] = levelset::kTombstone;
            uintE norm_v = L[v].level - level;
            resize_sizes[norm_v]++;
            removed++;
          } else if (L[v].level > level && L[v].level < desire_level) {
            // Another case: L[v].level < desire_level
            up.table[i] = levelset::kTombstone;
            uintE norm_v = L[v].level - level;
            resize_sizes[norm_v]++;
            removed++;
          }
        }
      }

      // Perform resizing.
      up.resize_down(removed);
      assert(down.size() > level);
      for (size_t i=0; i<resize_sizes_seq.size(); i++) {
        auto l = level + i;
        down[l].resize(resize_sizes[i]);
      }

      uintE inserted = 0;
      for (size_t i=0; i<all_up.size(); i++) {
        uintE v = all_up[i];
        if (L[v].level == level && L[v].desire_level == kNotMoving) {
          down[level].insert(v);
          inserted++;
        } else if (L[v].level > level && L[v].level < desire_level) {
          down[L[v].level].insert(v);
          inserted++;
        }
      }
      assert(removed == inserted);
    }

    inline double group_degree(size_t group, double eps) const {
      return pow((1 + eps), group);
    }

    inline bool upper_invariant(const size_t levels_per_group, double upper_constant,
             double eps) const {
      uintE group = level / levels_per_group;
      uintE up_degree = upper_constant * group_degree(group, eps);
      return up.num_elms() <= up_degree;
    }

    inline bool lower_invariant(const size_t levels_per_group, double eps) const {
      if (level == 0) return true;
      uintE lower_group = (level - 1) / levels_per_group;
      //auto up_size = up.num_elms();
      //auto prev_level_size = down[level - 1].num_elms();
      size_t our_group_degree = static_cast<size_t>(group_degree(lower_group, eps));
      return num_up_star_neighbors() >= our_group_degree;
    }

    // Method for computing whether the updegree from cur_level is enough to
    // satisfy Invariant 2.
    inline bool upstar_degree_satisfies_invariant(uintE cur_level, const size_t
            levels_per_group, uintE num_up_star_neighbors, uintE num_up_neighbors,
            double upper_constant, double eps) const {
        if (cur_level == 0) return true;
        uintE group = (cur_level - 1) / levels_per_group;
        size_t our_group_degree = static_cast<size_t>(group_degree(group, eps));
        uintE upper_group = cur_level / levels_per_group;
        size_t upper_group_degree = static_cast<size_t>(group_degree(upper_group, eps));

        return num_up_star_neighbors >= our_group_degree && num_up_neighbors
            < upper_constant * upper_group_degree;
    }

    inline bool is_dirty(const size_t levels_per_group, double upper_constant, double eps) const {
      bool upper = upper_invariant(levels_per_group, upper_constant, eps);
      bool lower = lower_invariant(levels_per_group, eps);
      return !(upper && lower);
    }

  };

  size_t n;  // number of vertices
  size_t levels_per_group;  // number of inner-levels per group,  O(\log n) many.
  parlay::sequence<LDSVertex> L_seq;
  LDSVertex* L;

  LDS(size_t n) : n(n) {
    levels_per_group = ceil(log(n) / log(OnePlusEps));
    L_seq = parlay::sequence<LDSVertex>(n);
    L = L_seq.begin();
  }

  LDS(size_t _n, double _eps, double _delta) : n(_n), eps(_eps), delta(_delta){
    OnePlusEps = (1 + _eps);
    UpperConstant = 2 + ((double) 3 / _delta);
    levels_per_group = ceil(log(n) / log(OnePlusEps));
    L_seq = parlay::sequence<LDSVertex>(_n);
    L = L_seq.begin();
  }

  uintE get_level(uintE ngh) {
    return L[ngh].level;
  }

  bool edge_exists(edge_type e) {
    auto[u, v] = e;
    auto l_u = L[u].level;
    auto l_v = L[v].level;
    if (l_u < l_v) {  // look in up(u)
      if (L[u].up.contains(v)) {
        assert(L[v].down[l_u].contains(u));
      }
      return L[u].up.contains(v);
    } else {  // look in up(v)
      if (L[v].up.contains(u)) {
          assert(L[u].up.contains(v) || L[u].down[l_v].contains(v));
      }
      return L[v].up.contains(u);
    }
  }

  // Invariant checking for an edge e that we expect to exist
  bool check_both_directions(edge_type e) {
    auto [u, v] = e;
    auto l_u = L[u].level;
    auto l_v = L[v].level;
    bool ok = true;
    if (l_u < l_v) {  // look in up(u)
      ok &= L[u].up.contains(v); assert(ok);
      ok &= L[v].down[l_u].contains(u); assert(ok);
    } else if (l_v < l_u) {  // look in up(v)
      ok &= L[v].up.contains(u); assert(ok);
      ok &= L[u].down[l_v].contains(v); assert(ok);
    } else {  // (l_v == l_u)
      ok &= L[v].up.contains(u); assert(ok);
      ok &= L[u].up.contains(v); assert(ok);
    }
    return ok;
  }


  // Input: sequence of vertex_ids
  template <class Seq, class Levels>
  void update_levels(Seq&& possibly_dirty, Levels& levels) {
    using level_and_vtx = std::pair<uintE, uintE>;

    // Compute dirty vertices, which have either lower / upper threshold
    // violated. Output this vertex as dirty only if it is not already in the
    // bucketing structure.
    auto level_and_vtx_seq = parlay::delayed_seq<level_and_vtx>(possibly_dirty.size(), [&] (size_t i) {
      uintE v = possibly_dirty[i];
      uintE level = UINT_E_MAX;
      if (L[v].is_dirty(levels_per_group, UpperConstant, eps)) {
        auto our_level = L[v].level;
        // Return only if the vertex is not in the bucketing structure.
        if (our_level >= levels.size() || !levels[our_level].contains(v)) {
          level = our_level;
        }
      }
      return std::make_pair(level, v);
    });
    auto dirty = parlay::filter(level_and_vtx_seq, [&] (const level_and_vtx& lv) {
      return lv.first != UINT_E_MAX;
    });

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


    assert(level_starts.size() >= 2);
    uintE max_current_level = dirty[level_starts[level_starts.size() - 2]].first + 1;
    if (levels.size() < max_current_level) {
      levels.resize(max_current_level);
    }

    parallel_for(0, level_starts.size() - 1, [&] (size_t i) {
      uintE idx = level_starts[i];
      uintE level = dirty[idx].first;
      uintE num_in_level = level_starts[i+1] - idx;

      auto stuff_in_level = parlay::delayed_seq<uintE>(num_in_level, [&] (size_t j) {
        return dirty[idx + j].second;
      });

      levels[level].append(stuff_in_level);
    });
  }

  // todo: move into LDSVertex?
  // Neighbors is a slice of sorted (level, neighbor_id) pairs.
  template <class Neighbors>
  void insert_neighbors(uintE vtx, Neighbors neighbors) {
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
  }

  // Resize the tables
  // Can probably be optimized more
  template <class Neighbors>
  void delete_neighbors(uintE vtx, Neighbors neighbors) {
       // Delete neighbors from the adjacency list structures in parallel
       parallel_for(0, neighbors.size(), [&] (size_t i){
          auto [level_id, v] = neighbors[i];
          if (level_id != kUpLevel) {
              bool deleted = L[vtx].down[level_id].remove(v);
              assert(deleted);
          } else {
              bool deleted = L[vtx].up.remove(v);
              assert(deleted);
          }
      });

      parallel_for(0, neighbors.size(), [&] (size_t i){
        if ((i == 0) || neighbors[i].first != neighbors[i-1].first) { // start of a new level
            uintE level_id = neighbors[i].first;
            size_t j = i;
            for (; j < neighbors.size(); j++) {
                if (neighbors[j].first != level_id) break;
            }
            size_t num_deleted = j - i;
            if (level_id != kUpLevel) {
                L[vtx].down[level_id].resize_down(num_deleted);
            } else {
                L[vtx].up.resize_down(num_deleted);
            }
        }
      });
  }



  // returns the total number of moved vertices
  template <class Levels>
  size_t rebalance_insertions(Levels&& levels, size_t cur_level_id, size_t total_moved = 0) {
    if (cur_level_id >= levels.size())
      return total_moved;

    auto& cur_level = levels[cur_level_id];
    if (cur_level.num_elms() == 0)
      return rebalance_insertions(levels, cur_level_id+1, total_moved);

    // Figure out the desire_level for each vertex in cur_level.
    // Sets the desire level for each vertex with a valid desire_level in L.
    parallel_for(0, cur_level.size(), [&] (size_t i) {
      uintE v = cur_level.table[i];
      uintE desire_level = UINT_E_MAX;
      if (levelset::valid(v) && L[v].is_dirty(levels_per_group, UpperConstant, eps)) {
        assert(L[v].lower_invariant(levels_per_group, eps));
        assert(!L[v].upper_invariant(levels_per_group, UpperConstant, eps));

        desire_level = L[v].get_desire_level_upwards(v, L, levels_per_group, UpperConstant, eps);
        L[v].desire_level = desire_level;

        assert(L[v].level == cur_level_id);
        assert(desire_level > cur_level_id);
      }
    });

    auto vertex_seq = parlay::delayed_seq<uintE>(cur_level.size(), [&] (size_t i) {
      uintE u = cur_level.table[i];
      if (levelset::valid(u)) {
        if (L[u].desire_level == kNotMoving) {
          u = UINT_E_MAX;  // filter out, otherwise leave in
        }
      }
      return u;
    });
    // Dirty now contains (dl(v), v) pairs for all v \in the current level
    // moving to dl(v) (their desire level).
    auto dirty_seq = parlay::filter(vertex_seq, [&] (const uintE& u) {
      return u != UINT_E_MAX;
    });
    auto dirty = dirty_seq.begin();

    // todo: can get away with a uintE seq for dirty.
    parallel_for(0, dirty_seq.size(), [&] (size_t i) {
      uintE v = dirty[i];
      uintE desire_level = L[v].desire_level;
      L[v].down.resize(desire_level);
    });

    // Consider the neighbors N(u) of a vertex u moving from cur_level to dl(u) > cur_level.
    // - {v \in N(u) | dl(u) < l(v)}: move u from cur_level -> dl(u) in v's Down
    // - {v \in N(u) | cur_level < l(v) <= dl(u)}:
    //     remove u from v's Down and insert into v's Up
    // - {v \in N(u) | l(v) == cur_level}:
    //     if dl(v) <= dl(u): unaffected (u is already in v's Up)
    //     else (dl(u) < dl(v)): remove u from v's Up and insert into v's Down
    //
    // Let's emit a sequence of (u,v) edge tuples for all up_neighbors, and sort
    // by the recieving endpoint. The steps the recieving vertex does will
    // depend on the calculation above. Note that the following steps are very
    // similar to the batch_insertion code. Is there some clean way of merging?

    // Compute the number of neighbors we affect.
    auto degrees = parlay::map(parlay::make_slice(dirty_seq), [&] (auto v) {
      return L[v].num_up_neighbors();
    });
    size_t sum_degrees = parlay::scan_inplace(parlay::make_slice(degrees));

    // Write the affected (flipped) edges into an array of
    //   (affected_neighbor, moved_id)
    // edge pairs.
    auto flipped = sequence<edge_type>::uninitialized(sum_degrees);
    parallel_for(0, dirty_seq.size(), [&] (size_t i) {
      uintE v = dirty[i];
      size_t offset = degrees[i];
      size_t end_offset = (i == dirty_seq.size()-1) ? sum_degrees : degrees[i+1];

      auto output = flipped.cut(offset, end_offset);
      size_t written = L[v].emit_up_neighbors(output, v);
      assert(written == (end_offset - offset));

      // Remove neighbors in up with level < desire_level that are not moving.
      L[v].filter_up_neighbors(v, L);
    });

    // Sort based on the affected_neighbor. Note that there are no dup edges.
    auto compare_tup = [&] (const edge_type& l, const edge_type& r) { return l < r; };
    parlay::sort_inplace(parlay::make_slice(flipped), compare_tup);

    // Compute the starts of each (modified) vertex's new edges.
    auto bool_seq = parlay::delayed_seq<bool>(flipped.size() + 1, [&] (size_t i) {
      return (i == 0) || (i == flipped.size()) || (std::get<0>(flipped[i-1]) != std::get<0>(flipped[i]));
    });
    auto starts = parlay::pack_index(bool_seq);

    // Save the vertex ids (we will use this in update_levels).
    auto affected = sequence<uintE>::from_function(starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE vtx_id = flipped[idx].first;
      return vtx_id;
    });

    // Map over the vertices being modified. Overwrite their incident edges with
    // (level_id, neighbor_id) pairs, sort by level_id, resize each level to the
    // correct size, and then insert in parallel.
    parallel_for(0, starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE u = std::get<0>(flipped[idx]);
      uintE incoming_degree = starts[i+1] - starts[i];
      auto neighbors = parlay::make_slice(flipped.begin() + idx,
          flipped.begin() + idx + incoming_degree);

      uintE l_u = L[u].level;

      // (1) vtx (u) is a vertex moving from the current level.
      if (l_u == cur_level_id && L[u].desire_level != UINT_E_MAX) {
        uintE dl_u = L[u].desire_level;
        assert(dl_u != UINT_E_MAX);
        assert(dl_u > l_u);

        // Map the incident edges to (level, neighbor_id).
        parallel_for(0, incoming_degree, [&] (size_t off) {
          auto [u_dup, v] = neighbors[off];
          assert(u == u_dup);
          uintE dl_v = L[v].desire_level;  // using desire_level not level.
          if (dl_v >= dl_u) { dl_v = kUpLevel; }  // stay in up
          neighbors[off] = {dl_v, v};  // send to dl_v
        });

        // Sort neighbors by level.
        parlay::sort_inplace(neighbors);

        auto level_seq = parlay::delayed_seq<uintE>(neighbors.size(), [&] (size_t i) {
          return neighbors[i].first;
        });
        // first key greater than or equal to kUpLevel
        size_t upstart = parlay::internal::binary_search(level_seq, kUpLevel, std::less<uintE>());

        if (upstart > 0) {  // stuff to delete from L[u].up
          parallel_for(0, upstart, [&] (size_t j) {
            uintE v = neighbors[j].second;
            bool removed = L[u].up.remove(v);
            assert(removed);
          });
          L[u].up.resize_down(upstart);

          // No need to update stuff in [upstart, end), since they are already
          // in u's Up.
          auto lower_neighbors = neighbors.cut(0, upstart);

          // Insert the neighbors to their new locations.
          insert_neighbors(u, lower_neighbors);
        }
      } else if (l_u > cur_level_id) {

        // Map the incident edges to (level, neighbor_id).
        parallel_for(0, incoming_degree, [&] (size_t off) {
          auto [u_dup, v] = neighbors[off];
          assert(u == u_dup);
          uintE dl_v = L[v].desire_level;  // using desire_level not level.
          if (dl_v >= l_u) { dl_v = kUpLevel; }  // move to up
          neighbors[off] = {dl_v, v};  // send to dl_v
        });

        // Sort neighbors by level.
        parlay::sort_inplace(neighbors);

        parallel_for(0, neighbors.size(), [&] (size_t i) {
          uintE v = neighbors[i].second;
          bool removed = L[u].down[cur_level_id].remove(v);
          assert(removed);
        });
        // Every updated edge is removed from cur_level.
        L[u].down[cur_level_id].resize_down(incoming_degree);

        // Insert the neighbors to their new locations.
        insert_neighbors(u, neighbors);
      } else {
        assert(l_u == cur_level_id && L[u].desire_level == UINT_E_MAX);
        // Don't have to do anything for these vertices. They are staying put at
        // l_u, but their neighbors (already in u's Up) are moving to higher
        // levels (and thus staying in u's Up).
      }
    });

    // Update current level for the dirty vertices, and reset
    // the desire_level.
    parallel_for(0, dirty_seq.size(), [&] (size_t i) {
      uintE v = dirty[i];
      L[v].level = L[v].desire_level;
      L[v].desire_level = UINT_E_MAX;
    });

    update_levels(std::move(affected), levels);

    // TODO: update total_moved properly (not necessary for correctness, but
    // interesting for logging / experimental evaluation).
    return rebalance_insertions(levels, cur_level_id + 1, total_moved);
  }

  // Used to balance the level data structure for deletions
  // Returns the total number of moved vertices
  template <class Levels>
  size_t rebalance_deletions(Levels&& levels, size_t cur_level_id, size_t total_moved = 0) {
      if (cur_level_id >= levels.size()) {
          return total_moved;
      }

      // Get vertices that have desire level equal to the current desire level

      auto nodes_to_move = gbbs::sparse_set<uintE>();

      // For deletion, we need to figure out all vertices that want to move to
      // the current level.
      //
      // TODO: this is probably inefficient. We probably need to optimize this
      // further. Right now it just figures out the desire level brute-force for
      // every dirty vertex.

      auto level_resizes = parlay::sequence<uintE>(levels.size(), (uintE) 0);
      auto level_size = parlay::sequence<size_t>(levels.size());

      parallel_for(0, levels.size(), [&] (size_t i) {
        auto elements_this_level = parlay::sequence<size_t>(levels[i].size(), (size_t) 0);
        parallel_for(0, levels[i].size(), [&] (size_t j) {
            uintE v = levels[i].table[j];

            uintE desire_level = UINT_E_MAX;

            if (levelset::valid(v) && L[v].is_dirty(levels_per_group, UpperConstant, eps)) {
                assert(L[v].upper_invariant(levels_per_group, UpperConstant, eps));
                assert(!L[v].lower_invariant(levels_per_group, eps));

                desire_level = L[v].get_desire_level_downwards(v, L, levels_per_group, UpperConstant,
                        eps);
                L[v].desire_level = desire_level;

                assert(L[v].level = i);
                assert(desire_level < i);

                if (desire_level == cur_level_id) {
                    elements_this_level[j] = 1;
                }
            }
        });
        size_t num_this_level = parlay::scan_inplace(parlay::make_slice(elements_this_level));
        level_size[i] = num_this_level;
      });


      size_t num_to_move = parlay::scan_inplace(parlay::make_slice(level_size));

      nodes_to_move.resize(num_to_move);

      auto outer_level_sizes = parlay::sequence<size_t>(levels.size(), (size_t) 0);
      parallel_for(0, levels.size(), [&] (size_t i) {
        auto inner_level_sizes = parlay::sequence<size_t>(levels[i].size(), (size_t) 0);
        parallel_for(0, levels[i].size(), [&] (size_t j) {
            uintE v = levels[i].table[j];

            if (levelset::valid(v) && L[v].is_dirty(levels_per_group, UpperConstant, eps)) {
                if (L[v].desire_level == cur_level_id) {
                    nodes_to_move.insert(v);
                    levels[i].remove(v);
                    inner_level_sizes[j] = 1;
                }
            }
        });
        outer_level_sizes[i] = parlay::scan_inplace(parlay::make_slice(inner_level_sizes));
      });

      parallel_for (0, outer_level_sizes.size(), [&] (size_t i){
          levels[i].resize_down(outer_level_sizes[i]);
      });

      // Turn nodes_to_move into a sequence
      auto nodes_to_move_seq = nodes_to_move.entries();

      // Compute the number of neighbors we affect
      auto degrees = parlay::map(parlay::make_slice(nodes_to_move_seq), [&] (auto v) {
          auto desire_level = L[v].desire_level;
          assert(desire_level < L[v].level);
          return L[v].num_neighbors_higher_than_level(desire_level);
      });
      size_t sum_degrees = parlay::scan_inplace(parlay::make_slice(degrees));

      // Write the affected neighbors into an array of (affected_neighbor,
      // moved_id) edge pairs.
      auto flipped = sequence<edge_type>::uninitialized(sum_degrees);
      parallel_for(0, nodes_to_move_seq.size(), [&] (size_t i){
        uintE v = nodes_to_move_seq[i];
        size_t offset = degrees[i];
        size_t end_offset = ((i == nodes_to_move_seq.size() - 1) ? sum_degrees : degrees[i + 1]);

        auto output = flipped.cut(offset, end_offset);
        size_t written = L[v].emit_neighbors_higher_than_level(output, v, L[v].desire_level);
        assert(written == (end_offset - offset));
      });

      // Sort by neighbor vertex index
      auto compare_tup = [&] (const edge_type& l, const edge_type& r) {return l < r;};
      parlay::sort_inplace(parlay::make_slice(flipped), compare_tup);

      // Compute the starts of the neighbor vertex indices
      auto bool_seq = parlay::delayed_seq<bool>(flipped.size() + 1, [&] (size_t i) {
        return (i == 0) || (i == flipped.size()) ||
            (std::get<0>(flipped[i-1]) != std::get<0>(flipped[i]));
      });
      auto starts = parlay::pack_index(bool_seq);

      // Save the vertex ids of the vertices which did not move
      auto affected = sequence<uintE>::from_function(starts.size() - 1, [&] (size_t i) {
        size_t idx = starts[i];
        uintE vtx_id = flipped[idx].first;
        return vtx_id;
      });

      // First, update the down-levels of vertices that do not move (not in
      // nodes_to_move). Then, all vertices which moved to the same level should
      // be both in each other's up adjacency lists.
      parallel_for(0, starts.size() - 1, [&] (size_t i) {
        size_t idx = starts[i];
        uintE u = std::get<0>(flipped[idx]);
        if (!nodes_to_move.contains(u) && L[u].level > cur_level_id) {
            uintE incoming_degree = starts[i+1] - starts[i];
            auto neighbors = parlay::make_slice(flipped.begin() + idx,
                    flipped.begin() + idx + incoming_degree);

            // All vertices in neighbors are moving to level cur_level
            // Update down[cur_level] to contain these vertices
            assert(cur_level_id <= L[u].level);
            L[u].down[cur_level_id].resize(neighbors.size());
            parallel_for(0, neighbors.size(), [&] (size_t q) {
                L[u].down[cur_level_id].insert(neighbors[q].second);
            });

            // Remove the vertices that moved from their previous levels in
            // L[u].down.
            //
            auto my_level = L[u].level;
            auto neighbor_levels = sequence<size_t>::uninitialized(neighbors.size());
            parallel_for(0, neighbors.size(), [&] (size_t i) {
                auto neighbor_id = neighbors[i].second;
                auto neighbor = L[neighbor_id];
                auto neighbor_level = neighbor.level;
                assert(neighbor_level >= cur_level_id);

                if (neighbor_level < my_level) {
                    assert(L[u].down[neighbor_level].contains(neighbor_id));
                    L[u].down[neighbor_level].remove(neighbor_id);
                    neighbor_levels[i] = neighbor_level;
                } else {
                    assert(L[u].up.contains(neighbor_id));
                    L[u].up.remove(neighbor_id);
                    neighbor_levels[i] = my_level;
                }
            });

            // Get the num deleted by sorting the levels of all neighbors
            auto compare_tup = [&] (const size_t l, const size_t r) { return l < r; };
            parlay::sort_inplace(parlay::make_slice(neighbor_levels), compare_tup);
            auto new_bool_seq = parlay::delayed_seq<bool>(neighbor_levels.size() + 1, [&] (size_t i) {
                      return (i == 0) || (i == neighbor_levels.size()) ||
                             (neighbor_levels[i-1] != neighbor_levels[i]);
                 });
            auto new_starts = parlay::pack_index(new_bool_seq);

            parallel_for (0, new_starts.size() - 1, [&] (size_t i){
                size_t idx = new_starts[i];
                size_t n_level = neighbor_levels[idx];
                size_t num_deleted = new_starts[i+1] - idx;
                if (n_level == my_level)
                    L[u].up.resize_down(num_deleted);
                else
                    L[u].down[n_level].resize_down(num_deleted);
            });
        }
      });

      // Move vertices in nodes_to_move to cur_level. Update the data structures
      // of each moved vertex and neighbors in flipped.
      //
      // Re-sort the flipped edges by endpoint which moved
      auto flipped_reverse = sequence<edge_type>::uninitialized(sum_degrees);
      parallel_for(0, flipped_reverse.size(), [&] (size_t i){
        flipped_reverse[i] = std::make_pair(flipped[i].second, flipped[i].first);
      });

      // Compute the starts of the reverse flipped edges
      auto bool_seq_reverse = parlay::delayed_seq<bool>(flipped_reverse.size() + 1, [&] (size_t i){
        return (i == 0) || (i == flipped_reverse.size()) || (std::get<0>(flipped_reverse[i-1])
                != std::get<0>(flipped_reverse[i]));
      });
      auto reverse_starts = parlay::pack_index(bool_seq_reverse);

      // Update the data structures (vertices kept at each level) of each vertex
      // that moved.
      //
      //NOTE: this should be a parallel for loop (as above) but I was running
      //into concurrency issues.
      //for (size_t i = 0; i < reverse_starts.size() - 1; i++) {
      parallel_for(0, reverse_starts.size(), [&] (size_t i) {
        size_t idx = reverse_starts[i];
        uintE moved_vertex_v = std::get<0>(flipped_reverse[idx]);
        size_t idx_plus = reverse_starts[i+1];
        uintE next_v = std::get<0>(flipped_reverse[idx_plus]);
        assert(moved_vertex_v < next_v);
        uintE num_neighbors = reverse_starts[i+1] - reverse_starts[i];

        auto neighbors = parlay::make_slice(flipped_reverse.begin() + idx,
                flipped_reverse.begin() + idx + num_neighbors);

        auto my_level = L[moved_vertex_v].level;
        auto move_up_size = sequence<size_t>::uninitialized(neighbors.size());
        //parallel_for (0, neighbors.size(), [&] (size_t j) {
        for (size_t j = 0; j < neighbors.size(); j++) {
            auto neighbor_id = neighbors[j].second;
            auto neighbor = L[neighbor_id];
            auto neighbor_level = neighbor.level;
            assert(neighbor_level >= cur_level_id);
            assert(moved_vertex_v == neighbors[j].first);

            if (neighbor_level < my_level) {
                move_up_size[j] = 1;
                assert(L[moved_vertex_v].down[neighbor_level].contains(neighbor_id));
            } else {
                move_up_size[j] = 0;
            }
        }
        //});

        size_t indegree_sum = parlay::scan_inplace(parlay::make_slice(move_up_size));

        L[moved_vertex_v].up.resize(indegree_sum);
        //parallel_for(0, neighbors.size(), [&] (size_t k){
        for (size_t k = 0; k < neighbors.size(); k++) {
            if (L[neighbors[k].second].level < my_level) {
                L[moved_vertex_v].up.insert(neighbors[k].second);
            }
        }
        //});
      });

      // Update current level of the moved vertices and reset the desired level
      parallel_for(0, nodes_to_move_seq.size(), [&] (size_t i){
        uintE v = nodes_to_move_seq[i];
        L[v].level = cur_level_id;
        L[v].desire_level = UINT_E_MAX;
        L[v].down.resize(cur_level_id);
        assert(L[v].upper_invariant(levels_per_group, UpperConstant, eps));
        assert(L[v].lower_invariant(levels_per_group, eps));
      });

      // update the levels with neighbors
      update_levels(std::move(affected), levels);

      return rebalance_deletions(levels, cur_level_id + 1, total_moved + nodes_to_move_seq.size());
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

    // Save the vertex ids. The next step will overwrite the edge pairs to store
    // (neighbor, current_level).  (saving is not necessary if we modify + sort
    // in a single parallel loop)
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

      // Insert moving neighbors.
      insert_neighbors(vtx, neighbors);

    }, 1);

    // New edges are done being deleted. Update the level structure.
    // Interface: supply vertex seq -> process will settle everything.

    using dirty_elts = sparse_set<uintE>;
    sequence<dirty_elts> levels;

    // Place the affected vertices into levels based on their current level.
    update_levels(std::move(affected), levels);

    // Update the level structure (basically a sparse bucketing structure).
    size_t total_moved = rebalance_insertions(std::move(levels), 0);
  }

  template <class Seq>
  void batch_deletion(const Seq& deletions_unfiltered) {
    // Remove edges that do not exist in the graph.
    auto deletions_filtered = parlay::filter(parlay::make_slice(deletions_unfiltered),
            [&] (const edge_type& e) { return edge_exists(e); });


    // Duplicate the edges in both directions and sort.
    auto deletions_dup = sequence<edge_type>::uninitialized(2*deletions_filtered.size());
    parallel_for(0, deletions_filtered.size(), [&] (size_t i) {
        auto [u, v] = deletions_filtered[i];
        deletions_dup[2*i] = {u, v};
        deletions_dup[2*i + 1] = {v, u};
    });
    auto compare_tup = [&] (const edge_type& l, const edge_type& r) {return l < r;};
    parlay::sort_inplace(parlay::make_slice(deletions_dup), compare_tup);

    // Remove duplicate deletions.
    auto not_dup_seq = parlay::delayed_seq<bool>(deletions_dup.size(), [&](size_t i){
        auto [u, v] = deletions_dup[i];
        bool not_self_loop = u != v;
        return not_self_loop && ((i == 0) || (deletions_dup[i] != deletions_dup[i-1]));
    });
    auto deletions = parlay::pack(parlay::make_slice(deletions_dup), not_dup_seq);

    // Compute the starts of each (modified) vertex's deleted edges.
    auto bool_seq = parlay::delayed_seq<bool>(deletions.size() + 1, [&] (size_t i) {
      return (i == 0) || (i == deletions.size()) ||
        (std::get<0>(deletions[i-1]) != std::get<0>(deletions[i]));
    });
    auto starts = parlay::pack_index(bool_seq);

    // Save the vertex ids. The next step will overwrite the edge pairs to store
    // (neighbor, current_level).  (saving is not necessary if we modify + sort
    // in a single parallel loop)
    auto affected = sequence<uintE>::from_function(starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE vtx_id = std::get<0>(deletions[idx]);
      return vtx_id;
    });

    // Map over the vertices being modified. Overwrite their incident edges with
    // (level_id, neighbor_id) pairs, sort by level_id, resize each level to the
    // correct size, and then insert in parallel.
    parallel_for(0, starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE vtx = std::get<0>(deletions[idx]);
      uintE our_level = L[vtx].level;

      // Number of edges deleted that are adjacent to vtx
      uintE outgoing_degree = starts[i+1] - starts[i];

      // Neighbors incident to deleted edges
      auto neighbors = parlay::make_slice(deletions.begin() + idx,
          deletions.begin() + idx + outgoing_degree);

      // Map the incident edges to (level, neighbor_id).
      parallel_for(0, outgoing_degree, [&] (size_t off) {
        auto [u, v] = neighbors[off];
        assert(vtx == u);
        //assert(edge_exists({vtx, v}) || edge_exists({v, vtx}));
        uintE neighbor_level = L[v].level;
        if (neighbor_level >= our_level) { neighbor_level = kUpLevel; }
        neighbors[off] = {neighbor_level, v};
      });

      // Sort neighbors by level.
      parlay::sort_inplace(neighbors);

      // Delete moving neighbors.
      delete_neighbors(vtx, neighbors);

    }, 1);

    // New edges are done being inserted. Update the level structure.
    // Interface: supply vertex seq -> process will settle everything.

    using dirty_elts = sparse_set<uintE>;

    // Maintains the dirty vertices by their level
    sequence<dirty_elts> levels;

    // Place the affected vertices into levels based on their current level.
    // For deletion, need to move only vertices that are going to the lowest
    // level and then re-update the dirty vertices and repeat.
    // First, start by finding all the dirty vertices via the below method that
    // are adjacent to an edge deletion.
    update_levels(std::move(affected), levels);

    // Update the level structure (basically a sparse bucketing structure).
    size_t total_moved = rebalance_deletions(std::move(levels), 0);
  }

  void check_invariants() {
    bool invs_ok = true;
    for (size_t i=0; i<n; i++) {
      bool upper_ok = L[i].upper_invariant(levels_per_group, UpperConstant, eps);
      bool lower_ok = L[i].lower_invariant(levels_per_group, eps);
      assert(upper_ok);
      assert(lower_ok);
      invs_ok &= upper_ok;
      invs_ok &= lower_ok;
    }
    assert(invs_ok);
  }

  inline uintE group_for_level(uintE level) const {
    return level / levels_per_group;
  }

  uintE core(uintE v) const {
    auto l = L[v].level;
    uintE group = group_for_level(l);
    if (l % levels_per_group != levels_per_group - 1 && group != 0) group--;
    return ceil(L[v].group_degree(group, eps));
  }


  uintE max_coreness() const {
    auto levels = parlay::delayed_seq<uintE>(n, [&] (size_t i) {
        return L[i].level;
    });
    uintE max_level = pbbslib::reduce_max(levels);
    uintE max_group = group_for_level(max_level);
    return ceil(L[0].group_degree(max_group, eps));
  }
};

template <class Graph>
inline void RunLDS(Graph& G) {
  // using W = typename Graph::weight_type;
  size_t n = G.n;
  auto layers = LDS(n);

  auto edges = G.edges();

  size_t num_batches = 1000;
  size_t batch_size = edges.size() / num_batches;
  for (size_t i=0; i<num_batches; i++) {
    size_t start = batch_size*i;
    size_t end = std::min(start + batch_size, edges.size());

    auto batch = parlay::delayed_seq<LDS::edge_type>(end - start, [&] (size_t i) {
      uintE u = std::get<0>(edges[start + i]);
      uintE v = std::get<1>(edges[start + i]);
      return std::make_pair(u, v);
    });

    layers.batch_insertion(batch);
  }

//  for (size_t i = 0; i < n; i++) {
//    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
//      if (u < v) {
//        layers.insert_edge({u, v});
//      }
//    };
//    G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
//  }

  layers.check_invariants();

  for (size_t i=0; i<num_batches; i++) {
    size_t start = batch_size*i;
    size_t end = std::min(start + batch_size, edges.size());

    auto batch = parlay::delayed_seq<LDS::edge_type>(end - start, [&] (size_t i) {
      uintE u = std::get<0>(edges[start + i]);
      uintE v = std::get<1>(edges[start + i]);
      return std::make_pair(u, v);
    });

    layers.batch_deletion(batch);

    /*for (size_t i=0; i<batch.size(); i++) {
        bool exists = layers.edge_exists(batch[i]);
        //bool ok = layers.check_both_directions(batch[i]);
        assert(!exists);
        //assert(ok);
    }*/

//    for (size_t i=0; i<batch.size(); i++) {
//      bool exists = layers.edge_exists(batch[i]);
//      bool ok = layers.check_both_directions(batch[i]);
//      assert(exists);
//      assert(ok);
//    }
  }

  layers.check_invariants();
}

template <class W>
inline void RunLDS (BatchDynamicEdges<W>& batch_edge_list, long batch_size, bool compare_exact,
        LDS& layers) {
    auto batch = batch_edge_list.edges;
    for (size_t i = 0; i < batch.size(); i += batch_size) {
        timer t; t.start();
        auto end_size = std::min(i + batch_size, batch.size());
        auto insertions = parlay::filter(parlay::make_slice(batch.begin() + i,
                    batch.begin() + end_size), [&] (const DynamicEdge<W>& edge){
            return edge.insert;
        });

        auto deletions = parlay::filter(parlay::make_slice(batch.begin() + i,
                    batch.begin() + end_size), [&] (const DynamicEdge<W>& edge){
            return !edge.insert;
        });

        auto batch_insertions = parlay::delayed_seq<std::pair<uintE, uintE>>(insertions.size(),
                [&] (size_t i) {
            uintE vert1 = insertions[i].from;
            uintE vert2 = insertions[i].to;
            return std::make_pair(vert1, vert2);
        });


        auto batch_deletions = parlay::delayed_seq<std::pair<uintE, uintE>>(deletions.size(),
            [&] (size_t i) {
            uintE vert1 = deletions[i].from;
            uintE vert2 = deletions[i].to;
            return std::make_pair(vert1, vert2);
        });

        layers.batch_insertion(batch_insertions);

        layers.batch_deletion(batch_deletions);

        double tt = t.stop();
        std::cout << "### Batch Running Time: " << tt << std::endl;
        std::cout << "### Batch Num: " << std::min(batch.size(), i + batch_size) << std::endl;
        std::cout << "### Coreness Estimate: " << layers.max_coreness() << std::endl;
        if (compare_exact) {
            auto graph = dynamic_edge_list_to_symmetric_graph(batch_edge_list, std::min(batch.size(),
                        i + batch_size));

            // Run kcore on graph
            auto cores = KCore(graph, 16);

            auto max_core = parlay::reduce(cores, parlay::maxm<uintE>());
            std::cout << "### Coreness Exact: " << max_core << std::endl;

            // Compare cores[v] to layers.core(v)
            auto approximation_error = parlay::delayed_seq<float>(batch_edge_list.max_vertex,
                    [&] (size_t j) -> float {
                auto exact_core = j >= graph.n ? 0 : cores[j];
                auto approx_core = layers.core(j);
                if (exact_core == 0 || approx_core == 0) {
                    return 0;
                }
                return (exact_core > approx_core) ? (float) exact_core / (float) approx_core :
                       (float) approx_core / (float) exact_core;
            });

            double mult_appx = (2 + 2*layers.eps);
            float bad = parlay::reduce(parlay::delayed_seq<float>(batch_edge_list.max_vertex, [&](size_t j) -> float{
                auto true_core = j >= graph.n ? 0 : cores[j];
                auto appx_core = layers.core(j);
                return (appx_core > (mult_appx * true_core)) + (appx_core < (true_core/mult_appx));
            }), parlay::addm<float>());

            // Output min, max, and average error
            float sum_error = parlay::reduce(approximation_error, parlay::addm<float>());
            float max_error = parlay::reduce(approximation_error, parlay::maxm<float>());
            float min_error = parlay::reduce(approximation_error,
              parlay::make_monoid([](float l, float r){
                if (l == 0) return r;
                if (r == 0) return l;
                return std::min(r, l);
            }, (float) 0));
            float denominator = parlay::reduce(parlay::delayed_seq<float>(batch_edge_list.max_vertex,
                        [&] (size_t j) -> float{
                auto exact_core = j >= graph.n ? 0 : cores[j];
                auto approx_core = layers.core(j);
                return (exact_core != 0) && (approx_core != 0);
            }), parlay::addm<float>());
            auto avg_error = (denominator == 0) ? 0 : sum_error / denominator;
            std::cout << "### Num Bad: " << bad << std::endl;
            std::cout << "### Per Vertex Average Coreness Error: " << avg_error << std::endl; fflush(stdout);
            std::cout << "### Per Vertex Min Coreness Error: " << min_error << std::endl; fflush(stdout);
            std::cout << "### Per Vertex Max Coreness Error: " << max_error << std::endl; fflush(stdout);
        }
    }
}

template <class Graph, class W>
inline void RunLDS(Graph& G, BatchDynamicEdges<W> batch_edge_list, long batch_size,
        bool compare_exact, double eps, double delta) {
    uintE max_vertex = std::max(uintE{G.n}, batch_edge_list.max_vertex);
    auto layers = LDS(max_vertex, eps, delta);
    if (G.n > 0) RunLDS(G);
    if (batch_edge_list.max_vertex > 0) RunLDS(batch_edge_list, batch_size, compare_exact, layers);
}

}  // namespace gbbs
