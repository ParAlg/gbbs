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
    inline uintE get_desire_level_upwards(uintE vtx_id, Levels& L, const size_t levels_per_group) const {
      assert(!upper_invariant(levels_per_group));
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
            uintE up_degree = kUpperConstant*group_degree(group);

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

    inline uintE num_up_neighbors() const {
      return up.num_elms();
    }

    template <class OutputSeq>
    inline uintE emit_up_neighbors(OutputSeq output, uintE our_id) const {
      // Could filter using filter, but let's do this seq. for now...
      size_t off = 0;
      for (size_t i=0; i<up.table_seq.size(); i++) {
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

    inline double group_degree(size_t group) const {
      return pow(kOnePlusEpsilon, group);
    }

    inline bool upper_invariant(const size_t levels_per_group) const {
      uintE group = level / levels_per_group;
      uintE up_degree = kUpperConstant * group_degree(group);
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
  parlay::sequence<LDSVertex> L_seq;
  LDSVertex* L;

  LDS(size_t n) : n(n) {
    levels_per_group = ceil(log(n) / log(kOnePlusEpsilon));
    L_seq = parlay::sequence<LDSVertex>(n);
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
      return L[u].up.contains(v);
    } else {  // look in up(v)
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
      if (L[v].is_dirty(levels_per_group)) {
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
  }


  // returns the total number of moved vertices
  template <class Levels>
  size_t rebalance_insertions(Levels&& levels, size_t cur_level_id, size_t total_moved = 0) {
    // std::cout << "[RebalanceInsertions]: cur_level = " << cur_level_id << " levels_available = " << levels.size() << std::endl;
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
      if (levelset::valid(v) && L[v].is_dirty(levels_per_group)) {
        assert(L[v].lower_invariant(levels_per_group));
        assert(!L[v].upper_invariant(levels_per_group));

        desire_level = L[v].get_desire_level_upwards(v, L, levels_per_group);
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
      // std::cout << "Resized v = " << v << "'s down to " << desire_level << std::endl;
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

      // Insert moving neighbors.
      insert_neighbors(vtx, neighbors);

    }, 1);

    // New edges are done being inserted. Update the level structure.
    // Interface: supply vertex seq -> process will settle everything.

    using dirty_elts = sparse_set<uintE>;
    sequence<dirty_elts> levels;

    // Place the affected vertices into levels based on their current level.
    update_levels(std::move(affected), levels);

    // Update the level structure (basically a sparse bucketing structure).
    size_t total_moved = rebalance_insertions(std::move(levels), 0);

    // std::cout << "During insertions, " << total_moved << " vertices moved." << std::endl;
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

  auto edges = G.edges();

  size_t num_batches = 100;
  size_t batch_size = edges.size() / num_batches;
  for (size_t i=0; i<num_batches; i++) {
    std::cout << "===== Starting batch i = " << i << std::endl;
    size_t start = batch_size*i;
    size_t end = std::min(start + batch_size, edges.size());

    auto batch = parlay::delayed_seq<LDS::edge_type>(end - start, [&] (size_t i) {
      uintE u = std::get<0>(edges[start + i]);
      uintE v = std::get<1>(edges[start + i]);
      return std::make_pair(u, v);
    });

    layers.batch_insertion(batch);

//    for (size_t i=0; i<batch.size(); i++) {
//      bool exists = layers.edge_exists(batch[i]);
//      bool ok = layers.check_both_directions(batch[i]);
//      assert(exists);
//      assert(ok);
//    }
  }

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




}  // namespace gbbs
