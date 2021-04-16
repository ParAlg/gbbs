#pragma once

#include <unordered_set>
#include <stack>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "benchmarks/EdgeOrientation/ParallelLDS/LDS.h"
#include "benchmarks/EdgeOrientation/ParallelLDS/sparse_set.h"

namespace gbbs {

struct DynamicColoring {

  using color_palette = parlay::sequence<uintE>;
  using color = uintE;
  using lds = gbbs::LDS;
  using random = parlay::random;
  using edge_type = std::pair<uintE, uintE>;
  using levelset = gbbs::sparse_set<uintE>;

  static constexpr uintE kColorNotUsed = UINT_E_MAX;

  struct Palette {
    uintE v_; // Vertex for which this palette belongs
    uintE color_; // Vertex v's color
    uintE size_; // Size of the palette
    uintE first_color_id_; // Color id of the first color in the palette
    color_palette colors_; // Palette containing occupied and free colors

    Palette(): v_(0), size_(0), color_(0), first_color_id_(0) {}

    Palette(uintE vertex_id, uintE color_id,
            uintE palette_size, uintE first_color_id): v_(vertex_id), color_(color_id),
        size_(palette_size), first_color_id_(first_color_id) {
        colors_ = parlay::sequence<uintE>(size_);
    }

    inline uintE mark_color(uintE color_id, size_t num_occupied) {
        if (color_id >= first_color_id_) {
            uintE palette_index = color_id - first_color_id_;
            if (palette_index < size_) {
                colors_[palette_index] += num_occupied;
            }
        }
    }

    inline uintE unmark_color(uintE color_id, size_t num_unoccupied) {
            if (color_id >= first_color_id_)  {
                    uintE palette_index = color_id - first_color_id_;
                    if (palette_index < size_) {
                        if (colors_[palette_index] > num_unoccupied)
                            colors_[palette_index] -= num_unoccupied;
                        else
                            colors_[palette_index] = 0;
                    }
            }
    }

    inline void unmark_all_colors() {
        parallel_for(0, colors_.size(), [&] (size_t i){
            colors_[i] = 0;
        });
    }

    inline uintE find_random_color(random r) {
        auto bool_colors = parlay::delayed_seq<bool>(colors_.size(), [&] (size_t i){
            return colors_[i] == 0;
        });

        auto free_colors = parlay::pack_index(bool_colors);

        auto num_free_colors = free_colors.size();

        auto r_v = r.fork(v_);
        auto random_index = r_v.rand();
        auto mask = (1 << parlay::log2_up(num_free_colors)) - 1;
        random_index = random_index & mask;
        while (random_index >= num_free_colors) {
                r = r.next();
                r_v = r.fork(v_);
                random_index = r_v.rand();
                random_index = random_index & mask;
        }
        /* Below could be useful later but not right now because we look at all
         * up neighbors again each time we move levels.
        if (color_ >= first_color_id) {
            colors_[color_ - first_color_id_] -= 1;
            colors_[free_colors[random_index]] += 1;
        }*/
        color_ = free_colors[random_index] + first_color_id_;
        return color_;
    }

    // TODO(qqliu): Maybe don't empty the palette for these two methods later.
    inline void shrink_palette(size_t shrinked_size) {
        colors_.resize(shrinked_size);
        // For now we empty the entire palette
        unmark_all_colors();
        size_ = shrinked_size;
    }

    inline void enlarge_palette(size_t enlarged_size) {
        colors_.resize(enlarged_size);
        // For now we empty the entire palette
        unmark_all_colors();
        size_ = enlarged_size;
    }

    inline uintE get_color() {
        return color_;
    }

    inline uintE get_size() {
        return size_;
    }
  };

  size_t num_vertices_;
  parlay::sequence<Palette> vertex_color_seq;
  Palette* VC;
  lds layers_of_vertices_;

  DynamicColoring(size_t num_vertices, lds layers): num_vertices_(num_vertices),
    layers_of_vertices_(layers) {
    vertex_color_seq = parlay::sequence<Palette>(num_vertices);
    VC = vertex_color_seq.begin();
    auto levels_per_group = layers_of_vertices_.levels_per_group;
    auto one_plus_epsilon = layers_of_vertices_.OnePlusEps;
    auto upper_constant = layers_of_vertices_.UpperConstant;

    parallel_for(0, num_vertices, [&] (size_t i){
        VC[i].v_ = i;
        VC[i].size_ = get_color_palette_size(0, levels_per_group, one_plus_epsilon, upper_constant);
        VC[i].colors_ = parlay::sequence<uintE>(VC[i].size_);
    });
  }

  size_t get_color_palette_size(size_t level, size_t levels_per_group,
          double one_plus_epsilon, double upper_constant) {
    auto group = level / levels_per_group;
    return 2 * upper_constant * pow(one_plus_epsilon, group);
  }

  size_t get_first_color(size_t level, size_t levels_per_group,
          double one_plus_epsilon, double upper_constant) {
    auto group = level / levels_per_group;
    auto additional_levels = level - group * levels_per_group;
    auto first_color = 0;
    if (group > 0) {
        first_color += 2 * levels_per_group * upper_constant * ((pow(one_plus_epsilon, group) - 1)/(one_plus_epsilon - 1));
    }
    first_color += additional_levels * get_color_palette_size(level, levels_per_group,
            one_plus_epsilon, upper_constant);
    return first_color;
  }

  // It is fine for layers higher than v to have the same color as v. Thus, we
  // can do 'lazy' coloring where if a vertex moves up to a higher layer, we
  // don't need to recolor it unless it conflicts with another color. Thus, the
  // invariant that no vertex v has neighbors below level(v) with colors from
  // v's palette is maintained.
  //
  // TODO(qqliu): This may need to be changed if we're doing a version with edge flips.
  template <class Seq>
  void batch_insertion(const Seq& insertions_unfiltered) {
    auto insertions = layers_of_vertices_.batch_insertion(insertions_unfiltered);

    // Finds the set of conflicts resulting from edge insertions.
    auto bool_conflict = parlay::delayed_seq<bool>(insertions.size(), [&] (size_t i) {
        auto v_1 = std::get<0>(insertions[i]);
        auto v_2 = std::get<1>(insertions[i]);
        if (v_1 < num_vertices_ && v_2 < num_vertices_) {
            return VC[v_1].color_ == VC[v_2].color_;
        } else {
            return false;
        }
    });

    // Edge insertions that result in conflicts.
    auto color_edges = parlay::pack(parlay::make_slice(insertions), bool_conflict);
    auto color_conflicts = parlay::sequence<uintE>(2*color_edges.size(), (uintE)0);
    // Gets set of vertices with conflicts.
    parallel_for(0, color_edges.size(), [&] (size_t i) {
        color_conflicts[i] = std::get<0>(color_edges[i]);
        color_conflicts[i + color_edges.size()] = std::get<1>(color_edges[i]);
    });

    // Remove duplicate vertices.
    parlay::integer_sort_inplace(parlay::make_slice(color_conflicts));
    auto dedup_conflicts = parlay::delayed_seq<bool>(color_conflicts.size(), [&] (size_t i){
        return (i == 0) || (color_conflicts[i] != color_conflicts[i-1]);
    });
    auto conflict_vertices = parlay::pack(parlay::make_slice(color_conflicts), dedup_conflicts);

    // Initialize random to randomly select vertices.
    auto r = random();
    while (conflict_vertices.size() > 0) {
        parallel_for(0, conflict_vertices.size(), [&] (size_t i) {
            auto v = conflict_vertices[i];
            auto level = layers_of_vertices_.get_level(v);
            auto levels_per_group = layers_of_vertices_.levels_per_group;
            auto one_plus_eps = layers_of_vertices_.OnePlusEps;
            auto upper_constant = layers_of_vertices_.UpperConstant;

            auto new_palette_size = get_color_palette_size(level, levels_per_group, one_plus_eps, upper_constant);
            auto first_color_palette = get_first_color(level, levels_per_group, one_plus_eps, upper_constant);

            VC[v].unmark_all_colors();
            if (first_color_palette > VC[v].first_color_id_) {
                VC[v].enlarge_palette(new_palette_size);
                VC[v].first_color_id_ = first_color_palette;
            }

            auto up_neighbors = layers_of_vertices_.L[v].up;
            assert(up_neighbors.num_elms() <= (VC[v].size_/2));
            auto possible_colors = parlay::sequence<uintE>(up_neighbors.size());
            parallel_for(0, up_neighbors.size(), [&] (size_t i){
                auto neighbor = up_neighbors.table[i];
                if (neighbor < num_vertices_) {
                    auto neighbor_color = VC[neighbor].color_;
                    possible_colors[i] = neighbor_color;
                } else {
                    possible_colors[i] = kColorNotUsed;
                }
            });

            // Filter out the real colors
            auto neighbor_colors = parlay::filter(possible_colors, [&] (const uintE& color) {
                return color != kColorNotUsed;
            });

            if (neighbor_colors.size() > 0) {
                // Sort based on the color
                parlay::sort_inplace(parlay::make_slice(neighbor_colors));
                auto bool_seq = parlay::delayed_seq<bool>(neighbor_colors.size() + 1, [&] (size_t i) {
                    return (i == 0) || (i == neighbor_colors.size()) || (neighbor_colors[i-1] != neighbor_colors[i]);
                });

                auto color_starts = parlay::pack_index(bool_seq);
                auto num_color = parlay::sequence<std::pair<uintE, uintE>>(color_starts.size() - 1);

                // Insert all neighbor colors and mark them in color array held by v
                parallel_for (0, color_starts.size() - 1, [&] (size_t q){
                    uintE color_index = color_starts[q];
                    uintE color = neighbor_colors[color_index];
                    num_color[q] = std::make_pair(color, color_starts[q+1] - color_index);
                });

                parallel_for(0, num_color.size(), [&] (size_t q){
                    uintE color = std::get<0>(num_color[q]);
                    uintE number = std::get<1>(num_color[q]);
                    VC[v].mark_color(color, number);
                });

                VC[v].find_random_color(r);
            }
        });

        auto next_conflict_vertices = parlay::filter(parlay::make_slice(conflict_vertices),
                [&] (const uintE& v){
            auto up_neighbors = layers_of_vertices_.L[v].up;
            auto bool_conflict = parlay::sequence<bool>(up_neighbors.size(), (bool) false);
            parallel_for(0, up_neighbors.size(), [&] (size_t i){
                auto neighbor = up_neighbors.table[i];
                if (neighbor < num_vertices_) {
                    auto neighbor_color = VC[neighbor].color_;
                    if (neighbor_color == VC[v].color_) {
                        bool_conflict[i] = true;
                    }
                }
            });
            auto any_conflict = parlay::filter(bool_conflict, [&] (const bool& conflict){
                return conflict;
            });
            if (any_conflict.size() > 0)
                return true;
            else
                return false;
        });
        r.next();
        conflict_vertices = next_conflict_vertices;
    }
  }

  template <class Seq>
  void batch_deletion(const Seq& deletions_unfiltered) {
    auto deletions = layers_of_vertices_.batch_deletion(deletions_unfiltered);

    // Compute all nodes that moved.
    auto moved_nodes = parlay::sequence<uintE>(2*deletions.size(), (uintE)0);
    parallel_for(0, deletions.size(), [&] (size_t i) {
        moved_nodes[i] = std::get<0>(deletions[i]);
        moved_nodes[i + deletions.size()] = std::get<1>(deletions[i]);
    });
    parlay::integer_sort_inplace(parlay::make_slice(moved_nodes));
    auto dedup_moved_nodes = parlay::delayed_seq<bool>(moved_nodes.size(), [&] (size_t i){
        return (i == 0) || (moved_nodes[i] != moved_nodes[i-1]);
    });
    auto reset_nodes = parlay::pack(parlay::make_slice(moved_nodes), dedup_moved_nodes);

    parallel_for(0, reset_nodes.size(), [&] (size_t i){
        uintE v = reset_nodes[i];
        auto old_first_color = VC[v].first_color_id_;
        auto new_level = layers_of_vertices_.get_level(v);
        auto one_plus_eps = layers_of_vertices_.OnePlusEps;
        auto upper_constant = layers_of_vertices_.UpperConstant;
        auto levels_per_group = layers_of_vertices_.levels_per_group;

        auto new_size = get_color_palette_size(new_level, levels_per_group, one_plus_eps, upper_constant);
        auto new_first_color = get_first_color(new_level, levels_per_group, one_plus_eps, upper_constant);
        // TODO(qqliu): We actually can't do lazy recoloring here, because the
        // invariant that no vertices lower than you can contain colors from
        // your palette is false. NEED TO FIX IMMEDIATELY.
        if (new_first_color < old_first_color){
            VC[v].size_ = new_size;
            VC[v].shrink_palette(new_size);
            VC[v].first_color_id_ = new_first_color;
        }
    });
    // TODO(qqliu): check that you don't need to recolor even after you
    // resize the palette.
  }

  void check_invariants() {
    bool invs_ok = true;
    for (size_t i = 0; i < num_vertices_; i++) {
      auto upper_neighbors = layers_of_vertices_.L[i].up;
      bool upper_ok = true;
      for (size_t j = 0; j < upper_neighbors.size(); j++) {
        auto neighbor = upper_neighbors.table[j];
        if (neighbor < num_vertices_) {
            if (VC[neighbor].color_ == VC[i].color_) {
                upper_ok = false;
            }
        }
      }

      bool lower_ok = true;
      auto lower_levels = layers_of_vertices_.L[i].down;
      for (size_t level = 0; level < lower_levels.size(); level++) {
        auto lower_neighbors = lower_levels[level];
        for (size_t j = 0; j < lower_neighbors.size(); j++) {
            auto neighbor = lower_neighbors.table[j];
            if (neighbor < num_vertices_) {
                if (VC[neighbor].color_ == VC[i].color_) {
                    lower_ok = false;
                }
            }
        }
      }
      assert(upper_ok);
      assert(lower_ok);
      invs_ok &= upper_ok;
      invs_ok &= lower_ok;
    }
    assert(invs_ok);
  }

  size_t num_unique_colors() {
    levelset colors = sparse_set<uintE>();
    for (size_t i = 0; i < num_vertices_; i++) {
        if (!colors.contains(VC[i].color_)) {
            colors.resize(1);
            colors.insert(VC[i].color_);
        }
    }
    return colors.num_elms();
  }
};

template <class Graph>
inline void RunDynamicColoring(Graph& G, bool optimized_deletion) {
  // using W = typename Graph::weight_type;
  size_t n = G.n;
  auto layers = LDS(n, optimized_deletion);
  auto coloring = DynamicColoring(n, layers);

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

    coloring.batch_insertion(batch);
    layers.check_invariants();
    coloring.check_invariants();
  }

  for (size_t i=0; i<num_batches; i++) {
    size_t start = batch_size*i;
    size_t end = std::min(start + batch_size, edges.size());

    auto batch = parlay::delayed_seq<LDS::edge_type>(end - start, [&] (size_t i) {
      uintE u = std::get<0>(edges[start + i]);
      uintE v = std::get<1>(edges[start + i]);
      return std::make_pair(u, v);
    });

    coloring.batch_deletion(batch);
    layers.check_invariants();
    coloring.check_invariants();
  }
}

template <class W>
inline void RunDynamicColoring (uintE n, BatchDynamicEdges<W>& batch_edge_list, long batch_size, bool compare_exact,
        LDS& layers, bool optimized_insertion, size_t offset) {
    auto coloring = DynamicColoring(n, layers);
    auto batch = batch_edge_list.edges;
    // First, insert / delete everything up to offset
    if (offset != 0) {
      for (size_t i = 0; i < offset; i += batch_size) {
        auto end_size = std::min(i + batch_size, offset);
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
        coloring.batch_insertion(batch_insertions);
        layers.check_invariants();
        coloring.check_invariants();

        coloring.batch_deletion(batch_deletions);
        layers.check_invariants();
        coloring.check_invariants();
      }
    }

    for (size_t i = offset; i < batch.size(); i += batch_size) {
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

        coloring.batch_insertion(batch_insertions);
        layers.check_invariants();
        coloring.check_invariants();

        double insertion_time = t.stop();

        t.start();
        coloring.batch_deletion(batch_deletions);
        layers.check_invariants();
        coloring.check_invariants();

        double deletion_time = t.stop();
        double tt = insertion_time + deletion_time;
        std::cout << "### Batch Running Time: " << tt << std::endl;
        std::cout << "### Insertion Running Time: " << insertion_time << std::endl;
        std::cout << "### Deletion Running Time: " << deletion_time << std::endl;
        std::cout << "### Batch Num: " << end_size - offset << std::endl;
        std::cout << "### Coreness Estimate: " << layers.max_coreness() << std::endl;
        std::cout << "### Number of Unique Colors: " << coloring.num_unique_colors() << std::endl;
    }
}

template <class Graph, class W>
inline void RunDynamicColoring(Graph& G, BatchDynamicEdges<W> batch_edge_list, long batch_size,
        bool compare_exact, double eps, double delta, bool optimized_insertion, size_t offset) {
    uintE max_vertex = std::max(uintE{G.n}, batch_edge_list.max_vertex);
    auto layers = LDS(max_vertex, eps, delta, optimized_insertion);
    if (G.n > 0) RunDynamicColoring(G, optimized_insertion);
    if (batch_edge_list.max_vertex > 0) RunDynamicColoring(max_vertex,
            batch_edge_list, batch_size, compare_exact, layers, optimized_insertion, offset);
}

}  // namespace gbbs
