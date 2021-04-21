#pragma once

#include <unordered_set>
#include <stack>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "benchmarks/EdgeOrientation/ParallelLDS/sparse_set.h"

namespace gbbs {

// Note: the following code is not always guaranteed to terminate.
// Perhaps we can find some bounds for random inputs but it is not guaranteed to
// terminate for worst-case inputs.
//
// palette_size_factor: determines whether we are choosing from a palette of
// size 2 * \delta + 1 or \delta + 1...etc.
struct RandomBlankColor{
    using random = parlay::random;
    using edge_type = std::pair<uintE, uintE>;
    using neighborset = gbbs::sparse_set<uintE>;
    using neighborsequence = gbbs::sequence<uintE>;

    static constexpr uintE kNotVertex = UINT_E_MAX;

    struct Vertex {
        uintE v_;
        uintE color_;
        neighborset neighbors_;

        Vertex(): color_(0) {
            neighbors_ = gbbs::sparse_set<uintE>();
        }

        Vertex(gbbs::sequence<uintE> new_neighbors, uintE color): color_(color) {
            neighbors_ = gbbs::sparse_set<uintE>();
            neighbors_.resize(new_neighbors.size());
            parallel_for(0, new_neighbors.size(), [&] (size_t i){
                neighbors_.insert(new_neighbors[i]);
            });
        }

        void remove_neighbors(neighborsequence neighbors) {
            parallel_for(0, neighbors.size(), [&] (size_t i){
                neighbors_.remove(neighbors[i]);
            });
            neighbors_.resize_down(neighbors.size());
        }

        void add_neighbors(neighborsequence neighbors) {
            neighbors_.resize(neighbors.size());
            parallel_for(0, neighbors.size(), [&] (size_t i){
                neighbors_.insert(neighbors[i]);
            });
        }

        bool receive_color(uintE color, random r, size_t palette_size_factor) {
            return color == color_;
        }

        bool inform_neighbors(Vertex neighbor, random r, size_t palette_size_factor) {
            return neighbor.receive_color(color_, r, palette_size_factor);
        }

        void find_new_color(random r, size_t palette_size_factor, neighborset unique_color) {
            parlay::sequence<bool> palette_colors = parlay::sequence<bool>(palette_size_factor * neighbors_.size() + 1);
            parallel_for(0, palette_size_factor * neighbors_.size() + 1, [&] (size_t i){
                if (!unique_color.contains(i))
                    palette_colors[i] = true;
                else
                    palette_colors[i] = false;
            });
            auto blank_colors = parlay::pack_index(palette_colors);
            auto r_v = r.fork(v_);
            auto random_index = r_v.rand();
            auto mask = (1 << parlay::log2_up(blank_colors.size())) - 1;
            random_index = random_index & mask;
            while (random_index >= blank_colors.size()) {
                r = r.next();
                r_v = r.fork(v_);
                random_index = r_v.rand();
                random_index = random_index & mask;
            }
            color_ = blank_colors[random_index];
        }
    };

    size_t palette_size_factor_;
    parlay::sequence<Vertex> vertex_colors;
    Vertex* VC;

    RandomBlankColor(size_t n, size_t palette_size_factor = 1): palette_size_factor_(palette_size_factor) {
        vertex_colors = parlay::sequence<Vertex>(n);
        VC = vertex_colors.begin();
    }

    bool edge_exists(edge_type e) {
        auto u = e.first;
        auto v = e.second;
        if (u < vertex_colors.size() && v < vertex_colors.size()) {
            if (VC[u].neighbors_.contains(v)) {
                assert(VC[v].neighbors_.contains(u));
                return true;
            } else {
                assert(!VC[v].neighbors_.contains(u));
                return false;
            }
        } else {
            return false;
        }
    }

    template <class Seq>
    void batch_insertion(const Seq& insertions_unfiltered) {
        auto insertions_filtered = parlay::filter(parlay::make_slice(insertions_unfiltered),
                        [&] (const edge_type& e) { return !edge_exists(e); });
        auto insertions_dup = sequence<edge_type>::uninitialized(2*insertions_filtered.size());
        parallel_for(0, insertions_filtered.size(), [&] (size_t i) {
                auto [u, v] = insertions_filtered[i];
                insertions_dup[2*i] = {u, v};
                insertions_dup[2*i + 1] = {v, u};
        });
        auto compare_tup = [&] (const edge_type& l, const edge_type& r) { return l < r; };
        parlay::sort_inplace(parlay::make_slice(insertions_dup), compare_tup);
        auto not_dup_seq = parlay::delayed_seq<bool>(insertions_dup.size(), [&] (size_t i) {
            auto [u, v] = insertions_dup[i];
            bool not_self_loop = u != v;
            return not_self_loop && ((i == 0) || (insertions_dup[i] != insertions_dup[i-1]));
        });
        auto insertions = parlay::pack(parlay::make_slice(insertions_dup), not_dup_seq);

        auto insertion_vertices = parlay::sequence<uintE>(2*insertions.size(), (uintE)0);
        parallel_for(0, insertions.size(), [&] (size_t i) {
            insertion_vertices[i] = std::get<0>(insertions[i]);
            insertion_vertices[i + insertions.size()] = std::get<1>(insertions[i]);
        });
        parlay::integer_sort_inplace(parlay::make_slice(insertion_vertices));
        auto bool_seq = parlay::delayed_seq<bool>(insertion_vertices.size() + 1, [&] (size_t i) {
            return (i == 0) || (i == insertion_vertices.size()) || (insertion_vertices[i-1] != insertion_vertices[i]);
        });
        auto vertex_starts = parlay::pack_index(bool_seq);
        // TODO(qqliu): Check to make sure this is once in either direction.
        // I.e. an edge insertion (u, v) should show up as (u, v) or (v, u) but
        // not both.
        auto num_insertions_per_vertex = parlay::sequence<std::pair<uintE, uintE>>(vertex_starts.size() - 1);
        parallel_for (0, vertex_starts.size() - 1, [&] (size_t q){
            uintE vertex_index = vertex_starts[q];
            uintE vertex = insertion_vertices[vertex_index];

            VC[vertex].neighbors_.resize(vertex_starts[q+1] - vertex_index);
        });

        parallel_for(0, insertions.size(), [&] (size_t i){
            auto edge = insertions[i];
            auto u = edge.first;
            auto v = edge.second;
            VC[u].neighbors_.insert(v);
            VC[v].neighbors_.insert(u);
        });

        auto color_edges = parlay::filter(parlay::make_slice(insertions),
            [&] (const edge_type& e) {
                auto u = e.first;
                auto v = e.second;
                return VC[u].color_ == VC[v].color_;
        });

        auto color_conflicts = parlay::sequence<uintE>(2*color_edges.size(), (uintE)0);
        parallel_for(0, color_edges.size(), [&] (size_t i) {
            color_conflicts[i] = std::get<0>(color_edges[i]);
            color_conflicts[i + color_edges.size()] = std::get<1>(color_edges[i]);
        });
        parlay::integer_sort_inplace(parlay::make_slice(color_conflicts));
        auto dedup_conflicts = parlay::delayed_seq<bool>(color_conflicts.size(), [&] (size_t i){
            return (i == 0) || (color_conflicts[i] != color_conflicts[i-1]);
        });
        auto conflict_vertices = parlay::pack(parlay::make_slice(color_conflicts), dedup_conflicts);

        auto r = random();
        while (conflict_vertices.size() > 0) {
            parallel_for(0, conflict_vertices.size(), [&] (size_t i){
                uintE v = conflict_vertices[i];
                neighborsequence occupied_colors = gbbs::sequence<uintE>(VC[v].neighbors_.size());
                parallel_for(0, VC[v].neighbors_.size(), [&] (size_t j){
                    auto neighbor = VC[v].neighbors_.table[j];
                    if (neighbor < vertex_colors.size()) {
                        occupied_colors[j] = VC[neighbor].color_;
                    } else {
                        occupied_colors[j] = kNotVertex;
                    }
                });
                parlay::integer_sort_inplace(parlay::make_slice(occupied_colors));
                auto bool_seq = parlay::delayed_seq<bool>(occupied_colors.size(),
                        [&] (size_t j) {
                    return ((j == 0) || (occupied_colors[j] != occupied_colors[j-1])) &&
                        (occupied_colors[j] != kNotVertex);
                });
                auto unique_color_starts = parlay::pack_index(bool_seq);
                neighborset unique_color = gbbs::sparse_set<uintE>();
                unique_color.resize(unique_color_starts.size());
                parallel_for(0, unique_color_starts.size(), [&] (size_t j){
                    auto color = occupied_colors[unique_color_starts[j]];
                    if (color < kNotVertex) {
                        unique_color.insert(color);
                    }
                });

                VC[v].find_new_color(r, palette_size_factor_, unique_color);
            });

            // TODO(qqliu): maybe this part can be more efficient.
            auto new_conflict_vertices = parlay::filter(parlay::make_slice(conflict_vertices),
                [&] (const uintE& u) {
                auto neighbors = VC[u].neighbors_;
                auto bool_conflict = parlay::sequence<bool>(neighbors.size(), (bool) false);
                parallel_for(0, neighbors.size(), [&] (size_t q){
                    auto neighbor = neighbors.table[q];
                    if (neighbor < vertex_colors.size()) {
                        auto neighbor_color = VC[neighbor].color_;
                        if (neighbor_color == VC[u].color_)
                            bool_conflict[q] = true;
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
            conflict_vertices = new_conflict_vertices;
            r = r.next();
        }
    }

    template <class Seq>
    void batch_deletion(const Seq& deletions_unfiltered) {
        auto deletions_filtered = parlay::filter(parlay::make_slice(deletions_unfiltered),
                        [&] (const edge_type& e) { return edge_exists(e); });
        auto deletions_dup = sequence<edge_type>::uninitialized(2*deletions_filtered.size());
        parallel_for(0, deletions_filtered.size(), [&] (size_t i) {
                auto [u, v] = deletions_filtered[i];
                deletions_dup[2*i] = {u, v};
                deletions_dup[2*i + 1] = {v, u};
        });
        auto compare_tup = [&] (const edge_type& l, const edge_type& r) { return l < r; };
        parlay::sort_inplace(parlay::make_slice(deletions_dup), compare_tup);
        auto not_dup_seq = parlay::delayed_seq<bool>(deletions_dup.size(), [&] (size_t i) {
            auto [u, v] = deletions_dup[i];
            bool not_self_loop = u != v;
            return not_self_loop && ((i == 0) || (deletions_dup[i] != deletions_dup[i-1]));
        });
        auto deletions = parlay::pack(parlay::make_slice(deletions_dup), not_dup_seq);
        parallel_for(0, deletions.size(), [&] (size_t i){
            edge_type edge = deletions[i];
        });

        parallel_for(0, deletions.size(), [&] (size_t i){
            auto edge = deletions[i];
            auto u = edge.first;
            auto v = edge.second;
            VC[u].neighbors_.remove(v);
            VC[v].neighbors_.remove(u);
        });

        auto deletion_vertices = parlay::sequence<uintE>(2*deletions.size(), (uintE)0);
        parallel_for(0, deletions.size(), [&] (size_t i) {
            deletion_vertices[i] = std::get<0>(deletions[i]);
            deletion_vertices[i + deletions.size()] = std::get<1>(deletions[i]);
        });
        parlay::integer_sort_inplace(parlay::make_slice(deletion_vertices));
        auto bool_seq = parlay::delayed_seq<bool>(deletion_vertices.size() + 1, [&] (size_t i) {
            return (i == 0) || (i == deletion_vertices.size()) || (deletion_vertices[i-1] != deletion_vertices[i]);
        });
        auto vertex_starts = parlay::pack_index(bool_seq);
        auto num_deletions_per_vertex = parlay::sequence<std::pair<uintE, uintE>>(vertex_starts.size() - 1);
        parallel_for(0, vertex_starts.size() - 1, [&] (size_t q){
            uintE vertex_index = vertex_starts[q];
            uintE vertex = deletion_vertices[vertex_index];

            VC[vertex].neighbors_.resize_down(vertex_starts[q+1] - vertex_index);
        });
    }

    void check_invariants() {
        bool invs_ok = true;
        for (size_t i = 0; i < vertex_colors.size(); i++) {
            auto neighbors = VC[i].neighbors_;
            parallel_for(0, neighbors.size(), [&] (size_t j){
                auto neighbor = neighbors.table[j];
                if (neighbor < vertex_colors.size()) {
                    if (VC[neighbor].color_ == VC[i].color_) {
                        invs_ok = false;
                    }
                }
            });
            assert(invs_ok);
        }
        assert(invs_ok);
    }

    // Probably should make the below parallel for faster runtime
    size_t num_unique_colors() {
        sparse_set<uintE> colors = sparse_set<uintE>();
        for (size_t i = 0; i < vertex_colors.size(); i++) {
            if (!colors.contains(VC[i].color_)) {
                colors.resize(1);
                colors.insert(VC[i].color_);
            }
        }
        return colors.num_elms();
    }
};

template <class Graph>
inline void RunDynamicColoring(Graph& G, size_t palette_size_factor) {
  // using W = typename Graph::weight_type;
  size_t n = G.n;
  auto coloring = RandomBlankColor(n, palette_size_factor);
  auto edges = G.edges();

  size_t num_batches = 1000;
  size_t batch_size = edges.size() / num_batches;
  for (size_t i=0; i<num_batches; i++) {
    size_t start = batch_size*i;
    size_t end = std::min(start + batch_size, edges.size());

    auto batch = parlay::delayed_seq<RandomBlankColor::edge_type>(end - start, [&] (size_t i) {
      uintE u = std::get<0>(edges[start + i]);
      uintE v = std::get<1>(edges[start + i]);
      return std::make_pair(u, v);
    });

    coloring.batch_insertion(batch);
    coloring.check_invariants();
  }

  for (size_t i=0; i<num_batches; i++) {
    size_t start = batch_size*i;
    size_t end = std::min(start + batch_size, edges.size());

    auto batch = parlay::delayed_seq<RandomBlankColor::edge_type>(end - start, [&] (size_t i) {
      uintE u = std::get<0>(edges[start + i]);
      uintE v = std::get<1>(edges[start + i]);
      return std::make_pair(u, v);
    });

    coloring.batch_deletion(batch);
    coloring.check_invariants();
  }
}

template <class W>
inline void RunDynamicColoring (uintE n, BatchDynamicEdges<W>& batch_edge_list, long batch_size,
        size_t offset, size_t palette_size_factor) {
    auto coloring = RandomBlankColor(n, palette_size_factor);
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
        coloring.check_invariants();

        coloring.batch_deletion(batch_deletions);
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
        coloring.check_invariants();

        double insertion_time = t.stop();

        t.start();
        coloring.batch_deletion(batch_deletions);
        coloring.check_invariants();

        double deletion_time = t.stop();
        double tt = insertion_time + deletion_time;
        std::cout << "### Batch Running Time: " << tt << std::endl;
        std::cout << "### Insertion Running Time: " << insertion_time << std::endl;
        std::cout << "### Deletion Running Time: " << deletion_time << std::endl;
        std::cout << "### Batch Num: " << end_size - offset << std::endl;
        std::cout << "### Number of Unique Colors: " << coloring.num_unique_colors() << std::endl;
    }
}

template <class Graph, class W>
inline void RunDynamicColoring(Graph& G, BatchDynamicEdges<W> batch_edge_list, long batch_size,
        size_t offset, size_t palette_size_factor) {
    uintE max_vertex = std::max(uintE{G.n}, batch_edge_list.max_vertex);
    if (G.n > 0) RunDynamicColoring(G, palette_size_factor);
    if (batch_edge_list.max_vertex > 0) RunDynamicColoring(max_vertex,
            batch_edge_list, batch_size, offset, palette_size_factor);
}
}  // namespace gbbs
