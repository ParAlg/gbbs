#pragma once

#include <unordered_set>
#include <stack>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "benchmarks/EdgeOrientation/ParallelLDS/LDS.h"

namespace gbbs {

struct Coloring {

  using color_palette = parlay::sequence<uintE>;
  using color = uintE;
  using lds = gbbs::LDS;
  using random = pbbs::random;

  static constexpr uintE kColorNotUsed = UINT_E_MAX;

  struct Palette {
    uintE v_; // Vertex for which this palette belongs
    uintE color_; // Vertex v's color
    uintE cur_color_index_; // Vertex v's current color index
    uintE size_; // Size of the palette
    uintE first_color_id_; // Color id of the first color in the palette
    color_palette colors_; // Palette containing occupied and free colors

    Palette(): v_(0), size_(0), color_(0) {}

    Palette(uintE vertex_id, uintE color_id,
            uintE palette_size, uintE first_color_id): v_(vertex_id), color_(color_id),
        size_(palette_size), first_color_id_(first_color_id) {
            colors_ = parlay::sequence<uintE>(size_);
        }

    inline uintE mark_color_occupied(uintE color_id, size_t num_occupied) {
        if (color_id >= first_color_id) {
            uintE palette_index = color_id - first_color_id_;
            if (palette_index <  size_) {
                colors_[palette_index] += num_occupied;
            }
        }
    }

    inline uintE unmark_color(uintE color_id, size_t num_unoccupied) {
            if (color_id >= first_color_id)  {
                    uintE palette_index = color_id - first_color_id_;
                    if (palette_index < size_) {
                        if (colors_[palette_index] > num_unoccupied)
                            colors_[palette_index] -= num_unoccupied;
                        else
                            colors[palette_index] = 0;
                    }
            }
    }

    inline uintE find_random_color(random r) {
        auto free_colors = parlay::filter(colors_, [&] (const uintE& num_occupied) {
            return num_occupied == 0;
        });

        auto num_free_colors = free_colors.size();

        auto r_v = r.fork(v_);
        auto random_index = r_v.rand();
        auto mask = (1 << pbbs::log2_up(num_free_colors)) - 1;

        random_index = random_index & mask;
        while (random_index >= num_free_colors) {
                r = r.next();
                r_v = r.fork(v_);
                random_index = r_v.rand();
                random_index = random_index & mask;
        }
        color_ = colors_[random_index];
        colors_[cur_color_index_] -= 1;
        colors_[random_index] += 1;
        return color_;
    }

    inline void shrink_palette(size_t shrinked_size) {
        colors_.resize(shrinked_size);
        size_ = shrinked_size;
    }

    inline void enlarge_palette(size_t enlarged_size) {
        colors_.resize(enlarged_size);
        parallel_for(size_, enlarged_size - 1, [&] (size_t i){
            colors_[i] = 0;
        });
        size_ = enlarged_size;
    }

    inline uintE get_color() {
        return color_;
    }

    inline uintE get_size() {
        return size_;
    }
  };

  size_t num_vertices;
  parlay::sequence<Palette> vertex_color_seq;
  Palette* VC;
  lds layers_of_vertices;

  Coloring(): num_vertices_(0) {}
  Coloring(size_t num_vertices, lds LDS): num_vertices_(num_vertices),
    layers_of_vertices_(LDS) {
    vertex_color_seq = parlay::sequence<LDSVertex>(num_vertices);
    VC = vertex_color_seq.begin();
  }

  template <class Seq>
  void batch_insertion(const Seq& insertions_unfiltered) {

  }

  template <class Seq>
  void batch_deletion(const Seq& deletions_unfiltered) {

  }

  void check_invariants() {
    bool invs_ok = true;
    for (size_t i=0; i<n; i++) {
      bool upper_ok = L[i].upper_invariant(levels_per_group, UpperConstant, eps, optimized_insertion);
      bool lower_ok = L[i].lower_invariant(levels_per_group, eps);
      assert(upper_ok);
      assert(lower_ok);
      invs_ok &= upper_ok;
      invs_ok &= lower_ok;
    }
    assert(invs_ok);
  }
};

template <class Graph>
inline void RunColoring(Graph& G, bool optimized_deletion) {
  // using W = typename Graph::weight_type;
  size_t n = G.n;
  auto layers = LDS(n, optimized_deletion);

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
        LDS& layers, bool optimized_insertion, size_t offset) {
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
        layers.batch_insertion(batch_insertions);
        layers.batch_deletion(batch_deletions);
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

        layers.batch_insertion(batch_insertions);
        double insertion_time = t.stop();

        t.start();
        layers.batch_deletion(batch_deletions);

        double deletion_time = t.stop();
        double tt = insertion_time + deletion_time;
        std::cout << "### Batch Running Time: " << tt << std::endl;
        std::cout << "### Insertion Running Time: " << insertion_time << std::endl;
        std::cout << "### Deletion Running Time: " << deletion_time << std::endl;
        std::cout << "### Batch Num: " << end_size - offset << std::endl;
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
        bool compare_exact, double eps, double delta, bool optimized_insertion, size_t offset) {
    uintE max_vertex = std::max(uintE{G.n}, batch_edge_list.max_vertex);
    auto layers = LDS(max_vertex, eps, delta, optimized_insertion);
    if (G.n > 0) RunLDS(G, optimized_insertion);
    if (batch_edge_list.max_vertex > 0) RunLDS(batch_edge_list, batch_size, compare_exact, layers, optimized_insertion, offset);
}

}  // namespace gbbs
