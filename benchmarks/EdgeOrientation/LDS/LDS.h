#pragma once

#include <stack>
#include <unordered_set>
#include <vector>
#include <stdlib.h>
#include <math.h>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

namespace gbbs {

struct Level {
  std::unordered_set<uintE> set;
  std::vector<uintE> vector;
  static constexpr int high_boundary = 1000;
  static constexpr int low_boundary = 500;
  bool use_vector;

  Level() : use_vector(true) {}

  void set_to_vector() {
    for (const auto& v : set) {
      vector.push_back(v);
    }
    set.clear();
    use_vector = true;
  }

  void vector_to_set() {
    for (const auto& v : vector) {
      set.insert(v);
    }
    vector.clear();
    use_vector = false;
  }

  void erase(uintE v) {
    if (use_vector) vector.erase(std::find(vector.begin(), vector.end(), v));
    else {
      set.erase(set.find(v));
      if (set.size() < low_boundary) set_to_vector();
    }
  }

  template <class Iterator>
  void erase(Iterator it) {
    if (use_vector) vector.erase(it);
    else {
      set.erase(it);
      if (set.size() < low_boundary) set_to_vector();
    }
  }

  void insert(uintE v) {
    if (use_vector) {
      vector.push_back(v);
      if (vector.size() > high_boundary) vector_to_set();
    } else set.insert(v);
  }

  template <class F>
  void iterate(F f) {
    if (use_vector) {
      for (const auto& v : vector) {
        f(v);
      }
    } else {
      for (const auto& v : set) {
        f(v);
      }
    }
  }

  template <class F>
  void special_iterate(F f) {
    if (use_vector) {
      auto end = set.end();
      for (auto it = vector.begin(); it != vector.end();) {
        f(it, end);
      }
    } else {
      auto end = vector.end();
      for (auto it = set.begin(); it != set.end();) {
        f(end, it);
      }
    }
  }

  size_t size() const {
    if (use_vector) return vector.size();
    return set.size();
  }
};

inline double group_degree(size_t group, double epsilon) {
  return pow(1 + epsilon, group);
}

inline double upper_constant(double delta, bool optimized) {
  if (optimized)
    return 1.1;
  else
    return (2 + static_cast<double>(3) / delta);
}

struct LDS {
  size_t n;  // number of vertices
  double delta = 9.0;
  double epsilon = 3.0;
  bool optimized_insertion = false;

  size_t total_work;

  //using level = std::unordered_set<uintE>;
  using down_neighbors = std::vector<Level>;
  using up_neighbors = Level;
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
        down[l_v].erase(v);
      } else {
        up.erase(v);
      }
    }


    inline bool upper_invariant(const size_t levels_per_group, double epsilon, double delta,
            bool optimized) const {
      uintE group = level / levels_per_group;
      return up.size() <=
             static_cast<size_t>(upper_constant(delta, optimized) * group_degree(group, epsilon));
    }

    inline bool lower_invariant(const size_t levels_per_group, double epsilon) const {
      if (level == 0) return true;
      uintE lower_group = (level - 1) / levels_per_group;
      auto up_size = up.size();
      auto prev_level_size = down[level - 1].size();
      return (up_size + prev_level_size) >=
             static_cast<size_t>(
                 group_degree(lower_group, epsilon));
    }
  };

  // number of inner-levels per group,  O(\log n) many.
  size_t levels_per_group;
  parlay::sequence<LDSVertex> L;
  std::stack<uintE> Dirty;

  LDS(size_t _n, double _eps, double _delta, bool _optimized) : n(_n), epsilon(_eps),
    delta(_delta), optimized_insertion(_optimized) {
    levels_per_group = ceil(log(n) / log(1 + epsilon));
    // levels_per_group = parlay::log2_up(n);
    L = parlay::sequence<LDSVertex>(n);
  }

  uintE get_level(uintE ngh) { return L[ngh].level; }

  // Moving u from level to level - 1.
  template <class Levels>
  void level_decrease(uintE u, Levels& L) {
    total_work++;
    uintE level = L[u].level;
    auto& up = L[u].up;
    assert(level > 0);
    auto& prev_level = L[u].down[level - 1];

    prev_level.iterate([&](const uintE& ngh) {
      up.insert(ngh);
    });
    L[u].down.pop_back();  // delete the last level in u's structure.

    up.iterate([&](const uintE& ngh) {
    //for (const auto& ngh : up) {
      if (get_level(ngh) == level) {
        L[ngh].up.erase(u);
        L[ngh].down[level - 1].insert(u);

      } else if (get_level(ngh) > level) {
        L[ngh].down[level].erase(u);
        L[ngh].down[level - 1].insert(u);

        if (get_level(ngh) == level + 1) {
          Dirty.push(ngh);
        }
      } else {
        // u is still "up" for this ngh, no need to update.
        assert(get_level(ngh) == (level - 1));
      }
    });
    L[u].level--;  // decrease level
  }

  // Moving u from level to level + 1.
  template <class Levels>
  void level_increase(uintE u, Levels& L) {
    total_work++;
    uintE level = L[u].level;
    std::vector<uintE> same_level;
    auto& up = L[u].up;

    up.special_iterate([&](std::vector<uintE>::iterator& vec_it,
      std::unordered_set<uintE>::iterator& set_it) {
      bool use_vec = (vec_it != up.vector.end());
      uintE ngh = use_vec ? *vec_it : *set_it; //*it;
      if (L[ngh].level == level) {
        same_level.emplace_back(ngh);
        if (use_vec) vec_it = up.vector.erase(vec_it);
        else set_it = up.set.erase(set_it);
        // u is still "up" for this ngh, no need to update.
      } else {
        if (use_vec) vec_it++;
        else set_it++;
        // Must update ngh's accounting of u.
        if (L[ngh].level > level + 1) {
          L[ngh].down[level].erase(u);
          L[ngh].down[level + 1].insert(u);
        } else {
          assert(L[ngh].level == level + 1);
          L[ngh].down[level].erase(u);
          L[ngh].up.insert(u);

          Dirty.push(ngh);
        }
      }
    });
    // We've now split L[u].up into stuff in the same level (before the
    // update) and stuff in levels >= level + 1. Insert same_level elms
    // into down.
    auto& down = L[u].down;
    down.emplace_back(Level());
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
      if (!L[u].upper_invariant(levels_per_group, epsilon, delta, optimized_insertion)) {
        // Move u to level i+1.
        level_increase(u, L);
        Dirty.push(u);  // u might need to move up more levels.
      } else if (!L[u].lower_invariant(levels_per_group, epsilon)) {
        level_decrease(u, L);
        Dirty.push(u);  // u might need to move down more levels.
      }
    }
  }

  bool insert_edge(edge_type e) {
    auto[u, v] = e;
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
    auto[u, v] = e;
    auto l_u = L[u].level;
    auto l_v = L[v].level;
    L[u].remove_neighbor(v, l_v);
    L[v].remove_neighbor(u, l_u);

    Dirty.push(u);
    Dirty.push(v);
    fixup();
    return true;
  }

  void check_invariants() {
    bool invs_ok = true;
    for (size_t i = 0; i < n; i++) {
      bool upper_ok = L[i].upper_invariant(levels_per_group, epsilon, delta, optimized_insertion);
      bool lower_ok = L[i].lower_invariant(levels_per_group, epsilon);
      assert(upper_ok);
      assert(lower_ok);
      invs_ok &= upper_ok;
      invs_ok &= lower_ok;
    }
  }

  uintE max_coreness() {
    auto levels = parlay::delayed_seq<uintE>(n, [&] (size_t i) { return L[i].level; });
    uintE max_level = pbbslib::reduce_max(levels);
    uintE max_group = group_for_level(max_level);
    return group_degree(max_group, epsilon);
  }

  uintE core(uintE v) {
    auto l = L[v].level;
    uintE group = group_for_level(l);
    // If l is not the highest level in the group, we drop to the previous group
    if (l % levels_per_group != levels_per_group - 1 && group != 0) group--;
    return ceil(group_degree(group, epsilon));
  }

  inline uintE group_for_level(uintE level) const {
    return level / levels_per_group;
  }

  size_t get_size() {
    size_t size = 0;
    size += sizeof(delta) + sizeof(epsilon) + sizeof(optimized_insertion)
        + sizeof(levels_per_group) + sizeof(n);

    for (size_t i = 0; i < n; i++) {
        auto vertex = L[i];
        size += sizeof(vertex.level);
        for (size_t j = 0; j < vertex.down.size(); j++) {
            auto level = vertex.down[j];
            for (size_t k = 0; k < level.vector.size(); k++)
                size += sizeof(level.vector[k]);
            for (const auto& elem: level.set)
                size += sizeof(elem);
        }

        auto level = vertex.up;
        for (size_t k = 0; k < level.vector.size(); k++)
            size += sizeof(level.vector[k]);
        for (const auto& elem: level.set)
            size += sizeof(elem);
    }
    return size;
  }
};

template <class Graph>
inline void RunLDS(Graph& G, LDS& layers) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  for (size_t i = 0; i < n; i++) {
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (u < v) {
        layers.insert_edge({u, v});
      }
    };
    G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
  }

  std::cout << "Finished all insertions!" << std::endl;
  layers.check_invariants();
  std::cout << "Coreness estimate = " << layers.max_coreness() << std::endl;

  for (size_t i = 0; i < n; i++) {
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (u < v) {
        layers.delete_edge({u, v});
      }
    };
    G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
  }

  std::cout << "Finished all deletions!" << std::endl;
  layers.check_invariants();

  std::cout << "Coreness estimate = " << layers.max_coreness() << std::endl;

  size_t sum_lev = 0;
  for (size_t i = 0; i < G.n; i++) {
    sum_lev += layers.L[i].level;
  }
  std::cout << "sum_lev = " << sum_lev << std::endl;

  std::cout << "Total level increases and decreases: " << layers.total_work
            << std::endl;
}

template <class W>
inline void RunLDS(BatchDynamicEdges<W>& batch_edge_list, long batch_size,
  bool compare_exact, LDS& layers, size_t offset, bool get_size) {
  auto batch = batch_edge_list.edges;
  if (offset != 0) {
    for (size_t i = 0; i < offset; i++) {
      if (batch[i].insert) layers.insert_edge({batch[i].from, batch[i].to});
      else layers.delete_edge({batch[i].from, batch[i].to});
    }
  }

  for (size_t i = offset; i < batch.size(); i += batch_size) {
    timer t; t.start();
    for (size_t j = i; j < std::min(batch.size(), i + batch_size); j++) {
      if (batch[j].insert) layers.insert_edge({batch[j].from, batch[j].to});
      else layers.delete_edge({batch[j].from, batch[j].to});
    }

    double tt = t.stop();
    std::cout << "### Batch Running Time: " << tt << std::endl;
    std::cout << "### Batch Num: " <<
      std::min(batch.size(), i + batch_size) - offset << std::endl;
    std::cout << "### Coreness Estimate: " << layers.max_coreness() << std::endl;
    if (get_size) {
        auto size = layers.get_size();
        std::cout << "### Size: " << size << std::endl;
    }
    if (compare_exact) {
      auto graph = dynamic_edge_list_to_symmetric_graph(batch_edge_list,
        std::min(batch.size(), i + batch_size));

      // Run k-core on the graph
      auto cores = KCore(graph, 16);

      auto max_core = parlay::reduce(cores, parlay::maxm<uintE>());
      std::cout << "### Coreness Exact: " << max_core << std::endl;

      // Compare cores[v] to layers.core(v)
      auto approximation_error = parlay::delayed_seq<float>(
        batch_edge_list.max_vertex, [&](size_t j) -> float {
        auto exact_core = j >= graph.n ? 0 : cores[j];
        auto approx_core = layers.core(j);
        if (exact_core == 0 || approx_core == 0) {
          return 0;
        }
        return (exact_core > approx_core) ?
          (float) exact_core / (float) approx_core :
          (float) approx_core / (float) exact_core;
      });

      double mult_appx = (2 + 2*layers.epsilon);
      float bad = parlay::reduce(
        parlay::delayed_seq<float>(
          batch_edge_list.max_vertex, [&](size_t j) -> float{
            auto true_core = j >= graph.n ? 0 : cores[j];
            auto appx_core = layers.core(j);
            return (appx_core > (mult_appx*true_core)) +
              (appx_core < (true_core / mult_appx));
      }), parlay::addm<float>());


      // Output min, max, and average error
      float sum_error = parlay::reduce(approximation_error, parlay::addm<float>());
      float max_error = parlay::reduce(approximation_error, parlay::maxm<float>());
      float min_error = parlay::reduce(
        approximation_error, parlay::make_monoid([](float l, float r){
        if (l == 0) return r;
        if (r == 0) return l;
        return std::min(r, l);
      }, (float) 0));
      float denominator = parlay::reduce(parlay::delayed_seq<float>(batch_edge_list.max_vertex, [&](size_t j) -> float{
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
inline void RunLDS(Graph& G, BatchDynamicEdges<W> batch_edge_list,
  long batch_size, bool compare_exact, double eps, double delta,
  bool optimized_insertion, size_t offset, bool get_size) {
  uintE max_vertex = std::max(uintE{G.n}, batch_edge_list.max_vertex);
  auto layers = LDS(max_vertex, eps, delta, optimized_insertion);
  if (G.n > 0) RunLDS(G, layers);
  if (batch_edge_list.max_vertex > 0)
    RunLDS(batch_edge_list, batch_size, compare_exact, layers, offset, get_size);
}

}  // namespace gbbs
