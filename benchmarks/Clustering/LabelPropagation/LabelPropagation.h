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
//
// This file provides a parallel implementation of (weighted) LabelPropagation.

#pragma once

#include <math.h>

#include "benchmarks/GraphColoring/Hasenplaugh14/GraphColoring.h"
#include "gbbs/gbbs.h"

namespace gbbs {

using label_type = int64_t;
constexpr label_type kInvalidLabel = std::numeric_limits<label_type>::max();

namespace internal {

template <class Graph>
label_type compute_new_color(Graph& G,
                             const parlay::sequence<label_type>& cur_labels,
                             gbbs::uintE node_id) {
  using Weight = typename Graph::weight_type;
  auto node = G.get_vertex(node_id);
  size_t degree = node.out_degree();

  // Collect the set of neighboring labels, and weights to them.
  using Elt = std::pair<label_type, double>;
  parlay::sequence<Elt> label_and_weight(degree);
  auto map_f = [&](const gbbs::uintE& u, const gbbs::uintE& v,
                   const Weight& weight, size_t index) {
    if constexpr (std::is_same_v<Weight, gbbs::empty>) {
      label_and_weight[index] = Elt{cur_labels[v], double{1}};
    } else {
      label_and_weight[index] = Elt{cur_labels[v], double{weight}};
    }
  };
  node.out_neighbors().map_with_index(map_f);

  // Sum up the weights.
  parlay::sort_inplace(label_and_weight);
  auto copy_f = [](Elt a, Elt b) {
    if (a.first == b.first) {
      return Elt{a.first, a.second + b.second};
    }
    return b;
  };
  Elt identity = {std::numeric_limits<label_type>::max(), 0};
  auto copy_monoid = parlay::make_monoid(copy_f, identity);
  parlay::scan_inclusive_inplace(label_and_weight, copy_monoid);

  // Collect the total weight per-label.
  auto label_ends =
      parlay::pack_index(parlay::delayed_seq<bool>(degree, [&](size_t i) {
        return (i == degree - 1) |
               (label_and_weight[i].first != label_and_weight[i + 1].first);
      }));
  auto weights_and_labels = parlay::map(label_ends, [&](size_t index) {
    return std::make_pair(label_and_weight[index].second,
                          label_and_weight[index].first);
  });
  parlay::sort_inplace(weights_and_labels, std::greater<>{});

  if (degree > 0) {
    return weights_and_labels[0].second;
  }
  return kInvalidLabel;
}

}  // namespace internal

// Expects an undirected (possibly weighted) graph.
template <class Graph>
parlay::sequence<label_type> LabelPropagation(
    Graph& G, const parlay::sequence<label_type>& initial_labels,
    size_t max_iters = 100,
    bool use_async = true,
    bool use_graph_coloring = false) {
  const uintE n = G.n;
  using Weight = typename Graph::weight_type;

  // Create two copies of the labels.
  parlay::sequence<label_type> cur_labels =
      parlay::sequence<label_type>(initial_labels);
  parlay::sequence<label_type> next_labels = cur_labels;

  parlay::sequence<bool> cur_active(n, true);
  parlay::sequence<bool> next_active(n, false);

  // Try the graph-coloring variant.
  using color = gbbs::uintE;
  parlay::sequence<color> coloring;
  parlay::sequence<std::pair<color, gbbs::uintE>> color_and_node;
  parlay::sequence<size_t> color_starts;
  if (use_graph_coloring) {
    // Note that if we are using graph coloring, we want to update the same
    // label set, so async should be set to true.
    use_async = true;
    coloring = Coloring(G);
    color_and_node = parlay::tabulate<std::pair<color, gbbs::uintE>>(
        coloring.size(), [&](gbbs::uintE i) { return std::make_pair(coloring[i], i); });
    parlay::sort_inplace(coloring);
    color_starts = parlay::pack_index(
        parlay::delayed_seq<bool>(color_and_node.size(), [&](size_t i) {
          return (i == 0) |
                 (color_and_node[i].first != color_and_node[i - 1].first);
        }));
  }

  std::cout << "Starting LabelPropagation. Parameters:" << std::endl;
  std::cout << "  - use_async = " << use_async << std::endl;
  std::cout << "  - use_graph_coloring = " << use_graph_coloring << std::endl;
  std::cout << "  - max_iters = " << max_iters << std::endl;

  size_t iter = 0;
  while (iter < max_iters) {
    std::cout << "Running iteration: " << iter << std::endl;
    std::atomic<bool> changed = false;

    auto process_node = [&](gbbs::uintE i) {
      // One of our neighbors changed in the last iteration (or this is
      // the first iteration). Recompute this node's label.
      if (cur_active[i]) {
        // Reset our flag for the next iteration.
        cur_active[i] = false;

        // Computes the next label based on cur_labels.
        label_type new_label = internal::compute_new_color(G, cur_labels, i);

        if (new_label != kInvalidLabel && cur_labels[i] != new_label) {
          // Set our label for the next iteration.
          if (use_async) {
            // If using async, just update the current label set.
            cur_labels[i] = new_label;
          } else {
            next_labels[i] = new_label;
          }
          // Mark our neighbors to be active in the next iteration.
          auto activate_f = [&](const gbbs::uintE& u, const gbbs::uintE& v,
                                const Weight& weight) {
            if (!next_active[v]) next_active[v] = true;
          };
          G.get_vertex(i).out_neighbors().map(activate_f);
          if (!changed) {
            changed = true;
          }
        }
      }
    };

    if (!use_graph_coloring) {
      // Map over all nodes.
      parlay::parallel_for(0, n, [&](size_t i) { process_node(i); });
    } else {
      // Map over each color set, one after the other.
      for (size_t c = 0; c < color_starts.size(); ++c) {
        size_t start_offset = color_starts[c];
        size_t end_offset = (c == color_starts.size() - 1)
                                ? color_and_node.size()
                                : color_starts[c + 1];
        // Map over all vertices of the same color in parallel.
        parlay::parallel_for(start_offset, end_offset, [&](size_t i) {
          gbbs::uintE node_id = color_and_node[i].second;
          process_node(node_id);
        });
      }
    }

    // Swap the current and next labels/active sets.
    if (!use_async) {
      std::swap(cur_labels, next_labels);
    }
    std::swap(cur_active, next_active);

    // Check convergence. If no labels changed in this iteration, quit.
    if (!changed) break;

    ++iter;
  }
  // for (size_t i=0; i < n; ++i) {
  //   std::cout << cur_labels[i] << std::endl;
  // }
  return cur_labels;
}

template <class Graph>
parlay::sequence<label_type> LabelPropagation(Graph& G) {
  auto initial_labels = parlay::tabulate(G.n, [&](label_type i) { return i; });
  return LabelPropagation(G, initial_labels);
}

}  // namespace gbbs
