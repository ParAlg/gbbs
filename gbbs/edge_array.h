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

#include "bridge.h"
#include "macros.h"

namespace gbbs {

template <class W>
struct EdgeUtils {
  using weight_type = W;
  using edge = std::tuple<uintE, uintE, W>;

  // Given a set of edges (potentially asymmetric), ensure that each edge (u,v)
  // in the input sequence appears as both (u,v) and (v,u). This function also
  // sorts the output edges by id.
  static parlay::sequence<edge> undirect_and_sort(const parlay::sequence<edge>& edges) {
    // Duplicate the edges.
    sequence<edge> edge_sequence(edges.size() * 2);
    parlay::parallel_for(0, edges.size(), [&] (size_t i) {
      auto [u, v, wgh] = edges[i];
      edge_sequence[2 * i] = {u, v, wgh};
      edge_sequence[2 * i + 1] = {v, u, wgh};
    });
    // Sort the edges.
    parlay::sort_inplace(
      parlay::make_slice(edge_sequence), [](const edge& left, const edge& right) {
        return std::tie(std::get<0>(left), std::get<1>(left)) <
               std::tie(std::get<0>(right), std::get<1>(right));
      });
    // Filter duplicates.
    return filter_index(parlay::make_slice(edge_sequence), [&] (const edge& e, size_t i) {
      return (i == 0) ||
             (std::get<0>(edge_sequence[i-1]) != std::get<0>(e)) ||
             (std::get<1>(edge_sequence[i-1]) != std::get<1>(e));
    });
  }

  // Given a set of sorted edges as input, compute offsets to the first
  // occurence of each vertex
  static parlay::sequence<size_t> compute_offsets(size_t n, const parlay::sequence<edge>& edges) {
    auto offsets = parlay::sequence<size_t>(n);
    size_t m = edges.size();
    parlay::parallel_for(0, m, [&](size_t i) {
      if (i == 0 || (std::get<0>(edges[i]) != std::get<0>(edges[i - 1]))) {
        offsets[std::get<0>(edges[i])] = i;
      }
      if (i != (m - 1)) {
        size_t our_vtx = std::get<0>(edges[i]);
        size_t next_vtx = std::get<0>(edges[i + 1]);
        if (our_vtx != next_vtx && (our_vtx + 1 != next_vtx)) {
          parlay::parallel_for(our_vtx + 1, next_vtx,
             [&](size_t k) { offsets[k] = i + 1; });
        }
      }
      if (i == (m - 1)) { /* last edge */
        parlay::parallel_for(std::get<0>(edges[i]) + 1, offsets.size(),
           [&](size_t j) { offsets[j] = m; });
      }
    });
    return offsets;
  }

  template <class vertex>
  static typename vertex::neighbor_type* get_neighbors(const parlay::sequence<edge>& edges) {
    using neighbor_type = typename vertex::neighbor_type;
    size_t m = edges.size();
    auto neighbors = gbbs::new_array_no_init<neighbor_type>(m);
    parlay::parallel_for(0, m, [&] (size_t i) {
      neighbors[i] = {std::get<1>(edges[i]), std::get<2>(edges[i])};
    });
    return neighbors;
  }

};

// Edge Array Representation
template <class W>
struct edge_array {
  using weight_type = W;
  using edge = std::tuple<uintE, uintE, W>;

  // A sequence of edge tuples.
  sequence<edge> E;

  size_t n;  // num vertices.

  edge_array(sequence<edge>&& _E, size_t _n) : E(_E), n(_n) {}

  edge_array() {}

  size_t size() { return E.size(); }

  // Clears the edge array.
  sequence<edge>&& to_seq() {
    n = 0;
    return std::move(E);
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    size_t m = size();
    parallel_for(0, m, [&](size_t i) {
      uintE u, v;
      W w;
      std::tie(u, v, w) = E[i];
      f(u, v, w);
    });
  }
};

template <class W, class Graph>
inline edge_array<W> to_edge_array(Graph& G) {
  using edge = std::tuple<uintE, uintE, W>;

  size_t n = G.n;
  auto sizes = sequence<uintT>::uninitialized(n);
  parallel_for(0, n,
               [&](size_t i) { sizes[i] = G.get_vertex(i).out_degree(); });
  size_t m = parlay::scan_inplace(make_slice(sizes));
  assert(m == G.m);

  auto arr = sequence<edge>::uninitialized(m);
  parallel_for(0, n, [&](size_t i) {
    size_t idx = 0;
    uintT offset = sizes[i];
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      arr[offset + idx] = std::make_tuple(u, v, wgh);
      idx++;
    };
    G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
  });
  return edge_array<W>(std::move(arr), n);
}

}  // namespace gbbs
