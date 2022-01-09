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

#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"

namespace gbbs {
template <class Graph>
double WorkEfficientDensestSubgraph(Graph& G, double epsilon = 0.001) {
  const size_t n = G.n;
  auto em = hist_table<uintE, uintE>(std::make_tuple(UINT_E_MAX, 0),
                                     (size_t)G.m / 50);

  double density_multiplier =
      (1 + epsilon);  // note that this is not (2+eps), since the density we
                      // compute includes edges in both directions already.

  auto D = sequence<uintE>::from_function(
      n, [&](size_t i) { return G.get_vertex(i).out_degree(); });
  //  auto vertices_remaining = sequence<uintE>(n, [&] (size_t i) { return i;
  //  });
  auto vertices_remaining =
      parlay::delayed_seq<uintE>(n, [&](size_t i) { return i; });
  auto alive = sequence<bool>::from_function(n, [&](size_t i) { return true; });

  size_t round = 1;
  sequence<uintE> A;
  size_t remaining_offset = 0;
  size_t num_vertices_remaining = n;

  double max_density = 0.0;

  // First round
  {
    size_t edges_remaining = G.m;
    // Update density
    double current_density =
        ((double)edges_remaining) / ((double)vertices_remaining.size());
    double target_density = (density_multiplier * ((double)edges_remaining)) /
                            ((double)vertices_remaining.size());
    debug(std::cout << "Target density on round " << round << " is "
                    << target_density << " erm = " << edges_remaining
                    << " vrm = " << vertices_remaining.size() << std::endl;
          std::cout << "Current density on round " << round << " is "
                    << current_density << std::endl;);
    if (current_density > max_density) {
      max_density = current_density;
    }

    auto keep_seq = parlay::delayed_seq<bool>(
        n, [&](size_t i) { return !(D[i] <= target_density); });

    auto splits = parlay::split_two(vertices_remaining, keep_seq);
    A = std::move(splits.first);
    size_t num_removed = splits.second;
    debug(std::cout << "removing " << num_removed << " vertices" << std::endl;);

    auto removed = sequence<uintE>::uninitialized(num_removed);
    parallel_for(0, num_removed, [&](size_t i) {
      auto v = A[i];
      removed[i] = v;
      alive[v] = false;
    });
    auto vs = vertexSubset(n, std::move(removed));

    auto cond_f = [&](const uintE& u) { return alive[u]; };

    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const std::optional<std::tuple<uintE, uintE> > {
          uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
          D[v] -= edgesRemoved;
          return std::nullopt;
        };

    nghCount(G, vs, cond_f, apply_f, em, no_output);

    round++;
    remaining_offset = num_removed;
    num_vertices_remaining -= num_removed;
  }

  while (num_vertices_remaining > 0) {
    auto vtxs_remaining =
        A.cut(remaining_offset, remaining_offset + num_vertices_remaining);

    auto degree_f = [&](size_t i) {
      uintE v = vtxs_remaining[i];
      return static_cast<size_t>(D[v]);
    };
    auto degree_seq =
        parlay::delayed_seq<size_t>(vtxs_remaining.size(), degree_f);
    long edges_remaining = parlay::reduce(degree_seq);

    // Update density
    double current_density =
        ((double)edges_remaining) / ((double)vtxs_remaining.size());
    double target_density = (density_multiplier * ((double)edges_remaining)) /
                            ((double)vtxs_remaining.size());
    debug(std::cout << "Target density on round " << round << " is "
                    << target_density << " erm = " << edges_remaining
                    << " vrm = " << vtxs_remaining.size() << std::endl;
          std::cout << "Current density on round " << round << " is "
                    << current_density << std::endl;);
    if (current_density > max_density) {
      max_density = current_density;
    }

    auto keep_seq = parlay::delayed_seq<bool>(
        vtxs_remaining.size(),
        [&](size_t i) { return !(D[vtxs_remaining[i]] <= target_density); });

    auto split_vtxs_m = parlay::split_two(vtxs_remaining, keep_seq);
    A = std::move(split_vtxs_m.first);
    size_t num_removed = split_vtxs_m.second;
    debug(std::cout << "removing " << num_removed << " vertices" << std::endl;);

    auto removed = sequence<uintE>::uninitialized(num_removed);
    parallel_for(0, num_removed, [&](size_t i) {
      auto v = A[i];
      alive[v] = false;
      removed[i] = v;
    });
    auto vs = vertexSubset(n, std::move(removed));

    num_vertices_remaining -= num_removed;
    if (num_vertices_remaining > 0) {
      auto apply_f = [&](const std::tuple<uintE, uintE>& p)
          -> const std::optional<std::tuple<uintE, uintE> > {
            uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
            D[v] -= edgesRemoved;
            return std::nullopt;
          };

      auto cond_f = [&](const uintE& u) { return alive[u]; };
      nghCount(G, vs, cond_f, apply_f, em, no_output);
    }

    round++;
    remaining_offset = num_removed;
  }

  std::cout << "### Density of (2(1+\eps))-Densest Subgraph is: " << max_density
            << std::endl;
  return max_density;
}
}  // namespace gbbs
