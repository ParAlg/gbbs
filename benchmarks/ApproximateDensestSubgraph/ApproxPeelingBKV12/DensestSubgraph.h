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

#include "ligra/edge_map_reduce.h"
#include "ligra/ligra.h"

template <template <typename W> class vertex, class W>
void WorkEfficientDensestSubgraph(graph<vertex<W> >& GA, double epsilon = 0.001) {
  const size_t n = GA.n;
  auto em = EdgeMap<uintE, vertex, W>(GA, std::make_tuple(UINT_E_MAX, 0), (size_t)GA.m / 15);

  double density_multiplier = (1+epsilon); // note that this is not (2+eps), since the density we compute includes edges in both directions already.

  auto D = sequence<uintE>(n, [&](size_t i) { return GA.V[i].getOutDegree(); });
//  auto vertices_remaining = sequence<uintE>(n, [&] (size_t i) { return i; });
  auto vertices_remaining = pbbs::delayed_seq<uintE>(n, [&] (size_t i) { return i; });

  size_t round = 1;
  uintE* last_arr = nullptr;
  size_t remaining_offset = 0;
  size_t num_vertices_remaining = n;

  double max_density = 0.0;

  // First round
  {
    size_t edges_remaining = GA.m;
    // Update density
    double current_density = ((double)edges_remaining) / ((double)vertices_remaining.size());
    double target_density = (density_multiplier*((double)edges_remaining)) / ((double)vertices_remaining.size());
    debug(std::cout << "Target density on round " << round << " is " << target_density << " erm = " << edges_remaining << " vrm = " << vertices_remaining.size() << std::endl;
    std::cout << "Current density on round " << round << " is " << current_density << std::endl;);
    if (current_density > max_density) {
      max_density = current_density;
    }

    auto keep_seq = pbbs::delayed_seq<bool>(n, [&] (size_t i) {
      return !(D[i] <= target_density);
    });

    auto split_vtxs_m = pbbs::split_two(vertices_remaining, keep_seq);
    uintE* this_arr = split_vtxs_m.first.to_array();
    size_t num_removed = split_vtxs_m.second;
    auto vs = vertexSubset(n, num_removed, this_arr);
    debug(std::cout << "removing " << num_removed << " vertices" << std::endl;);

    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const Maybe<std::tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
      D[v] -= edgesRemoved;
      return Maybe<std::tuple<uintE,uintE>>();
    };

    auto moved = em.template edgeMapCount<uintE>(vs, apply_f);
    moved.del();

    round++;
    last_arr = this_arr;
    remaining_offset = num_removed;
    num_vertices_remaining -= num_removed;
    if (vs.dense()) {
      pbbs::free_array(vs.d);
    }
  }

  while (num_vertices_remaining > 0) {
    uintE* start = last_arr + remaining_offset;
    uintE* end = start + num_vertices_remaining;
    auto vtxs_remaining = pbbs::make_range(start, end);

    auto degree_f = [&] (size_t i) {
      uintE v = vtxs_remaining[i];
      return static_cast<size_t>(D[v]);
    };
    auto degree_seq = pbbslib::make_sequence<size_t>(vtxs_remaining.size(), degree_f);
    long edges_remaining = pbbslib::reduce_add(degree_seq);

    // Update density
    double current_density = ((double)edges_remaining) / ((double)vtxs_remaining.size());
    double target_density = (density_multiplier*((double)edges_remaining)) / ((double)vtxs_remaining.size());
    debug(std::cout << "Target density on round " << round << " is " << target_density << " erm = " << edges_remaining << " vrm = " << vtxs_remaining.size() << std::endl;
    std::cout << "Current density on round " << round << " is " << current_density << std::endl;);
    if (current_density > max_density) {
      max_density = current_density;
    }

    auto keep_seq = pbbs::delayed_seq<bool>(vtxs_remaining.size(), [&] (size_t i) {
      return !(D[vtxs_remaining[i]] <= target_density);
    });

    auto split_vtxs_m = pbbs::split_two(vtxs_remaining, keep_seq);
    uintE* this_arr = split_vtxs_m.first.to_array();
    size_t num_removed = split_vtxs_m.second;
    auto vs = vertexSubset(n, num_removed, this_arr);
    debug(std::cout << "removing " << num_removed << " vertices" << std::endl;);

    num_vertices_remaining -= num_removed;
    if (num_vertices_remaining > 0) {
      auto apply_f = [&](const std::tuple<uintE, uintE>& p)
          -> const Maybe<std::tuple<uintE, uintE> > {
        uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
        D[v] -= edgesRemoved;
        return Maybe<std::tuple<uintE,uintE>>();
      };

      auto moved = em.template edgeMapCount<uintE>(vs, apply_f, no_output);
      moved.del();
    }

    round++;
    pbbs::free_array(last_arr);
    last_arr = this_arr;
    remaining_offset = num_removed;
    if (vs.dense()) {
      pbbs::free_array(vs.d);
    }
  }

  if (last_arr) {
    pbbs::free_array(last_arr);
  }
  cout << "### Density of (2(1+\eps))-Densest Subgraph is: " << max_density << endl;
}
