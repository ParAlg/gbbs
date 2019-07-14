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

#include "edge_map_reduce.h"
#include "ligra.h"
#include "KCore.h"

// (2+2\epsilon)-appx DS for undirected graph (Bahmani et al.)
// DS = \nullset
// While some vertices remain:
//   1. peel vertices with degree < density
//   2. update density
//   3. If density(G_cur) > density(DS), update DS to G_cur.
// return DS

// Implements a work-inefficient version of the Bahmani et al. appx algorithm.
// The algorithm scans all vertices each round. The total work is O(m+n\log n).
template <template <typename W> class vertex, class W>
void WorkInefficientDensestSubgraph(graph<vertex<W> >& GA, double epsilon = 0.001) {
  const size_t n = GA.n;
  auto em = EdgeMap<uintE, vertex, W>(GA, std::make_tuple(UINT_E_MAX, 0), (size_t)GA.m / 15);

  double density_multiplier = (1+epsilon);

  auto bits = sequence<bool>(n, true);
  auto D = sequence<uintE>(n, [&](size_t i) { return GA.get_vertex(i).getOutDegree(); });

  long vertices_remaining = n;
  size_t round = 1;
  double max_density = 0.0;
  while (vertices_remaining > 0) {

    auto degree_f = [&] (size_t i) {
      return bits[i] ? static_cast<size_t>(D[i]) : static_cast<size_t>(0);
    };
    auto degree_seq = pbbslib::make_sequence<size_t>(n, degree_f);
    long edges_remaining = pbbslib::reduce_add(degree_seq);

    // update density
    double density = ((double)edges_remaining) / ((double)vertices_remaining);
    double target_density = (density_multiplier*((double)edges_remaining)) / ((double)vertices_remaining);
    debug(std::cout << "target density on round " << round << " is " << target_density << " erm = " << edges_remaining << " vrm = " << vertices_remaining << std::endl;
    std::cout << "density on round " << round << " is " << density << std::endl;);
    if (density > max_density) {
      max_density = density;
    }

    // filter out peeled vertices
    auto in_f = [&] (size_t i) {
      return i;
    };
    auto in_seq = pbbslib::make_sequence<uintE>(n, in_f);
    auto filter_low_deg = [&] (uintE i) {
      bool active = bits[i];
      if (active && D[i] <= target_density) {
        bits[i] = false;
        return true;
      };
      return false;
    };
    auto peeled = pbbslib::filter(in_seq, filter_low_deg);

    size_t vertices_removed = peeled.size();
    auto vs = vertexSubset(n, vertices_removed, peeled.to_array());
    debug(std::cout << "removing " << vertices_removed << " vertices" << std::endl;);
    vertices_remaining -= vertices_removed;

    if (vertices_remaining > 0) {
      auto apply_f = [&](const std::tuple<uintE, uintE>& p)
          -> const Maybe<std::tuple<uintE, uintE> > {
        uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
          D[v] -= edgesRemoved;
        return Maybe<std::tuple<uintE,uintE>>();
      };

      auto moved = em.template edgeMapCount<uintE>(vs, apply_f);
      moved.del();
    }

    vs.del();
    round++;
  }
  cout << "### Density of (2(1+\eps))-Densest Subgraph is: " << max_density << endl;
}

template <template <typename W> class vertex, class W>
void WorkEfficientDensestSubgraph(graph<vertex<W> >& GA, double epsilon = 0.001) {
  const size_t n = GA.n;
  auto em = EdgeMap<uintE, vertex, W>(GA, std::make_tuple(UINT_E_MAX, 0), (size_t)GA.m / 15);

  double density_multiplier = (1+epsilon); // note that this is not (2+eps), since the density we compute includes edges in both directions already.

  auto D = sequence<uintE>(n, [&](size_t i) { return GA.get_vertex(i).getOutDegree(); });
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

// Implements a parallel version of Charikar's 2-appx that runs in O(m+n)
// expected work and O(\rho\log n) depth w.h.p.
template <template <typename W> class vertex, class W>
void CharikarAppxDensestSubgraph(graph<vertex<W> >& GA) {
  // deg_ord = degeneracy_order(GA)
  // ## Now, density check for graph after removing each vertex, in the peeling-order.
  // Let S = stores 2*#edges to vertices > in degeneracy order. Note that 2* is
  //         needed since higher-ordered vertices don't have the edge to us.
  //
  // S = scan_add(S, fl_inplace | fl_reverse) ## reverse scan
  // density w/o vertex_i = S[i] / (n - i)
  // Compute the max over all v.

  size_t n = GA.n;
  auto degeneracy_order = DegeneracyOrder(GA);
  auto vtx_to_position = sequence<uintE>(n);

  parallel_for(0, n, [&] (size_t i) {
    uintE v = degeneracy_order.A[i];
    vtx_to_position[v] = i;
  });

  auto density_above = sequence<size_t>(n);

  par_for(0, n, 1, [&] (size_t i) {
    uintE pos_u = vtx_to_position[i];
    auto vtx_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      uintE pos_v = vtx_to_position[v];
      return pos_u < pos_v;
    };
#ifdef NVM
    density_above[pos_u] = 2*GA.get_vertex(i).countOutNgh(i, vtx_f, false);
#else
    density_above[pos_u] = 2*GA.get_vertex(i).countOutNgh(i, vtx_f);
#endif
  });

  size_t total_edges = pbbslib::scan_inplace(density_above.rslice(), pbbslib::addm<size_t>(),
      pbbslib::fl_inplace);
  if (total_edges != GA.m) {
    cout << "Assert failed: total_edges should be " << GA.m << " but is: " <<
      total_edges << endl;
    exit(0);
  }

  auto density_seq = pbbs::delayed_seq<double>(n, [&] (size_t i) {
    size_t dens = density_above[i];
    size_t rem = n - i;
    return static_cast<double>(dens) / static_cast<double>(rem);
  });
  double max_density = pbbslib::reduce_max(density_seq);
  cout << "### Density of 2-Densest Subgraph is: " << max_density << endl;
}


