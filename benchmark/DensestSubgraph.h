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

#include "bucket.h"
#include "edge_map_reduce.h"
#include "ligra.h"

// Our goal: O(m + n) work and polylog depth densest subgraph.
// Q: is there a linear-work sequential algorithm?
//
// Charikar's sequential greedy 2-apx works as follows:
// DS = \nullset
// While some vertices remain:
//   1. peel the lowest degree vertex
//   2. update density
//   3. If density(G_cur) > density(DS), update DS to G_cur.
//
// Can implement this algorithm in linear time using Charikar's ideas:
// - store bucket data structure (n buckets)
// - min-bkt can drop by one in each iteration
// - cost of traversing empty buckets is charged to the edge that reduced the
//   min-bkt, or is at most O(n)
// therefore overall alg runs in O(m+n) time

// (2+2\epsilon)-appx DS for undirected graph (Bahmani et al.)
// DS = \nullset
// While some vertices remain:
//   1. peel vertices with degree < density
//   2. update density
//   3. If density(G_cur) > density(DS), update DS to G_cur.
// return DS
template <template <typename W> class vertex, class W>
void WorkInefficientDensestSubgraph(graph<vertex<W> >& GA, double epsilon = 0.001) {
  const size_t n = GA.n;
  auto em = EdgeMap<uintE, vertex, W>(GA, std::make_tuple(UINT_E_MAX, 0), (size_t)GA.m / 20);

  double density_multiplier = (2*(1+epsilon));

  auto bits = sequence<bool>(n, true);
  auto D = sequence<uintE>(n, [&](size_t i) { return GA.V[i].getOutDegree(); });

  long vertices_remaining = n;
  size_t round = 1;
  while (vertices_remaining > 0) {

    auto degree_f = [&] (size_t i) {
      return bits[i] ? static_cast<size_t>(D[i]) : static_cast<size_t>(0);
    };
    auto degree_seq = pbbslib::make_sequence<size_t>(n, degree_f);
    long edges_remaining = pbbslib::reduce_add(degree_seq);

    // update density
    double target_density = (density_multiplier*((double)edges_remaining)) / ((double)vertices_remaining);
    double density = ((double)edges_remaining) / ((double)vertices_remaining);
    std::cout << "target density on round " << round << " is " << target_density << " erm = " << edges_remaining << " vrm = " << vertices_remaining << std::endl;
    std::cout << "density on round " << round << " is " << density << std::endl;

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
    std::cout << "removing " << vertices_removed << " vertices" << std::endl;

    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const Maybe<std::tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
        D[v] -= edgesRemoved;
      return Maybe<std::tuple<uintE,uintE>>();
    };

    auto moved = em.template edgeMapCount<uintE>(vs, apply_f);
    moved.del();

    vertices_remaining -= vertices_removed;
    vs.del();
    round++;
  }
}

template <template <typename W> class vertex, class W>
void WorkEfficientDensestSubgraph(graph<vertex<W> >& GA, double epsilon = 0.001) {
  const size_t n = GA.n;
  auto em = EdgeMap<uintE, vertex, W>(GA, std::make_tuple(UINT_E_MAX, 0), (size_t)GA.m / 20);

  double density_multiplier = (2*(1+epsilon));

  auto D = sequence<uintE>(n, [&](size_t i) { return GA.V[i].getOutDegree(); });
  auto vertices_remaining = sequence<uintE>(n, [&] (size_t i) { return i; });

  size_t round = 1;
  while (vertices_remaining.size() > 0) {

    // Reduce over the remaining vertices. Note that we can skip this
    // computation on the first round.
    auto degree_f = [&] (size_t i) {
      uintE v = vertices_remaining[i];
      return static_cast<size_t>(D[v]);
    };
    auto degree_seq = pbbslib::make_sequence<size_t>(vertices_remaining.size(), degree_f);
    long edges_remaining = pbbslib::reduce_add(degree_seq);

    // Update density
    double target_density = (density_multiplier*((double)edges_remaining)) / ((double)vertices_remaining.size());
    double current_density = ((double)edges_remaining) / ((double)vertices_remaining.size());
    std::cout << "Target density on round " << round << " is " << target_density << " erm = " << edges_remaining << " vrm = " << vertices_remaining.size() << std::endl;
    std::cout << "Current density on round " << round << " is " << current_density << std::endl;

    // Filter out peeled vertices
    auto remove_seq = [&] (uintE v) {
      return D[v] <= target_density;
    };
    auto keep_seq = [&] (uintE v) {
      return D[v] > target_density;
    };
    // Can use split2; figure out later.
    auto peeled = pbbs::filter(vertices_remaining, remove_seq);
    auto vtx2 = pbbs::filter(vertices_remaining, keep_seq);
    vertices_remaining = std::move(vtx2);

    size_t peeled_size = peeled.size();
    auto vs = vertexSubset(n, peeled_size, peeled.to_array());
    std::cout << "removing " << peeled_size << " vertices" << std::endl;

    auto apply_f = [&](const std::tuple<uintE, uintE>& p)
        -> const Maybe<std::tuple<uintE, uintE> > {
      uintE v = std::get<0>(p), edgesRemoved = std::get<1>(p);
      D[v] -= edgesRemoved;
      return Maybe<std::tuple<uintE,uintE>>();
    };

    auto moved = em.template edgeMapCount<uintE>(vs, apply_f);
    moved.del();

    vs.del();
    round++;
  }
}
