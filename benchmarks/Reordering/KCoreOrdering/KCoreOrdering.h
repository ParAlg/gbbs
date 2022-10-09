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

#include "gbbs/gbbs.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/Reordering/reorder_cost.h"

namespace gbbs {

  // Given a graph G and an ordering pi, dump the graph in format X
  // using this ordering.
  // X = {adj, byte, bytepd}

  template <class Graph>
  void KCoreOrdering(const Graph& G, size_t num_buckets) {
    auto order = DegeneracyOrder(G, num_buckets);

    auto perm = parlay::sequence<uintE>(G.num_vertices());
    parlay::parallel_for(0, G.num_vertices(), [&] (size_t i) {
      uintE id = order[i];
      perm[id] = i;
    });

    PrintOrderingCosts(G, perm);
  }

}  // namespace gbbs
