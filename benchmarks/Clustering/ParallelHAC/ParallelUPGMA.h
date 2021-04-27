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

#include <queue>
#include <unordered_set>
#include <vector>

#include "gbbs/gbbs.h"
#include "ClusteredGraph.h"

namespace gbbs {
namespace clustering {

// Algorithm (high-level):
// 1. Bucket the edges into O(log_{1+epsilon}(W_max / W_min)) edge classes (buckets).
// 2. For the i-th bucket of edges:  // process edges from highest to smallest bucket
//    a. random-mate step, performing (1+epsilon)-close merges.
//    b. subset of edges are no longer in this bucket---move to a lower bucket
//    c. repeat until this bucket is fully contracted.

// Questions:
//   - How to represent dendrogram output? Array, as before?

  template <class Weights,
          // provides get_weight : () -> Weights::weight_type which is the
          // datatype that is stored for each edge incident to a _cluster_. This
          // could involve more than simply storing the underlying weight, or
          // could internally be a representation like gbbs::empty.
    template <class WW> class w_vertex, class IW>  // the weight type of the underlying graph
auto HAC(symmetric_graph<w_vertex, IW>& G, Weights& weights) {

  using clustered_graph =
      gbbs::clustering::clustered_graph<Weights, IW, w_vertex>;

}

}  // namespace clustering
}  // namespace gbbs
