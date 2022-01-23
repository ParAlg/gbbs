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
#include <stack>
#include <vector>
#include "ClusteredGraph.h"
#include "gbbs/gbbs.h"

namespace gbbs {
namespace nn_chain {

template <class Weights, class ClusterGraph>
void run_chain(Weights& weights, ClusterGraph& CG, std::stack<uintE>& chain,
               bool* on_stack) {
  assert(chain.size() > 0);

  debug(std::cout << std::endl;
        std::cout << "Starting chain from " << chain.top() << std::endl;);

  while (chain.size() > 0) {
    uintE top = chain.top();
    assert(on_stack[top]);

    auto edge_opt = CG.clusters[top].highest_priority_edge();
    assert(
        edge_opt.has_value());  // Must have at least one edge incident to it.

    auto nn = std::get<0>(*edge_opt);
    auto top_weight __attribute__((unused)) = std::get<1>(*edge_opt);

    debug(std::cout << "top is currently = " << top
                    << " nearest neighbor is nn = " << nn << std::endl;);
    if (on_stack[nn]) {
      // Found a reciprocal nearest neighbor.
      assert(chain.size() > 1);  // nn is definitely on the stack.
      chain.pop();               // Remove top.
      uintE r_nn = chain.top();  // Get our reciprocal nearest neighbor.

      // Multiple equal weight edges that go further back in the chain. It must
      // be the case that there is an equal weight edge to the previous
      // neighbor (r_nn), so assert this as a sanity check, and set nn to r_nn.
      if (r_nn != nn) {
        auto weight_opt __attribute__((unused)) =
            CG.clusters[top].neighbors.find(r_nn);
        nn = r_nn;
        // using ngh_map = typename ClusterGraph::neighbor_map;
        // auto top_entries = ngh_map::entries(CG.clusters[top].neighbors);
        assert(top_weight == *weight_opt);
      }
      chain.pop();  // remove nn.

      // Merge clusters top and nn.
      uintE merged_id __attribute__((unused)) = CG.unite(top, nn, top_weight);

      debug(std::cout << "Found reciprocal edge between top " << top << " ("
                      << CG.clusters[top].get_current_id() << ") and nn " << nn
                      << " (" << CG.clusters[nn].get_current_id()
                      << ") with weight: " << Weights::AsString(top_weight)
                      << std::endl;
            std::cout << "Done unite. Merged into " << merged_id << std::endl;);

      // Remove merged_id from on_stack.
      assert(merged_id == nn || merged_id == top);
      on_stack[top] = false;
      on_stack[nn] = false;
    } else {
      // nn not yet in chain. Push onto chain.
      chain.push(nn);
      on_stack[nn] = true;

      debug(std::cout << "Pushing onto stack from " << top << " ("
                      << CG.clusters[top].get_current_id() << ") and nn " << nn
                      << " (" << CG.clusters[nn].get_current_id()
                      << ") with weight: " << Weights::AsString(top_weight)
                      << std::endl;);
    }
  }
}

template <class Weights,
          // provides get_weight : () -> Weights::weight_type which is the
          // datatype that is stored for each edge incident to a _cluster_. This
          // could involve more than simply storing the underlying weight, or
          // could internally be a representation like gbbs::empty.
          template <class WW> class w_vertex,
          class IW>  // the weight type of the underlying graph
auto HAC(symmetric_graph<w_vertex, IW>& G, Weights& weights) {
  using W =
      typename Weights::weight_type;  // potentially a more complex type than IW

  using pq_elt = std::tuple<uintE, uintE, W>;
  using edge = std::tuple<uintE, W>;

  using clustered_graph =
      gbbs::clustering::clustered_graph<Weights, IW, w_vertex>;

  size_t n = G.n;

  // This object stores a representation of both the original clusters and the
  // clusters formed throughout the clustering process.
  auto CG = clustered_graph(G, weights);

  // PQ stores O(n) values --- one per cluster. The PQ values are (cluster,
  // cluster, weight) triples.
  auto pq_cmp = [](const pq_elt& l, const pq_elt& r) {
    return std::get<2>(l) < std::get<2>(r);
  };
  std::priority_queue<pq_elt, std::vector<pq_elt>, decltype(pq_cmp)> pq(pq_cmp);

  // Compute the min-weight edge incident to each original vertex
  parlay::sequence<edge> min_neighbors(
      n, std::make_tuple(UINT_E_MAX, std::numeric_limits<W>::max()));
  for (size_t i = 0; i < n; i++) {
    auto edge_option = CG.clusters[i].highest_priority_edge();
    if (edge_option.has_value()) {
      min_neighbors[i] = edge_option.value();
    }
  }

  sequence<bool> on_stack(n, false);

  debug(std::cout << "Starting clustering" << std::endl;);
  std::stack<uintE> chain;
  for (size_t v = 0; v < n; v++) {
    // Cluster with non-zero number of outgoing edges.
    if (CG.is_active(v) && CG.clusters[v].size() > 0) {
      debug(std::cout << "Starting new chain from v = " << v << std::endl;);
      chain.push(v);
      assert(!on_stack[v]);
      on_stack[v] = true;
      run_chain(weights, CG, chain, on_stack.begin());
    }
  }
  debug(std::cout << "Finished clustering" << std::endl;);

  return CG.get_dendrogram();
}

}  // namespace nn_chain
}  // namespace gbbs
