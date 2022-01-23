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

#include "bench_utils.h"
#include "uf_utils.h"

#include "empty_starting_graph.h"
#include "subsample_starting_graph.h"

namespace gbbs {
namespace connectit {
template <class Graph, bool provides_initial_graph>
void jayanti_find_twotrysplit(Graph& G, size_t n,
                              sequence<incremental_update>& updates,
                              size_t batch_size, size_t insert_to_query,
                              size_t rounds, commandLine P) {
  run_multiple_jayanti_alg<Graph, find_twotrysplit, provides_initial_graph>(
      G, n, updates, batch_size, insert_to_query, rounds, P);
}

template <class Graph, bool provides_initial_graph>
void jayanti_find_simple(Graph& G, size_t n,
                         sequence<incremental_update>& updates,
                         size_t batch_size, size_t insert_to_query,
                         size_t rounds, commandLine P) {
  run_multiple_jayanti_alg<Graph, find_simple, provides_initial_graph>(
      G, n, updates, batch_size, insert_to_query, rounds, P);
}
}  // namespace connectit

template <class Graph, bool provides_initial_graph>
void run_all_tests(Graph& G, size_t n, sequence<incremental_update>& updates,
                   size_t batch_size, size_t insert_to_query, size_t rounds,
                   commandLine P) {
  run_tests<Graph, provides_initial_graph>(
      G, n, updates, batch_size, insert_to_query, rounds, P,
      connectit::jayanti_find_twotrysplit<Graph, provides_initial_graph>,
      {
          connectit::jayanti_find_twotrysplit<Graph, provides_initial_graph>,
          //        connectit::jayanti_find_simple<Graph,
          //        provides_initial_graph>
      });
}
}  // namespace gbbs

#include "empty_starting_graph.h"
#include "subsample_starting_graph.h"
