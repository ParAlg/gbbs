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
#include "lt_utils.h"

#include "empty_starting_graph.h"
#include "subsample_starting_graph.h"

namespace connectit {
  template <class Graph, bool provides_initial_graph>
  void liutarjan_R(Graph& G, size_t n, pbbs::sequence<incremental_update>& updates, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P) {
    run_multiple_liu_tarjan_alg<Graph, parent_connect, root_update, shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
  }

  template <class Graph, bool provides_initial_graph>
  void liutarjan_RF(Graph& G, size_t n, pbbs::sequence<incremental_update>& updates, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P) {
    run_multiple_liu_tarjan_alg<Graph, parent_connect, root_update, full_shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
  }

  // (5)
  template <class Graph, bool provides_initial_graph>
  void liutarjan_RA(Graph& G, size_t n, pbbs::sequence<incremental_update>& updates, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P) {
    run_multiple_liu_tarjan_alg<Graph, simple_connect, root_update, shortcut, alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
  }

  // (6)
  template <class Graph, bool provides_initial_graph>
  void liutarjan_RFA(Graph& G, size_t n, pbbs::sequence<incremental_update>& updates, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P) {
    run_multiple_liu_tarjan_alg<Graph, simple_connect, root_update, full_shortcut, alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
  }

  // (3)
  template <class Graph, bool provides_initial_graph>
  void liutarjan_PRA(Graph& G, size_t n, pbbs::sequence<incremental_update>& updates, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P) {
    run_multiple_liu_tarjan_alg<Graph, parent_connect, root_update, shortcut, alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
  }

  // (4)
  template <class Graph, bool provides_initial_graph>
  void liutarjan_PRFA(Graph& G, size_t n, pbbs::sequence<incremental_update>& updates, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P) {
    run_multiple_liu_tarjan_alg<Graph, parent_connect, root_update, full_shortcut, alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
  }

} // namespace connectit

template <class Graph, bool provides_initial_graph, class F>
void run_tests(Graph& G, size_t n, pbbs::sequence<incremental_update>& updates, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P, F test, std::initializer_list<F> tests) {
  for (auto test : tests) {
#ifdef USE_PCM_LIB
  auto before_state = get_pcm_state();
  timer ot; ot.start();
#endif

  test(G, n, updates, batch_size, insert_to_query, rounds, P);

#ifdef USE_PCM_LIB
  double elapsed = ot.stop();
  auto after_state = get_pcm_state();
  cpu_stats stats = get_pcm_stats(before_state, after_state, elapsed, rounds);
  print_cpu_stats(stats, P);
#endif
  }
}


template <class Graph, bool provides_initial_graph>
void run_all_tests(Graph& G, size_t n, pbbs::sequence<incremental_update>& updates, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P) {
  run_tests<Graph, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P,
      connectit::liutarjan_R<Graph, provides_initial_graph>,
      {
        connectit::liutarjan_R<Graph, provides_initial_graph>,
        connectit::liutarjan_RF<Graph, provides_initial_graph>,
        connectit::liutarjan_RA<Graph, provides_initial_graph>,
        connectit::liutarjan_RFA<Graph, provides_initial_graph>,
        connectit::liutarjan_PRA<Graph, provides_initial_graph>,
        connectit::liutarjan_PRFA<Graph, provides_initial_graph>,
      });
}


#include "empty_starting_graph.h"
#include "subsample_starting_graph.h"
