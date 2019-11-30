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

#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"
#include "benchmarks/Connectivity/UnionFind/Connectivity.h"
#include "benchmarks/Connectivity/ShiloachVishkin/Connectivity.h"
#include "benchmarks/Connectivity/LabelPropagation/Connectivity.h"
#include "benchmarks/Connectivity/BFSCC/Connectivity.h"
#include "benchmarks/Connectivity/Framework/framework.h"
#include "bench_utils.h"

#include "ligra/graph_mutation.h"
#include "pbbslib/strings/string_basics.h"


// use by defining #EMPTY_STARTING_GRAPH
#include "empty_starting_graph.h"
// used otherwise
#include "subsample_starting_graph.h"








  template <class Graph, bool provides_initial_graph>
  double pick_test(Graph& G, size_t n, pbbs::sequence<std::tuple<uintE, uintE>>& updates, size_t id, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P) {
#ifdef USE_PCM_LIB
  auto before_state = get_pcm_state();
  timer ot; ot.start();
#endif

    switch (id) {
      /* {unite, unite_early, unite_nd} x {find_compress, find_naive, find_atomic_split, find_atomic_halve} */
      case 0:
        return run_multiple_uf_alg<Graph, unite, find_compress, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 1:
        return run_multiple_uf_alg<Graph, unite, find_naive, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 2:
        return run_multiple_uf_alg<Graph, unite, find_atomic_split, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 3:
        return run_multiple_uf_alg<Graph, unite, find_atomic_halve, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      case 4:
        return run_multiple_uf_alg<Graph, unite_early, find_compress, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 5:
        return run_multiple_uf_alg<Graph, unite_early, find_naive, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 6:
        return run_multiple_uf_alg<Graph, unite_early, find_atomic_split, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 7:
        return run_multiple_uf_alg<Graph, unite_early, find_atomic_halve, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      case 8:
        return run_multiple_uf_alg<Graph, unite_nd, find_compress, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 9:
        return run_multiple_uf_alg<Graph, unite_nd, find_naive, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 10:
        return run_multiple_uf_alg<Graph, unite_nd, find_atomic_split, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 11:
        return run_multiple_uf_alg<Graph, unite_nd, find_atomic_halve, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      /* unite_rem_lock x {find_naive, find_atomic_split, find_atomic_halve} */
      case 12:
        return run_multiple_uf_alg<Graph, unite_rem_lock, find_naive, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 13:
        return run_multiple_uf_alg<Graph, unite_rem_lock, find_atomic_split, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 14:
        return run_multiple_uf_alg<Graph, unite_rem_lock, find_atomic_halve, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      /* unite_rem_cas x {find_naive, find_atomic_split, find_atomic_halve} x {split_atomic_one, halve_atomic_one, splice_simple, splice_atomic} */
      case 15:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_split, split_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 16:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_split, halve_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 17:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_split, splice_simple, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 18:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_split, splice_atomic, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      case 19:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_halve, split_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 20:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_halve, halve_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 21:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_halve, splice_simple, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 22:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_halve, splice_atomic, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      case 23:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_naive, split_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 24:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_naive, halve_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 25:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_naive, splice_simple, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 26:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_naive, splice_atomic, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      /* Jayanti strategies */
      case 27:
        return run_multiple_jayanti_alg<Graph, find_twotrysplit, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds,  P);
      case 28:
        return run_multiple_jayanti_alg<Graph, find_simple, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds,  P);

      /* Liu-Tarjan algorithms */
      case 29:
        /* <parent_connect, update, shortcut> (Algorithm P) */
        return run_multiple_liu_tarjan_alg<Graph, parent_connect, simple_update, shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 30:
        /* <parent_connect, root_update, shortcut> (Algorithm R) */
        return run_multiple_liu_tarjan_alg<Graph, parent_connect, root_update, shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 31:
        /* <extended_connect, update, shortcut> (Algorithm E) */
        return run_multiple_liu_tarjan_alg<Graph, extended_connect, simple_update, shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 32:
        /* <parent_connect, update, full_shortcut> (PF) */
        return run_multiple_liu_tarjan_alg<Graph, parent_connect, simple_update, full_shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 33:
        /* <parent_connect, root_update, full_shortcut> (RF) */
        return run_multiple_liu_tarjan_alg<Graph, parent_connect, root_update, full_shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 34:
        /* <extended_connect, update, full_shortcut> (EF) */
        return run_multiple_liu_tarjan_alg<Graph, extended_connect, simple_update, full_shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      /* Shiloach-Vishkin */
      case 35:
        return run_multiple_shiloach_vishkin<Graph, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      default:
        std::cout << "Unknown test: " << id << std::endl;
        return 0.0 ;
    }
#ifdef USE_PCM_LIB
  double elapsed = ot.stop();
  auto after_state = get_pcm_state();
  cpu_stats stats = get_pcm_stats(before_state, after_state, elapsed, rounds);
  print_cpu_stats(stats, P);
#endif

  }
} // namespace connectit




/* Compatible algorithms:
 * - UF-variants
 * - Jayanti algorithms
 * - Shiloach-Vishkin
 * - Tarjan algorithms */

/* Modes:
 * starting from empty graph
 * starting with a full graph, sub-sample some % of it to use as updates */

constexpr int num_tests = 36; // update if new algorithm is added




// Two different modes: one starting without a base graph, and one starting with
// a base graph.
#ifdef EMPTY_STARTING_GRAPH


#else
generate_symmetric_once_main(Benchmark_nonempty_start, true);
#endif

