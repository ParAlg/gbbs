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


static timer bt;
using uchar = unsigned char;

#define time(_var,_body)    \
  bt.start();               \
  _body;		    \
  double _var = bt.stop();

/* ************************* Connectivity wrappers *************************** */
template <class S1, class S2>
inline void cc_check(S1& correct, S2& check);

namespace connectit {

  template <class Graph, class Alg, bool provides_initial_graph>
  auto run_abstract_alg(
      Graph& G,
      size_t n,
      pbbs::sequence<std::tuple<uintE, uintE>>& updates,
      size_t batch_size,
      size_t insert_to_query,
      Alg& alg) {
    /* compute initial components */
    auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });

    alg.initialize(parents);
    if constexpr (provides_initial_graph) {
      alg.template compute_components<no_sampling>(parents);
    }

    /* process batches */
    timer tt; tt.start();
    size_t m = updates.size();
    size_t n_batches = (m + batch_size - 1) / batch_size;
    std::vector<double> batch_times;
    for (size_t i=0; i<n_batches; i++) {
      size_t start = i*batch_size;
      size_t end = std::min((i+1)*batch_size, m);
      auto update = updates.slice(start, end);

      timer tt; tt.start();
      alg.process_batch(parents, update, insert_to_query);
      double batch_time = tt.stop();
      batch_times.emplace_back(batch_time);
      // (optional correctness check after each batch)
    }
    double t = tt.stop();
    double med_batch_time = median(batch_times);

    double throughput = static_cast<double>(updates.size()) / t;
    return std::make_tuple(t, med_batch_time, throughput);
  }


  template<
    class Graph,
    UniteOption    unite_option,
    FindOption     find_option,
    bool provides_initial_graph>
  bool run_multiple_uf_alg(
      Graph& G,
      size_t n,
      pbbs::sequence<std::tuple<uintE, uintE>>& updates,
      size_t batch_size,
      size_t insert_to_query,
      size_t rounds,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine P) {
      /* Create initial parents array */

      auto find = get_find_function<find_option>();
      auto unite = get_unite_function<unite_option, decltype(find)>(n, find);
      using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
      auto alg = UF(G, unite, find);

      return run_abstract_alg<Graph, decltype(alg), provides_initial_graph>(G, n, updates, batch_size, insert_to_query, alg);
    };

    auto name = uf_options_to_string<no_sampling, find_option, unite_option>();
    return run_multiple(G, rounds, name, P, test);
  }

  template<
    class Graph,
    UniteOption    unite_option,
    FindOption     find_option,
    SpliceOption   splice_option,
    bool provides_initial_graph>
  bool run_multiple_uf_alg(
      Graph& G,
      size_t n,
      pbbs::sequence<std::tuple<uintE, uintE>>& updates,
      size_t batch_size,
      size_t insert_to_query,
      size_t rounds,
      commandLine& P) {
    static_assert(unite_option == unite_rem_cas);

    auto test = [&] (Graph& G, commandLine P) {
      auto find = get_find_function<find_option>();
      auto splice = get_splice_function<splice_option>();
      auto unite = unite_variants::UniteRemCAS<decltype(splice), decltype(find), find_option>(splice, find);
      using UF = union_find::UFAlgorithm<decltype(find), decltype(unite), Graph>;
      auto alg = UF(G, unite, find);

      return run_abstract_alg<Graph, decltype(alg), provides_initial_graph>(G, n, updates, batch_size, insert_to_query, alg);
    };

    auto name = uf_options_to_string<no_sampling, find_option, unite_option>();
    return run_multiple(G, rounds, name, P, test);
  }

  template <
    class Graph,
    JayantiFindOption find_option,
    bool provides_initial_graph>
  bool run_multiple_jayanti_alg(
      Graph& G,
      size_t n,
      pbbs::sequence<std::tuple<uintE, uintE>>& updates,
      size_t batch_size,
      size_t insert_to_query,
      size_t rounds,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine& P) {
      auto find = get_jayanti_find_function<find_option>();
      using UF = jayanti_rank::JayantiTBUnite<Graph, decltype(find)>;
      auto alg = UF(G, n, find);
      return run_abstract_alg<Graph, decltype(alg), provides_initial_graph>(G, n, updates, batch_size, insert_to_query, alg);
    };

    auto name = jayanti_options_to_string<no_sampling, find_option>();
    return run_multiple(G, rounds, name, P, test);
  }

  template <
    class Graph,
    LiuTarjanConnectOption  connect_option,
    LiuTarjanUpdateOption   update_option,
    LiuTarjanShortcutOption shortcut_option,
    LiuTarjanAlterOption    alter_option,
    bool provides_initial_graph>
  bool run_multiple_liu_tarjan_alg(
      Graph& G,
      size_t n,
      pbbs::sequence<std::tuple<uintE, uintE>>& updates,
      size_t batch_size,
      size_t insert_to_query,
      size_t rounds,
      commandLine& P) {

    auto test = [&] (Graph& G, commandLine& P) {
      auto connect = lt::get_connect_function<connect_option>();
      auto update = lt::get_update_function<update_option>();
      auto shortcut = lt::get_shortcut_function<shortcut_option>();

      static_assert(alter_option == no_alter);

      using LT = lt::LiuTarjanAlgorithm<
        decltype(connect),
        connect_option,
        decltype(update),
        update_option,
        decltype(shortcut),
        shortcut_option,
        Graph>;
      auto alg = LT(G, n, connect, update, shortcut);
      return run_abstract_alg<Graph, decltype(alg), provides_initial_graph>(G, n, updates, batch_size, insert_to_query, alg);
    };

    auto name = liu_tarjan_options_to_string<no_sampling,connect_option,update_option,shortcut_option,alter_option>();
    return run_multiple(G, rounds, name, P, test);
  }


  template <class Graph, bool provides_initial_graph>
  bool run_multiple_shiloach_vishkin(
      Graph& G,
      size_t n,
      pbbs::sequence<std::tuple<uintE, uintE>>& updates,
      size_t batch_size,
      size_t insert_to_query,
      size_t rounds,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine& P) {
      auto alg = shiloachvishkin_cc::SVAlgorithm<Graph>(G);
      return run_abstract_alg<Graph, decltype(alg), provides_initial_graph>(G, n, updates, batch_size, insert_to_query, alg);
    };

    auto name = "shiloach_vishkin";
    return run_multiple(G, rounds, name, P, test);
  }



  template <class Graph, bool provides_initial_graph>
  double pick_test(Graph& G, size_t n, pbbs::sequence<std::tuple<uintE, uintE>>& updates, size_t id, size_t batch_size, size_t insert_to_query, size_t rounds, commandLine P) {
#ifdef USE_PCM_LIB
  auto before_state = get_pcm_state();
  timer ot; ot.start();
#endif

    switch (id) {
      case 1:
        return run_multiple_uf_alg<Graph, unite, find_compress, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 2:
        return run_multiple_uf_alg<Graph, unite, find_naive, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 3:
        return run_multiple_uf_alg<Graph, unite, find_atomic_split, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 4:
        return run_multiple_uf_alg<Graph, unite, find_atomic_halve, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      case 5:
        return run_multiple_uf_alg<Graph, unite_early, find_compress, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 6:
        return run_multiple_uf_alg<Graph, unite_early, find_naive, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 7:
        return run_multiple_uf_alg<Graph, unite_early, find_atomic_split, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 8:
        return run_multiple_uf_alg<Graph, unite_early, find_atomic_halve, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      case 9:
        return run_multiple_uf_alg<Graph, unite_nd, find_compress, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 10:
        return run_multiple_uf_alg<Graph, unite_nd, find_naive, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 11:
        return run_multiple_uf_alg<Graph, unite_nd, find_atomic_split, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 12:
        return run_multiple_uf_alg<Graph, unite_nd, find_atomic_halve, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      case 13:
        return run_multiple_uf_alg<Graph, unite_rem_lock, find_compress, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 14:
        return run_multiple_uf_alg<Graph, unite_rem_lock, find_naive, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 15:
        return run_multiple_uf_alg<Graph, unite_rem_lock, find_atomic_split, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 16:
        return run_multiple_uf_alg<Graph, unite_rem_lock, find_atomic_halve, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      /* Jayanti strategies */
      case 17:
        return run_multiple_jayanti_alg<Graph, find_twotrysplit, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds,  P);
      case 18:
        return run_multiple_jayanti_alg<Graph, find_simple, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds,  P);

      /* UF Rem-CAS strategies */
      case 19:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_compress, split_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 20:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_compress, halve_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 21:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_compress, splice_simple, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 22:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_compress, splice_atomic, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      case 23:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_split, split_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 24:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_split, halve_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 25:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_split, splice_simple, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 26:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_split, splice_atomic, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      case 27:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_halve, split_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 28:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_halve, halve_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 29:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_halve, splice_simple, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 30:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_atomic_halve, splice_atomic, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      /* Liu-Tarjan algorithms */
      case 31:
        /* <parent_connect, update, shortcut> (Algorithm P) */
        return run_multiple_liu_tarjan_alg<Graph, parent_connect, simple_update, shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 32:
        /* <parent_connect, root_update, shortcut> (Algorithm R) */
        return run_multiple_liu_tarjan_alg<Graph, parent_connect, root_update, shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 33:
        /* <extended_connect, update, shortcut> (Algorithm E) */
        return run_multiple_liu_tarjan_alg<Graph, extended_connect, simple_update, shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 34:
        /* <parent_connect, update, full_shortcut> (PF) */
        return run_multiple_liu_tarjan_alg<Graph, parent_connect, simple_update, full_shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 35:
        /* <parent_connect, root_update, full_shortcut> (RF) */
        return run_multiple_liu_tarjan_alg<Graph, parent_connect, root_update, full_shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 36:
        /* <extended_connect, update, full_shortcut> (EF) */
        return run_multiple_liu_tarjan_alg<Graph, extended_connect, simple_update, full_shortcut, no_alter, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      /* Shiloach-Vishkin strategies */
      case 37:
        return run_multiple_shiloach_vishkin<Graph, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);

      /* UF Rem-CAS strategies */
      case 38:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_naive, split_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 39:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_naive, halve_atomic_one, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 40:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_naive, splice_simple, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);
      case 41:
        return run_multiple_uf_alg<Graph, unite_rem_cas, find_naive, splice_atomic, provides_initial_graph>(G, n, updates, batch_size, insert_to_query, rounds, P);




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

constexpr int num_tests = 42; // update if new algorithm is added



template <class Graph>
double Benchmark_nonempty_start(Graph& G, commandLine P) {
  using W = typename Graph::weight_type;
  int test_num = P.getOptionIntValue("-t", -1);
  int rounds = P.getOptionIntValue("-r", 5);
  assert(symmetric == true);


  double update_pct = P.getOptionDoubleValue("-update_pct", 0.1); /* percentage of edges to sample for updates */
  double deletion_pct = P.getOptionDoubleValue("-deletion_pct", 0.6) + update_pct; /* percentage of edges to delete */

//  double insert_to_query_ratio = P.getOptionDoubleValue("-insert_to_query", 0.5); /* 2 ins/query */
//  assert(insert_to_query_ratio < 1.0);
  size_t insert_to_query = P.getOptionLongValue("-insert_to_query", 2);


  /* idea: hash every edge to a double in (0, 1). */

  size_t batch_size = P.getOptionLongValue("-batch_size", 1000000); /* batch size */

  /* 1) sample edges to use as updates */

  auto hash_to_double = [&] (const uintE& u, const uintE& v) -> double {
    auto min_v = std::min(u, v);
    auto max_v = std::max(u, v);
    size_t hashed_v = pbbs::hash64((static_cast<size_t>(min_v) << 32UL) + static_cast<size_t>(max_v));
    return static_cast<double>(hashed_v) / static_cast<double>(std::numeric_limits<size_t>::max());
  };

  auto update_pred =  [&] (const uintE& u, const uintE& v, const W& wgh) {
    return hash_to_double(u,v) < update_pct; /* return in sample */
  };
  auto updates_arr = sample_edges(G, update_pred);
  auto updates = pbbs::sequence((std::tuple<uintE, uintE>*)updates_arr.E, updates_arr.m);
  updates_arr.E = nullptr; /* relinquish memory */

  /* 2) call filter_graph to delete all deletions + updates from G */
  auto delete_pred =  [&] (const uintE& u, const uintE& v, const W& wgh) {
    return hash_to_double(u,v) > deletion_pct; /* keep */
  };
  /* Can also try using filter_edges */
  auto FG = filter_graph(G, delete_pred);
  std::cout << "### Initial graph size in edges = " << FG.m << std::endl;

  /* select a static_alg */

  size_t n = FG.n;
  if (test_num == -1) {
    std::cout << "test_num == -1" << std::endl;
    for (int i=1; i <= num_tests; i++) {
      connectit::pick_test<decltype(FG), true>(FG, n, updates, i, batch_size, insert_to_query, rounds, P);
    }
  } else {
    connectit::pick_test<decltype(FG), true>(FG, n, updates, test_num, batch_size, insert_to_query, rounds, P);
  }

  return 1.0;
}

// Two different modes: one starting without a base graph, and one starting with
// a base graph.

#ifdef EMPTY_STARTING_GRAPH

/* run synthetic coo */
int main(int argc, char* argv[]) {
  auto P = commandLine(argc, argv, "");
  int test_num = P.getOptionIntValue("-t", -1);
  int rounds = P.getOptionIntValue("-r", 5);

  auto in_file = P.getOptionValue("-in_file", "");
  if (in_file == "") {
    std::cout << "must specify valid coo input: " << in_file << std::endl;
    abort();
  }

  sequence<char> chars = pbbs::char_seq_from_file(in_file);
  auto tokens = pbbslib::tokenize(chars, [] (const char c) { return pbbs::is_space(c); });
  // parseback to ints

  assert(tokens.size() % 2 == 0); // m tuples, into two tokens each

  if (tokens.size() % 2 != 0) {
    std::cout << "malformed coo input" << std::endl;
    std::cout.flush();
    abort();
  }
  size_t m = tokens.size() / 2;
  auto updates = pbbs::sequence<std::tuple<uintE, uintE>>(m);


  uintE n = 0;
  parallel_for(0, m, [&] (size_t i) {
    uintE l = std::atoi(tokens[2*i]);
    uintE r = std::atoi(tokens[2*i + 1]);
    if (l > n) {
      pbbs::write_min<uintE>(&n, l, std::greater<uintE>());
    }
    if (r > n) {
      pbbs::write_min<uintE>(&n, r, std::greater<uintE>());
    }
    updates[i] = std::make_tuple(l, r);
  });
  n = n + 1; /* 0 indexed */

  size_t batch_size = P.getOptionLongValue("-batch_size", 1000000); /* batch size */

  size_t insert_to_query = P.getOptionLongValue("-insert_to_query", 2);

  auto FG = edge_array<pbbs::empty>();
  if (test_num == -1) {
    for (int i=1; i <= num_tests; i++) {
      connectit::pick_test<decltype(FG), false>(FG, n, updates, i, batch_size, insert_to_query, rounds, P);
    }
  } else {
    connectit::pick_test<decltype(FG), false>(FG, n, updates, test_num, batch_size, insert_to_query, rounds, P);
  }
}

#else
generate_symmetric_once_main(Benchmark_nonempty_start, true);
#endif

