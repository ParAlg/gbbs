template <class Graph>
double t_bfs_cc(Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
  time(t, auto CC = bfs_cc::CC(G));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

template <class Graph>
double t_sample_only_algorithm_cc(Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
  time(t, auto CC = labelprop_cc::CC</*use_permutation = */false>(G););
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

namespace connectit {

  template<
    class Graph,
    SamplingOption sampling_option,
    UniteOption    unite_option,
    FindOption     find_option>
  bool run_multiple_uf_alg(
      Graph& G,
      size_t rounds,
      pbbs::sequence<parent>& correct,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
      timer tt; tt.start();
      auto CC =
          run_uf_alg<
            Graph,
            sampling_option,
            find_option,
            unite_option>(G, P);
      double t = tt.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, CC);
      }
      return t;
    };
    auto name = uf_options_to_string<sampling_option, find_option, unite_option>();
    return run_multiple(G, rounds, correct, name, P, test);
  }

  template<
    class Graph,
    SamplingOption sampling_option,
    UniteOption    unite_option,
    FindOption     find_option,
    SpliceOption   splice_option>
  bool run_multiple_uf_alg(
      Graph& G,
      size_t rounds,
      pbbs::sequence<parent>& correct,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
      timer tt; tt.start();
      auto CC =
          run_uf_alg<
            Graph,
            sampling_option,
            find_option,
            unite_option,
            splice_option>(G, P);
      double t = tt.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, CC);
      }
      return t;
    };
    auto name = uf_options_to_string<sampling_option, find_option, unite_option, splice_option>();
    return run_multiple(G, rounds, correct, name, P, test);
  }

  template<
    class Graph,
    SamplingOption    sampling_option,
    JayantiFindOption find_option>
  bool run_multiple_jayanti_alg(
      Graph& G,
      size_t rounds,
      pbbs::sequence<parent>& correct,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
      timer tt; tt.start();
      auto CC =
          run_jayanti_alg<
            Graph,
            sampling_option,
            find_option>(G, P);
      double t = tt.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, CC);
      }
      return t;
    };
    auto name = jayanti_options_to_string<sampling_option, find_option>();
    return run_multiple(G, rounds, correct, name, P, test);
  }


  template<
    class Graph,
    SamplingOption          sampling_option,
    LiuTarjanConnectOption  connect_option,
    LiuTarjanUpdateOption   update_option,
    LiuTarjanShortcutOption shortcut_option,
    LiuTarjanAlterOption    alter_option>
  bool run_multiple_liu_tarjan_alg(
      Graph& G,
      size_t rounds,
      pbbs::sequence<parent>& correct,
      commandLine& P) {
    auto test = [&] (Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
      timer tt; tt.start();
      auto CC =
          run_liu_tarjan_alg<
            Graph,
            sampling_option,
            connect_option,
            update_option,
            shortcut_option,
            alter_option>(G, P);
      double t = tt.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, CC);
      }
      return t;
    };
    auto name = liu_tarjan_options_to_string<sampling_option, connect_option, update_option, shortcut_option, alter_option>();
    return run_multiple(G, rounds, correct, name, P, test);
  }

  template<
    class Graph,
    SamplingOption    sampling_option,
    template <class G> class Algorithm,
    AlgorithmType algorithm_type>
  bool run_multiple_sample_only_alg(
      Graph& G,
      size_t rounds,
      pbbs::sequence<parent>& correct,
      commandLine& P,
      std::string name) {
    auto test = [&] (Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
      timer tt; tt.start();
      auto CC = run_sample_only_alg<Graph, sampling_option, Algorithm, algorithm_type>(G, P);
      double t = tt.stop();
      if (P.getOptionValue("-check")) {
        cc_check(correct, CC);
      }
      return t;
    };
    auto test_name = name + "; " + sampling_to_string<sampling_option>();
    return run_multiple(G, rounds, correct, test_name, P, test);
  }


  template <class Graph>
  bool pick_test(Graph& G, size_t id, size_t rounds, commandLine P, pbbs::sequence<parent>& correct) {

    switch (id) {
      case 0:
        return run_multiple(G, rounds, correct, "gbbs_cc", P, t_gbbs_cc<Graph>);

//      /* Union Find Algorithms: {4 sampling} x {unite, unite_early, unite_nd} x {find_compress, find_naive, find_atomic_split, find_atomic_halve} = 48 */
//      /* + {4 sampling} x {unite_rem_cas} x {find_naive, find_atomic_split, find_atomic_halve}  = 12 */
//      case 1:
//        return run_multiple_uf_alg<Graph, no_sampling, unite, find_compress>(G, rounds, correct, P);
//      case 2:
//        return run_multiple_uf_alg<Graph, no_sampling, unite, find_naive>(G, 1, correct, P);
//      case 3:
//        return run_multiple_uf_alg<Graph, no_sampling, unite, find_atomic_split>(G, rounds, correct, P);
//      case 4:
//        return run_multiple_uf_alg<Graph, no_sampling, unite, find_atomic_halve>(G, rounds, correct, P);
//
//      case 5:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_early, find_compress>(G, rounds, correct, P);
//      case 6:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_early, find_naive>(G, 1, correct, P);
//      case 7:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_early, find_atomic_split>(G, rounds, correct, P);
//      case 8:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_early, find_atomic_halve>(G, rounds, correct, P);
//
//      case 9:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_nd, find_compress>(G, rounds, correct, P);
//      case 10:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_nd, find_naive>(G, 1, correct, P);
//      case 11:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_nd, find_atomic_split>(G, rounds, correct, P);
//      case 12:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_nd, find_atomic_halve>(G, rounds, correct, P);
//
//      case 13:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_lock, find_naive>(G, 1, correct, P);
//      case 14:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_lock, find_atomic_split>(G, rounds, correct, P);
//      case 15:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_lock, find_atomic_halve>(G, rounds, correct, P);
//
//      case 16:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite, find_compress>(G, rounds, correct, P);
//      case 17:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite, find_naive>(G, 1, correct, P);
//      case 18:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite, find_atomic_split>(G, rounds, correct, P);
//      case 19:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite, find_atomic_halve>(G, rounds, correct, P);
//
//      case 20:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_early, find_compress>(G, rounds, correct, P);
//      case 21:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_early, find_naive>(G, 1, correct, P);
//      case 22:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_early, find_atomic_split>(G, rounds, correct, P);
//      case 23:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_early, find_atomic_halve>(G, rounds, correct, P);
//
//      case 24:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_nd, find_compress>(G, rounds, correct, P);
//      case 25:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_nd, find_naive>(G, 1, correct, P);
//      case 26:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_nd, find_atomic_split>(G, rounds, correct, P);
//      case 27:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_nd, find_atomic_halve>(G, rounds, correct, P);
//
//      case 28:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_lock, find_naive>(G, 1, correct, P);
//      case 29:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_lock, find_atomic_split>(G, rounds, correct, P);
//      case 30:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_lock, find_atomic_halve>(G, rounds, correct, P);
//
//      case 31:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite, find_compress>(G, rounds, correct, P);
//      case 32:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite, find_naive>(G, 1, correct, P);
//      case 33:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite, find_atomic_split>(G, rounds, correct, P);
//      case 34:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite, find_atomic_halve>(G, rounds, correct, P);
//
//      case 35:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_early, find_compress>(G, rounds, correct, P);
//      case 36:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_early, find_naive>(G, 1, correct, P);
//      case 37:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_early, find_atomic_split>(G, rounds, correct, P);
//      case 38:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_early, find_atomic_halve>(G, rounds, correct, P);
//
//      case 39:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_nd, find_compress>(G, rounds, correct, P);
//      case 40:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_nd, find_naive>(G, 1, correct, P);
//      case 41:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_nd, find_atomic_split>(G, rounds, correct, P);
//      case 42:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_nd, find_atomic_halve>(G, rounds, correct, P);
//
//      case 43:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_lock, find_naive>(G, 1, correct, P);
//      case 44:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_lock, find_atomic_split>(G, rounds, correct, P);
//      case 45:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_lock, find_atomic_halve>(G, rounds, correct, P);
//
//

      case 46:
        return run_multiple_uf_alg<Graph, sample_kout, unite, find_compress>(G, rounds, correct, P);
//      case 47:
//        return run_multiple_uf_alg<Graph, sample_kout, unite, find_naive>(G, 1, correct, P);
//      case 48:
//        return run_multiple_uf_alg<Graph, sample_kout, unite, find_atomic_split>(G, rounds, correct, P);
//      case 49:
//        return run_multiple_uf_alg<Graph, sample_kout, unite, find_atomic_halve>(G, rounds, correct, P);
//
//      case 50:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_early, find_compress>(G, rounds, correct, P);
//      case 51:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_early, find_naive>(G, 1, correct, P);
//      case 52:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_early, find_atomic_split>(G, rounds, correct, P);
//      case 53:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_early, find_atomic_halve>(G, rounds, correct, P);
//
//      case 54:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_nd, find_compress>(G, rounds, correct, P);
//      case 55:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_nd, find_naive>(G, 1, correct, P);
//      case 56:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_nd, find_atomic_split>(G, rounds, correct, P);
//      case 57:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_nd, find_atomic_halve>(G, rounds, correct, P);
//
//      case 58:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_lock, find_naive>(G, 1, correct, P);
//      case 59:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_lock, find_atomic_split>(G, rounds, correct, P);
//      case 60:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_lock, find_atomic_halve>(G, rounds, correct, P);
//
//      /* Jayanti strategies {4 sampling} x {find_twotrysplit, find_simple} */
//      case 61:
//        return run_multiple_jayanti_alg<Graph, sample_kout, find_twotrysplit>(G, rounds, correct,  P);
//      case 62:
//        return run_multiple_jayanti_alg<Graph, sample_kout, find_simple>(G, rounds, correct,  P);
//      case 63:
//        return run_multiple_jayanti_alg<Graph, sample_bfs, find_twotrysplit>(G, rounds, correct,  P);
//      case 64:
//        return run_multiple_jayanti_alg<Graph, sample_bfs, find_simple>(G, rounds, correct,  P);
//      case 65:
//        return run_multiple_jayanti_alg<Graph, sample_ldd, find_twotrysplit>(G, rounds, correct,  P);
//      case 66:
//        return run_multiple_jayanti_alg<Graph, sample_ldd, find_simple>(G, rounds, correct,  P);
//      case 67:
//        return run_multiple_jayanti_alg<Graph, no_sampling, find_twotrysplit>(G, rounds, correct,  P);
//      case 68:
//        return run_multiple_jayanti_alg<Graph, no_sampling, find_simple>(G, rounds, correct,  P);
//
//      /* Label Propagation strategies: {4 sampling}*/
//      case 69:
//        return run_multiple_sample_only_alg<Graph, sample_kout, labelprop_cc::LPAlgorithm, label_prop_type>(G, rounds, correct, P, "label_prop");
//      case 70:
//        return run_multiple_sample_only_alg<Graph, sample_bfs, labelprop_cc::LPAlgorithm, label_prop_type>(G, rounds, correct, P, "label_prop");
//      case 71:
//        return run_multiple_sample_only_alg<Graph, sample_ldd, labelprop_cc::LPAlgorithm, label_prop_type>(G, rounds, correct, P, "label_prop");
//      case 72:
//        return run_multiple_sample_only_alg<Graph, no_sampling, labelprop_cc::LPAlgorithm, label_prop_type>(G, rounds, correct, P, "label_prop");
//
//      /* UF Rem-CAS strategies: {4 sampling} x {find_atomic_split, find_atomic_halve, find_naive} x {split_atomic_one, halve_atomic_one, splice_simple, splice_atomic} */
//      case 73:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_split, split_atomic_one>(G, rounds, correct, P);
//      case 74:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_split, halve_atomic_one>(G, rounds, correct, P);
//      case 75:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_split, splice_simple>(G, rounds, correct, P);
//      case 76:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_split, splice_atomic>(G, rounds, correct, P);
//
//      case 77:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_halve, split_atomic_one>(G, rounds, correct, P);
//      case 78:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_halve, halve_atomic_one>(G, rounds, correct, P);
//      case 79:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_halve, splice_simple>(G, rounds, correct, P);
//      case 80:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_atomic_halve, splice_atomic>(G, rounds, correct, P);
//
//      case 81:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_naive, split_atomic_one>(G, rounds, correct, P);
//      case 82:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_naive, halve_atomic_one>(G, rounds, correct, P);
//      case 83:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_naive, splice_simple>(G, rounds, correct, P);
//      case 84:
//        return run_multiple_uf_alg<Graph, no_sampling, unite_rem_cas, find_naive, splice_atomic>(G, rounds, correct, P);
//
//      case 85:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_split, split_atomic_one>(G, rounds, correct, P);
//
//      case 86:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_split, halve_atomic_one>(G, rounds, correct, P);
//      case 87:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_split, splice_simple>(G, rounds, correct, P);
//      case 88:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_split, splice_atomic>(G, rounds, correct, P);
//
//      case 89:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_halve, split_atomic_one>(G, rounds, correct, P);
//      case 90:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_halve, halve_atomic_one>(G, rounds, correct, P);
//      case 91:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_halve, splice_simple>(G, rounds, correct, P);
//      case 92:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_atomic_halve, splice_atomic>(G, rounds, correct, P);
//
//      case 93:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_naive, split_atomic_one>(G, rounds, correct, P);
//      case 94:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_naive, halve_atomic_one>(G, rounds, correct, P);
//      case 95:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_naive, splice_simple>(G, rounds, correct, P);
//      case 96:
//        return run_multiple_uf_alg<Graph, sample_kout, unite_rem_cas, find_naive, splice_atomic>(G, rounds, correct, P);
//
//      case 97:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_split, split_atomic_one>(G, rounds, correct, P);
//      case 98:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_split, halve_atomic_one>(G, rounds, correct, P);
//      case 99:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_split, splice_simple>(G, rounds, correct, P);
//      case 100:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_split, splice_atomic>(G, rounds, correct, P);
//
//      case 101:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_halve, split_atomic_one>(G, rounds, correct, P);
//      case 102:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_halve, halve_atomic_one>(G, rounds, correct, P);
//      case 103:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_halve, splice_simple>(G, rounds, correct, P);
//      case 104:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_atomic_halve, splice_atomic>(G, rounds, correct, P);
//
//      case 105:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_naive, split_atomic_one>(G, rounds, correct, P);
//      case 106:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_naive, halve_atomic_one>(G, rounds, correct, P);
//      case 107:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_naive, splice_simple>(G, rounds, correct, P);
//      case 108:
//        return run_multiple_uf_alg<Graph, sample_bfs, unite_rem_cas, find_naive, splice_atomic>(G, rounds, correct, P);
//
//      case 109:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_split, split_atomic_one>(G, rounds, correct, P);
//      case 110:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_split, halve_atomic_one>(G, rounds, correct, P);
//      case 111:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_split, splice_simple>(G, rounds, correct, P);
//      case 112:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_split, splice_atomic>(G, rounds, correct, P);
//
//      case 113:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_halve, split_atomic_one>(G, rounds, correct, P);
//      case 114:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_halve, halve_atomic_one>(G, rounds, correct, P);
//      case 115:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_halve, splice_simple>(G, rounds, correct, P);
//      case 116:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_atomic_halve, splice_atomic>(G, rounds, correct, P);
//
//      case 117:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_naive, split_atomic_one>(G, rounds, correct, P);
//      case 118:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_naive, halve_atomic_one>(G, rounds, correct, P);
//      case 119:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_naive, splice_simple>(G, rounds, correct, P);
//      case 120:
//        return run_multiple_uf_alg<Graph, sample_ldd, unite_rem_cas, find_naive, splice_atomic>(G, rounds, correct, P);
//
//
//      /* Liu-Tarjan algorithms */
//      case 121:
//        /* <parent_connect, update, shortcut> (Algorithm P) */
//        return run_multiple_liu_tarjan_alg<Graph, no_sampling, parent_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
//      case 122:
//        /* <parent_connect, root_update, shortcut> (Algorithm R) */
//        return run_multiple_liu_tarjan_alg<Graph, no_sampling, parent_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
//      case 123:
//        /* <extended_connect, update, shortcut> (Algorithm E) */
//        return run_multiple_liu_tarjan_alg<Graph, no_sampling, extended_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
//      case 124:
//        /* <parent_connect, update, full_shortcut> (PF) */
//        return run_multiple_liu_tarjan_alg<Graph, no_sampling, parent_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
//      case 125:
//        /* <parent_connect, root_update, full_shortcut> (RF) */
//        return run_multiple_liu_tarjan_alg<Graph, no_sampling, parent_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);
//      case 126:
//        /* <extended_connect, update, full_shortcut> (EF) */
//        return run_multiple_liu_tarjan_alg<Graph, no_sampling, extended_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
//
//      case 127:
//        /* <parent_connect, update, shortcut> (Algorithm P) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
//      case 128:
//        /* <parent_connect, root_update, shortcut> (Algorithm R) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
//      case 129:
//        /* <extended_connect, update, shortcut> (Algorithm E) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_kout, extended_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
//      case 130:
//        /* <parent_connect, update, full_shortcut> (PF) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
//      case 131:
//        /* <parent_connect, root_update, full_shortcut> (RF) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_kout, parent_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);
//      case 132:
//        /* <extended_connect, update, full_shortcut> (EF) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_kout, extended_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
//
//      case 133:
//        /* <parent_connect, update, shortcut> (Algorithm P) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, parent_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
//      case 134:
//        /* <parent_connect, root_update, shortcut> (Algorithm R) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, parent_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
//      case 135:
//        /* <extended_connect, update, shortcut> (Algorithm E) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, extended_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
//      case 136:
//        /* <parent_connect, update, full_shortcut> (PF) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, parent_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
//      case 137:
//        /* <parent_connect, root_update, full_shortcut> (RF) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, parent_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);
//      case 138:
//        /* <extended_connect, update, full_shortcut> (EF) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, extended_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
//
//      case 139:
//        /* <parent_connect, update, shortcut> (Algorithm P) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, parent_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
//      case 140:
//        /* <parent_connect, root_update, shortcut> (Algorithm R) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, parent_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
//      case 141:
//        /* <extended_connect, update, shortcut> (Algorithm E) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, extended_connect, simple_update, shortcut, no_alter>(G, 1, correct, P);
//      case 142:
//        /* <parent_connect, update, full_shortcut> (PF) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, parent_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
//      case 143:
//        /* <parent_connect, root_update, full_shortcut> (RF) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, parent_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);
//      case 144:
//        /* <extended_connect, update, full_shortcut> (EF) */
//        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, extended_connect, simple_update, full_shortcut, no_alter>(G, 1, correct, P);
//
//      case 145:
//        return run_multiple(G, rounds, correct, "bfs_cc", P, t_bfs_cc<Graph>);
//
//      /* Shiloach-Vishkin strategies */
//      case 146:
//        return run_multiple_sample_only_alg<Graph, sample_kout, shiloachvishkin_cc::SVAlgorithm, shiloach_vishkin_type>(G, rounds, correct, P, "shiloach_vishkin");
//      case 147:
//        return run_multiple_sample_only_alg<Graph, sample_bfs, shiloachvishkin_cc::SVAlgorithm, shiloach_vishkin_type>(G, rounds, correct, P, "shiloach_vishkin");
//      case 148:
//        return run_multiple_sample_only_alg<Graph, sample_ldd, shiloachvishkin_cc::SVAlgorithm, shiloach_vishkin_type>(G, rounds, correct, P, "shiloach_vishkin");
//      case 149:
//        return run_multiple_sample_only_alg<Graph, no_sampling, shiloachvishkin_cc::SVAlgorithm, shiloach_vishkin_type>(G, rounds, correct, P, "shiloach_vishkin");
//
//      case 150:
//        return run_multiple_uf_alg<Graph, no_sampling, unite, find_split>(G, rounds, correct, P);
//
//
//      case 151:
//        return run_multiple_liu_tarjan_alg<Graph, no_sampling, extended_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
//      case 152:
//        return run_multiple_liu_tarjan_alg<Graph, no_sampling, extended_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);
//      case 153:
//        return run_multiple_liu_tarjan_alg<Graph, sample_kout, extended_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
//      case 154:
//        return run_multiple_liu_tarjan_alg<Graph, sample_kout, extended_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);
//      case 155:
//        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, extended_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
//      case 156:
//        return run_multiple_liu_tarjan_alg<Graph, sample_bfs, extended_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);
//      case 157:
//        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, extended_connect, root_update, shortcut, no_alter>(G, 1, correct, P);
//      case 158:
//        return run_multiple_liu_tarjan_alg<Graph, sample_ldd, extended_connect, root_update, full_shortcut, no_alter>(G, 1, correct, P);

      default:
      std::cout << "Unknown test" << std::endl;
      return false;
    }
  }
} // namespace connectit


