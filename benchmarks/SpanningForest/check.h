#pragma once

#include "ligra/graph.h"
#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"
//#include "benchmarks/Connectivity/UnionFind/Connectivity.h"
//#include "benchmarks/SpanningForest/SDB14/SpanningForest.h"
//#include "benchmarks/SpanningForest/BFSSF/SpanningForest.h"

#include "common.h"

namespace spanning_forest {

  template <class Seq>
  inline size_t num_cc(Seq& labels) {
    size_t n = labels.size();
    auto flags = sequence<uintE>(n + 1, [&](size_t i) { return 0; });
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      if (!flags[labels[i]]) {
        flags[labels[i]] = 1;
      }
    });
    pbbslib::scan_add_inplace(flags);
    std::cout << "# n_cc = " << flags[n] << "\n";
    return flags[n];
  }

  template <class Seq>
  inline void RelabelDet(Seq& ids) {
    using T = typename Seq::value_type;
    size_t n = ids.size();
    auto component_map = pbbs::sequence<T>(n + 1, (T)0);
    T cur_comp = 0;
    for (size_t i=0; i<n; i++) {
      T comp = ids[i];
      if (component_map[comp] == 0) {
        component_map[comp] = cur_comp;
        ids[i] = cur_comp;
        cur_comp++;
      } else {
        ids[i] = component_map[comp];
      }
    }
  }

  template <class S1, class S2>
  inline void cc_check(S1& correct, S2& check) {
    RelabelDet(correct);
    RelabelDet(check);

    bool is_correct = true;
    parent max_cor = 0;
    parent max_chk = 0;
    for (size_t i=0; i<correct.size(); i++) {
      if ((correct[i] != check[i])) {
        is_correct = false;
        std::cout << "# at i = " << i << " cor = " << correct[i] << " got: " << check[i] << std::endl;
        std::cout.flush();
        assert(correct[i] == check[i]);
        abort();
      }
      if (correct[i] > max_cor) {
        pbbs::write_max(&max_cor, correct[i], std::less<parent>());
      }
      if (check[i] > max_chk) {
        pbbs::write_max(&max_chk, check[i], std::less<parent>());
      }
    }
    cout << "# correctness check: " << is_correct << endl;
    cout << "# max_cor = " << max_cor << " max_chk = " << max_chk << endl;
  }

  pbbs::sequence<std::tuple<uintE, uintE, pbbs::empty>> double_edges(pbbs::sequence<edge>& in) {
    using weighted_edge = std::tuple<uintE, uintE, pbbs::empty>;
    auto double_in = pbbs::sequence<weighted_edge>(in.size() * 2, [&] (size_t i) {
      size_t ind = i/2;
      auto [u, v] = in[ind];
      if (i % 2 == 0) {
        return std::make_tuple(u, v, pbbs::empty());
      } else {
        return std::make_tuple(v, u, pbbs::empty());
      }
    });
    return double_in;
  }

  void check_spanning_forest(size_t n, pbbs::sequence<edge>& correct, pbbs::sequence<edge>& check) {
    // check sizes
    if (correct.size() != check.size()) {
      std::cout << "## Correct forest has: " << correct.size() << " many edges, but supplied forest has: " << check.size() << " edges." << std::endl;
      exit(-1);
    }
    /* convert to graphs, and check connectivity induced by edges */

    auto double_check = double_edges(check);
    auto G_check = sym_graph_from_edges(double_check, n);
    auto conn_check = workefficient_cc::CC(G_check);
    num_cc(conn_check);

    auto double_correct = double_edges(correct);
    auto G_double = sym_graph_from_edges(double_correct, n);
    auto conn_correct = workefficient_cc::CC(G_double);
    num_cc(conn_correct);

    cc_check(conn_correct, conn_check);
  }

}

