#pragma once

#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"
#include "gbbs/graph.h"

#include "common.h"

namespace gbbs {
namespace spanning_forest {

template <class Seq>
inline size_t num_cc(Seq& labels) {
  size_t n = labels.size();
  auto flags =
      sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    if (!flags[labels[i]]) {
      flags[labels[i]] = 1;
    }
  });
  parlay::scan_inplace(flags);
  std::cout << "# n_cc = " << flags[n] << "\n";
  return flags[n];
}

template <class Seq>
inline void RelabelDet(Seq& ids) {
  using T = typename Seq::value_type;
  size_t n = ids.size();
  auto component_map = sequence<T>(n + 1, (T)0);
  T cur_comp = 0;
  for (size_t i = 0; i < n; i++) {
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
  for (size_t i = 0; i < correct.size(); i++) {
    if ((correct[i] != check[i])) {
      is_correct = false;
      std::cout << "# at i = " << i << " cor = " << correct[i]
                << " got: " << check[i] << std::endl;
      std::cout.flush();
      assert(correct[i] == check[i]);
      abort();
    }
    if (correct[i] > max_cor) {
      gbbs::write_max(&max_cor, correct[i], std::less<parent>());
    }
    if (check[i] > max_chk) {
      gbbs::write_max(&max_chk, check[i], std::less<parent>());
    }
  }
  std::cout << "# correctness check: " << is_correct << std::endl;
  std::cout << "# max_cor = " << max_cor << " max_chk = " << max_chk
            << std::endl;
}

inline sequence<std::tuple<uintE, uintE, gbbs::empty>> double_edges(
    sequence<edge>& in) {
  using weighted_edge = std::tuple<uintE, uintE, gbbs::empty>;
  auto double_in =
      sequence<weighted_edge>::from_function(in.size() * 2, [&](size_t i) {
        size_t ind = i / 2;
        auto[u, v] = in[ind];
        if (i % 2 == 0) {
          return std::make_tuple(u, v, gbbs::empty());
        } else {
          return std::make_tuple(v, u, gbbs::empty());
        }
      });
  return double_in;
}

inline void check_spanning_forest(size_t n, sequence<edge>& correct,
                                  sequence<edge>& check) {
  // check sizes
  if (correct.size() != check.size()) {
    std::cout << "## Correct forest has: " << correct.size()
              << " many edges, but supplied forest has: " << check.size()
              << " edges." << std::endl;
    exit(-1);
  }
  /* convert to graphs, and check connectivity induced by edges */

  auto double_correct = double_edges(correct);
  auto G_double = sym_graph_from_edges(double_correct, n);
  auto conn_correct = workefficient_cc::CC(G_double);
  num_cc(conn_correct);

  auto double_check = double_edges(check);
  auto G_check = sym_graph_from_edges(double_check, n);
  auto conn_check = workefficient_cc::CC(G_check);
  num_cc(conn_check);

  cc_check(conn_correct, conn_check);
}

}  // namespace spanning_forest
}  // namespace gbbs
