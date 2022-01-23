#pragma once

#include <math.h>
#include <numeric>

#include "gbbs/gbbs.h"

namespace gbbs {
template <class Graph>
struct Co_PR_F {
  using W = typename Graph::weight_type;
  double *p_curr, *p_next;
  Graph& G;
  Co_PR_F(double* _p_curr, double* _p_next, Graph& G)
      : p_curr(_p_curr), p_next(_p_next), G(G) {}
  inline bool update(
      const uintE& s, const uintE& d,
      const W& wgh) {  // update function applies PageRank equation
    p_next[d] += p_curr[s] / G.get_vertex(s).out_degree();
    return 1;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d,
                           const W& wgh) {  // atomic Update
    gbbs::fetch_and_add(&p_next[d], p_curr[s] / G.get_vertex(s).out_degree());
    return 1;
  }
  inline bool cond(intT d) { return cond_true(d); }
};

double inner_product(sequence<double>& arr1, sequence<double>& arr2) {
  auto prod = parlay::delayed_seq<double>(
      arr1.size(), [&](size_t i) { return arr1[i] * arr2[i]; });
  auto sum = parlay::reduce(prod);
  return sum;
}

template <class Graph>
void CoSimRank_edgeMap(Graph& G, uintE v, uintE u, double eps = 0.000001,
                       double c = 0.85, size_t max_iters = 100) {
  const uintE n = G.n;

  auto p_curr_v = sequence<double>(n, static_cast<double>(0));
  p_curr_v[v] = static_cast<double>(1);
  auto p_next_v = sequence<double>(n, static_cast<double>(0));
  auto frontier = sequence<bool>(n, true);

  auto p_curr_u = sequence<double>(n, static_cast<double>(0));
  p_curr_u[u] = static_cast<double>(1);
  auto p_next_u = sequence<double>(n, static_cast<double>(0));

  // read from special array of just degrees

  auto degrees = sequence<uintE>::from_function(
      n, [&](size_t i) { return G.get_vertex(i).out_degree(); });

  auto frontier_v = sequence<bool>(n, false);
  frontier_v[v] = true;
  auto frontier_u = sequence<bool>(n, false);
  frontier_u[u] = true;

  vertexSubset Frontier_u(n, n, std::move(frontier_u));
  vertexSubset Frontier_v(n, n, std::move(frontier_v));

  size_t iter = 0;
  double sim = u == v ? 1 : 0;
  while (iter++ < max_iters) {
    debug(timer t; t.start(););
    // SpMV
    auto Frontier_v_new =
        edgeMap(G, Frontier_v,
                Co_PR_F<Graph>(p_curr_v.begin(), p_next_v.begin(), G), 0);
    auto Frontier_u_new = edgeMap(
        G, Frontier_u, Co_PR_F<Graph>(p_curr_u.begin(), p_next_u.begin(), G),
        0);  //, no_output

    sim += ((double)pow(c, iter) * inner_product(p_next_u, p_next_v));

    Frontier_v = std::move(Frontier_v_new);
    Frontier_u = std::move(Frontier_u_new);

    // Check convergence: compute L1-norm between p_curr and p_next
    auto differences_v = parlay::delayed_seq<double>(
        n, [&](size_t i) { return fabs(p_curr_v[i] - p_next_v[i]); });
    double L1_norm_v = parlay::reduce(differences_v, parlay::addm<double>());

    auto differences_u = parlay::delayed_seq<double>(
        n, [&](size_t i) { return fabs(p_curr_u[i] - p_next_u[i]); });
    double L1_norm_u = parlay::reduce(differences_u, parlay::addm<double>());
    if (L1_norm_v < eps && L1_norm_u < eps) break;

    debug(std::cout << "L1_norm = " << L1_norm_v << ", " << L1_norm_u
                    << std::endl;);
    // Reset p_curr
    parallel_for(0, n, [&](size_t i) { p_curr_v[i] = static_cast<double>(0); });
    std::swap(p_curr_v, p_next_v);

    parallel_for(0, n, [&](size_t i) { p_curr_u[i] = static_cast<double>(0); });
    std::swap(p_curr_u, p_next_u);

    debug(t.stop(); t.next("iteration time"););
  }

  auto max_pr_v = parlay::reduce_max(p_next_v);
  auto max_pr_u = parlay::reduce_max(p_next_u);

  std::cout << "max_pr = " << max_pr_v << ", " << max_pr_u << std::endl;
  std::cout << "sim = " << sim << std::endl;
}

template <class Graph>
void CoSimRank(Graph& G, uintE v, uintE u, double eps = 0.000001,
               double c = 0.85, size_t max_iters = 100) {
  using W = typename Graph::weight_type;
  const uintE n = G.n;

  auto p_curr_v = sequence<double>(n, static_cast<double>(0));
  p_curr_v[v] = static_cast<double>(1);
  auto p_next_v = sequence<double>(n, static_cast<double>(0));
  auto frontier_v = sequence<std::tuple<bool, double>>::uninitialized(n);
  parallel_for(0, n,
               [&](size_t i) {
                 frontier_v[i] = std::make_tuple(false, static_cast<double>(0));
               },
               1000);
  frontier_v[v] = std::make_tuple(true, static_cast<double>(0));

  auto p_curr_u = sequence<double>(n, static_cast<double>(0));
  p_curr_u[u] = static_cast<double>(1);
  auto p_next_u = sequence<double>(n, static_cast<double>(0));
  auto frontier_u = sequence<std::tuple<bool, double>>::uninitialized(n);
  parallel_for(0, n,
               [&](size_t i) {
                 frontier_u[i] = std::make_tuple(false, static_cast<double>(0));
               },
               1000);
  frontier_u[u] = std::make_tuple(true, static_cast<double>(0));

  vertexSubsetData<double> Frontier_u(n, n, std::move(frontier_u));
  vertexSubsetData<double> Frontier_v(n, n, std::move(frontier_v));

  // read from special array of just degrees

  auto degrees = sequence<uintE>::from_function(
      n, [&](size_t i) { return G.get_vertex(i).out_degree(); });

  auto EM_v = EdgeMap<double, Graph>(
      G, std::make_tuple(UINT_E_MAX, static_cast<double>(0)),
      (size_t)G.m / 1000);
  auto EM_u = EdgeMap<double, Graph>(
      G, std::make_tuple(UINT_E_MAX, static_cast<double>(0)),
      (size_t)G.m / 1000);

  auto cond_f = [&](const uintE& x) { return true; };
  auto map_f_v = [&](const uintE& d, const uintE& s, const W& wgh) -> double {
    return p_curr_v[s] / static_cast<double>(G.get_vertex(s).out_degree());
  };
  auto reduce_f = [&](double l, double r) { return l + r; };
  auto apply_f_v = [&](std::tuple<uintE, double> k) {
    const uintE& w = std::get<0>(k);
    p_next_v[w] = std::get<1>(k);
    return std::optional<std::tuple<uintE, double>>(k);
  };

  auto map_f_u = [&](const uintE& d, const uintE& s, const W& wgh) -> double {
    return p_curr_u[s] / static_cast<double>(G.get_vertex(s).out_degree());
  };
  auto apply_f_u = [&](std::tuple<uintE, double> k) {
    const uintE& w = std::get<0>(k);
    p_next_u[w] = std::get<1>(k);
    return std::optional<std::tuple<uintE, double>>(k);
  };

  size_t iter = 0;
  double sim = u == v ? 1 : 0;
  while (iter++ < max_iters) {
    timer t;
    t.start();
    // SpMV
    timer tt;
    tt.start();
    auto Frontier_v_new = EM_v.template edgeMapReduce_dense<double, double>(
        Frontier_v, cond_f, map_f_v, reduce_f, apply_f_v, 0.0, 0);
    auto Frontier_u_new = EM_u.template edgeMapReduce_dense<double, double>(
        Frontier_u, cond_f, map_f_u, reduce_f, apply_f_u, 0.0,
        0);     //, no_output
    tt.stop();  // tt.next("em time");

    sim += ((double)pow(c, iter) * inner_product(p_next_u, p_next_v));

    Frontier_v = std::move(Frontier_v_new);
    Frontier_u = std::move(Frontier_u_new);

    // Check convergence: compute L1-norm between p_curr and p_next
    auto differences_v = parlay::delayed_seq<double>(n, [&](size_t i) {
      auto x = p_curr_v[i];
      p_curr_v[i] = 0;
      return fabs(x - p_next_v[i]);
    });
    double L1_norm_v = parlay::reduce(differences_v, parlay::addm<double>());

    auto differences_u = parlay::delayed_seq<double>(n, [&](size_t i) {
      auto x = p_curr_u[i];
      p_curr_u[i] = 0;
      return fabs(x - p_next_u[i]);
    });
    double L1_norm_u = parlay::reduce(differences_u, parlay::addm<double>());
    if (L1_norm_v < eps && L1_norm_u < eps) break;

    // Reset p_curr
    std::swap(p_curr_v, p_next_v);
    std::swap(p_curr_u, p_next_u);
    t.stop();  // t.next("iteration time");
  }
  auto max_pr_v = parlay::reduce_max(p_next_v);
  auto max_pr_u = parlay::reduce_max(p_next_u);

  std::cout << "max_pr = " << max_pr_v << ", " << max_pr_u << std::endl;
  std::cout << "sim = " << sim << std::endl;
}
}  // namespace gbbs
