#include "ligra/edge_map_reduce.h"
#include "ligra/ligra.h"

#include <math.h>

template <class Graph>
void CoSimRank(Graph& G, uintE v, uintE u, double eps = 0.000001, double c = 0.6, size_t max_iters = 100) {
  using W = typename Graph::weight_type;
  const uintE n = G.n;

  double one_over_n = 1/(double)n;

  auto p_curr_v = pbbs::sequence<double>(n, static_cast<double>(0));
  p_curr_v[v] = static_cast<double>(1);
  auto p_next_v = pbbs::sequence<double>(n, static_cast<double>(0));
  auto frontier_v = pbbs::sequence<bool>(n, false);
  frontier_v[v] = true;

  auto p_curr_u = pbbs::sequence<double>(n, static_cast<double>(0));
  p_curr_v[u] = static_cast<double>(1);
  auto p_next_u = pbbs::sequence<double>(n, static_cast<double>(0));
  auto frontier_u = pbbs::sequence<bool>(n, false);
  frontier_u[u] = true;

  auto p_div_v = pbbs::sequence<double>(n, [&] (size_t i) -> double {
    return p_curr_v[i] / static_cast<double>(G.get_vertex(i).getOutDegree());
  });

  auto p_div_u = pbbs::sequence<double>(n, [&] (size_t i) -> double {
    return p_curr_u[i] / static_cast<double>(G.get_vertex(i).getOutDegree());
  });

  // read from special array of just degrees

  auto degrees = pbbs::sequence<uintE>(n, [&] (size_t i) { return G.get_vertex(i).getOutDegree(); });

  vertexSubset Frontier_v(n,n,frontier_v.to_array());
  auto EM_v = EdgeMap<double, Graph>(G, std::make_tuple(UINT_E_MAX, static_cast<double>(0)), (size_t)G.m/1000);

  vertexSubset Frontier_u(n,n,frontier_u.to_array());
  auto EM_u = EdgeMap<double, Graph>(G, std::make_tuple(UINT_E_MAX, static_cast<double>(0)), (size_t)G.m/1000);

  auto cond_f = [&] (const uintE& v) { return true; };
  auto map_f_v = [&] (const uintE& d, const uintE& s, const W& wgh) -> double {
    return p_div_v[s];
  };
  auto reduce_f = [&] (double l, double r) { return l + r; };
  auto apply_f_v = [&] (std::tuple<uintE, double> k) {
    const uintE& w = std::get<0>(k);
    const double& contribution = std::get<1>(k);
    p_next_v[w] = contribution;
    p_div_v[w] = p_next_v[w]/static_cast<double>(degrees[w]);
    return Maybe<std::tuple<uintE, double>>();
  };

  auto map_f_u = [&] (const uintE& d, const uintE& s, const W& wgh) -> double {
    return p_div_u[s];
  };
  auto apply_f_u = [&] (std::tuple<uintE, double> k) {
    const uintE& w = std::get<0>(k);
    const double& contribution = std::get<1>(k);
    p_next_u[w] = contribution;
    p_div_u[w] = p_next_u[w]/static_cast<double>(degrees[w]);
    return Maybe<std::tuple<uintE, double>>();
  };

  size_t iter = 0;
  double sim = 0;
  while (iter++ < max_iters) {
    timer t; t.start();
    // SpMV
    timer tt; tt.start();
    EM_v.template edgeMapReduce_dense<double, double>(Frontier_v, cond_f, map_f_v, reduce_f, apply_f_v, 0.0, no_output);
    EM_u.template edgeMapReduce_dense<double, double>(Frontier_u, cond_f, map_f_u, reduce_f, apply_f_u, 0.0, no_output);
    tt.stop(); tt.reportTotal("em time");

    sim += (pow(c, iter) * std::inner_product(p_next_u.begin(), p_next_u.end(), p_next_v.begin(), 0));

    // Check convergence: compute L1-norm between p_curr and p_next
    auto differences_v = pbbs::delayed_seq<double>(n, [&] (size_t i) {
      auto d = p_curr_v[i];
      p_curr_v[i] = 0;
      return fabs(d-p_next_v[i]);
    });
    double L1_norm_v = pbbs::reduce(differences_v, pbbs::addm<double>());

    auto differences_u = pbbs::delayed_seq<double>(n, [&] (size_t i) {
      auto d = p_curr_u[i];
      p_curr_u[i] = 0;
      return fabs(d-p_next_u[i]);
    });
    double L1_norm_u = pbbs::reduce(differences_u, pbbs::addm<double>());
    if(L1_norm_v < eps && L1_norm_u < eps) break;
    debug(cout << "L1_norm = " << L1_norm_v << ", " << L1_norm_u << endl;);

    // Reset p_curr
    std::swap(p_curr_v,p_next_v);
    std::swap(p_curr_u,p_next_u);
    t.stop(); t.reportTotal("iteration time");
  }
  Frontier.del();
  auto max_pr_v = pbbslib::reduce_max(p_next_v);
  auto max_pr_u = pbbslib::reduce_max(p_next_u);
  cout << "max_pr = " << max_pr_v << ", " << max_pr_u << endl;
  cout << "sim = " << sim << endl;
}
