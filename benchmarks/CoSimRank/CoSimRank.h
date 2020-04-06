#include "ligra/edge_map_reduce.h"
#include "ligra/ligra.h"

#include <math.h>
#include <numeric>


template <class Graph>
struct PR_F {
  using W = typename Graph::weight_type;
  double* p_curr, *p_next;
  Graph& G;
  PR_F(double* _p_curr, double* _p_next, Graph& G) :
    p_curr(_p_curr), p_next(_p_next), G(G) {}
  inline bool update(const uintE& s, const uintE& d, const W& wgh){ //update function applies PageRank equation
    p_next[d] += p_curr[s]/G.get_vertex(s).getOutDegree();
    return 1;
  }
  inline bool updateAtomic (const uintE& s, const uintE& d, const W& wgh) { //atomic Update
    pbbs::fetch_and_add(&p_next[d],p_curr[s]/G.get_vertex(s).getOutDegree());
    return 1;
  }
  inline bool cond (intT d) { return cond_true(d); }
};

template<class T>
T inner_product(T* arr1, T* arr2, size_t len, T* tmp=nullptr) {
  pbbs::sequence<T> prod;
  if (tmp == nullptr) prod = pbbs::sequence<T>(len, static_cast<T>(0));
  else prod = pbbs::sequence<T>(tmp, len);

  parallel_for(0, len, [&] (size_t i) { prod[i] = arr1[i] * arr2[i]; });
  auto sum = pbbslib::reduce_add(prod);

  if (tmp != nullptr) prod.to_array();
  return sum;
}


template <class Graph>
void CoSimRank_edgeMap(Graph& G, uintE v, uintE u, double eps = 0.000001, double c = 0.9, size_t max_iters = 100) {
  const uintE n = G.n;

  auto p_curr_v = pbbs::sequence<double>(n, static_cast<double>(0));
  p_curr_v[v] = static_cast<double>(1);
  auto p_next_v = pbbs::sequence<double>(n, static_cast<double>(0));
  auto frontier = pbbs::sequence<bool>(n, true);

  auto p_curr_u = pbbs::sequence<double>(n, static_cast<double>(0));
  p_curr_u[u] = static_cast<double>(1);
  auto p_next_u = pbbs::sequence<double>(n, static_cast<double>(0));

  // read from special array of just degrees

  auto degrees = pbbs::sequence<uintE>(n, [&] (size_t i) { return G.get_vertex(i).getOutDegree(); });

  auto frontier_v = pbbs::sequence<bool>(n, false);
  frontier_v[v] = true;
  auto frontier_u = pbbs::sequence<bool>(n, false);
  frontier_u[u] = true;

  vertexSubset Frontier_u(n,n,frontier_u.to_array());
  vertexSubset Frontier_v(n,n,frontier_v.to_array());

  size_t iter = 0;
  double sim = u == v ? 1 : 0;
  while (iter++ < max_iters) {
    debug(timer t; t.start(););
    // SpMV
    auto Frontier_v_new = edgeMap(G,Frontier_v,PR_F<Graph>(p_curr_v.begin(),p_next_v.begin(),G), 0, no_output);
    auto Frontier_u_new = edgeMap(G,Frontier_u,PR_F<Graph>(p_curr_u.begin(),p_next_u.begin(),G), 0, no_output);

    sim += ((double) pow(c, iter) * inner_product<double>(p_next_u.begin(), p_next_v.begin(), n));

    Frontier_u.del(); Frontier_v.del();
    Frontier_v = Frontier_v_new;
    Frontier_u = Frontier_u_new;

    // Check convergence: compute L1-norm between p_curr and p_next
    auto differences_v = pbbs::delayed_seq<double>(n, [&] (size_t i) {
      return fabs(p_curr_v[i]-p_next_v[i]);
    });
    double L1_norm_v = pbbs::reduce(differences_v, pbbs::addm<double>());

    auto differences_u = pbbs::delayed_seq<double>(n, [&] (size_t i) {
      return fabs(p_curr_u[i]-p_next_u[i]);
    });
    double L1_norm_u = pbbs::reduce(differences_u, pbbs::addm<double>());
    if(L1_norm_v < eps && L1_norm_u < eps) break;

    debug(cout << "L1_norm = " << L1_norm_v << ", " << L1_norm_u << endl;);
    // Reset p_curr
    parallel_for(0, n, [&] (size_t i) { p_curr_v[i] = static_cast<double>(0); });
    std::swap(p_curr_v,p_next_v);

    parallel_for(0, n, [&] (size_t i) { p_curr_u[i] = static_cast<double>(0); });
    std::swap(p_curr_u,p_next_u);

    debug(t.stop(); t.reportTotal("iteration time"););
  }
  Frontier_u.del(); Frontier_v.del();

  auto max_pr_v = pbbslib::reduce_max(p_next_v);
  auto max_pr_u = pbbslib::reduce_max(p_next_u);

  cout << "max_pr = " << max_pr_v << ", " << max_pr_u << endl;
  cout << "sim = " << sim << endl;
}

template <class Graph>
void CoSimRank(Graph& G, uintE v, uintE u, double eps = 0.000001, double c = 0.9, size_t max_iters = 100) {
  using W = typename Graph::weight_type;
  const uintE n = G.n;

  auto p_curr_v = pbbs::sequence<double>(n, static_cast<double>(0));
  p_curr_v[v] = static_cast<double>(1);
  auto p_next_v = pbbs::sequence<double>(n, static_cast<double>(0));
  auto frontier_v = pbbs::sequence<std::tuple<bool, double>>::no_init(n);
  parallel_for(0, n, [&](size_t i) { frontier_v[i] = std::make_tuple(false, static_cast<double>(0)); },
                 1000);
  frontier_v[v] = std::make_tuple(true, static_cast<double>(0));

  auto p_curr_u = pbbs::sequence<double>(n, static_cast<double>(0));
  p_curr_u[u] = static_cast<double>(1);
  auto p_next_u = pbbs::sequence<double>(n, static_cast<double>(0));
  auto frontier_u = pbbs::sequence<std::tuple<bool, double>>::no_init(n);
  parallel_for(0, n, [&](size_t i) { frontier_u[i] = std::make_tuple(false, static_cast<double>(0)); },
                 1000);
  frontier_u[u] = std::make_tuple(true, static_cast<double>(0));

  vertexSubsetData<double> Frontier_u(n,n,frontier_u.to_array());
  vertexSubsetData<double> Frontier_v(n,n,frontier_v.to_array());

  // read from special array of just degrees

  auto degrees = pbbs::sequence<uintE>(n, [&] (size_t i) { return G.get_vertex(i).getOutDegree(); });

  auto EM_v = EdgeMap<double, Graph>(G, std::make_tuple(UINT_E_MAX, static_cast<double>(0)), (size_t)G.m/1000);
  auto EM_u = EdgeMap<double, Graph>(G, std::make_tuple(UINT_E_MAX, static_cast<double>(0)), (size_t)G.m/1000);

  auto cond_f = [&] (const uintE& x) { return true; };
  auto map_f_v = [&] (const uintE& d, const uintE& s, const W& wgh) -> double {
    return p_curr_v[s] / static_cast<double>(G.get_vertex(s).getOutDegree());
  };
  auto reduce_f = [&] (double l, double r) { return l + r; };
  auto apply_f_v = [&] (std::tuple<uintE, double> k) {
    const uintE& w = std::get<0>(k);
    p_next_v[w] = std::get<1>(k);
    return Maybe<std::tuple<uintE, double>>();
  };

  auto map_f_u = [&] (const uintE& d, const uintE& s, const W& wgh) -> double {
    return p_curr_u[s] / static_cast<double>(G.get_vertex(s).getOutDegree());
  };
  auto apply_f_u = [&] (std::tuple<uintE, double> k) {
    const uintE& w = std::get<0>(k);
    p_next_u[w] = std::get<1>(k);
    return Maybe<std::tuple<uintE, double>>();
  };

  size_t iter = 0;
  double sim = u == v ? 1 : 0;
  while (iter++ < max_iters) {
    timer t; t.start();
    // SpMV
    timer tt; tt.start();
    auto Frontier_v_new = EM_v.template edgeMapReduce_dense<double, double>(Frontier_v, cond_f, map_f_v, reduce_f, apply_f_v, 0.0, no_output);
    auto Frontier_u_new = EM_u.template edgeMapReduce_dense<double, double>(Frontier_u, cond_f, map_f_u, reduce_f, apply_f_u, 0.0, no_output);
    tt.stop(); tt.reportTotal("em time");

    sim += ((double) pow(c, iter) * inner_product<double>(p_next_u.begin(), p_next_v.begin(), n));

    Frontier_u.del(); Frontier_v.del();
    Frontier_v = Frontier_v_new;
    Frontier_u = Frontier_u_new;

    // Check convergence: compute L1-norm between p_curr and p_next
    auto differences_v = pbbs::delayed_seq<double>(n, [&] (size_t i) {
      auto x = p_curr_v[i];
      p_curr_v[i] = 0;
      return fabs(x-p_next_v[i]);
    });
    double L1_norm_v = pbbs::reduce(differences_v, pbbs::addm<double>());

    auto differences_u = pbbs::delayed_seq<double>(n, [&] (size_t i) {
      auto x = p_curr_u[i];
      p_curr_u[i] = 0;
      return fabs(x-p_next_u[i]);
    });
    double L1_norm_u = pbbs::reduce(differences_u, pbbs::addm<double>());
    if(L1_norm_v < eps && L1_norm_u < eps) break;

    // Reset p_curr
    std::swap(p_curr_v,p_next_v);
    std::swap(p_curr_u,p_next_u);
    t.stop(); t.reportTotal("iteration time");
  }
  Frontier_u.del(); Frontier_v.del();
  auto max_pr_v = pbbslib::reduce_max(p_next_v);
  auto max_pr_u = pbbslib::reduce_max(p_next_u);

  cout << "max_pr = " << max_pr_v << ", " << max_pr_u << endl;
  cout << "sim = " << sim << endl;
}
