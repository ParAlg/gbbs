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
  inline bool cond (intT d) { return cond_true(d); }};

//vertex map function to update its p value according to PageRank equation
struct PR_Vertex_F {
  double* p_curr;
  double* p_next;
  PR_Vertex_F(double* _p_curr, double* _p_next, intE n) :
    p_curr(_p_curr), p_next(_p_next) {}
  inline bool operator () (uintE i) {
    return 1;
  }
};

//resets p
struct PR_Vertex_Reset {
  double* p_curr;
  PR_Vertex_Reset(double* _p_curr) :
    p_curr(_p_curr) {}
  inline bool operator () (uintE i) {
    p_curr[i] = 0.0;
    return 1;
  }
};


template <class Graph>
void CoSimRank_edgeMap(Graph& G, uintE v, uintE u, double eps = 0.000001, double c = 0.6, size_t max_iters = 100) {
  const uintE n = G.n;

  auto p_curr_v = pbbs::sequence<double>(n, static_cast<double>(0));
  p_curr_v[v] = static_cast<double>(1);
  auto p_next_v = pbbs::sequence<double>(n, static_cast<double>(0));
  auto frontier = pbbs::sequence<bool>(n, true);
  //frontier_v[v] = true;

  auto p_curr_u = pbbs::sequence<double>(n, static_cast<double>(0));
  p_curr_u[u] = static_cast<double>(1);
  auto p_next_u = pbbs::sequence<double>(n, static_cast<double>(0));
  //auto frontier_u = pbbs::sequence<bool>(n, true);
  //frontier_u[u] = true;

  // read from special array of just degrees

  auto degrees = pbbs::sequence<uintE>(n, [&] (size_t i) { return G.get_vertex(i).getOutDegree(); });

  vertexSubset Frontier(n,n,frontier.to_array());

  //vertexSubset Frontier_u(n,n,frontier.to_array());

  size_t iter = 0;
  double sim = u == v ? 1 : 0;
  while (iter++ < max_iters) {
    debug(timer t; t.start(););
    // SpMV
    edgeMap(G,Frontier,PR_F<Graph>(p_curr_v.begin(),p_next_v.begin(),G), 0, no_output);
    //vertexMap(Frontier_v,PR_Vertex_F(p_curr_v.begin(),p_next_v.begin(),n));

    edgeMap(G,Frontier,PR_F<Graph>(p_curr_u.begin(),p_next_u.begin(),G), 0, no_output);
    //vertexMap(Frontier_u,PR_Vertex_F(p_curr_u.begin(),p_next_u.begin(),n));

    sim += ((double) pow(c, iter) * std::inner_product(p_next_u.begin(), p_next_u.end(), p_next_v.begin(), 0));

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

  for (size_t i=0; i<100; i++) {
    std::cout << p_next_u[i] << endl;
  }
  for (size_t i=0; i<100; i++) {
    std::cout << p_next_v[i] << endl;
  }
}

template <class Graph>
void CoSimRank(Graph& G, uintE v, uintE u, double eps = 0.000001, double c = 0.6, size_t max_iters = 100) {
  using W = typename Graph::weight_type;
  const uintE n = G.n;

  auto p_curr_v = pbbs::sequence<double>(n, static_cast<double>(0));
  p_curr_v[v] = static_cast<double>(1);
  auto p_next_v = pbbs::sequence<double>(n, static_cast<double>(0));
  auto frontier_v = pbbs::sequence<bool>(n, true);
  //frontier_v[v] = true;

  auto p_curr_u = pbbs::sequence<double>(n, static_cast<double>(0));
  p_curr_u[u] = static_cast<double>(1);
  auto p_next_u = pbbs::sequence<double>(n, static_cast<double>(0));
  auto frontier_u = pbbs::sequence<bool>(n, true);
  //frontier_u[u] = true;

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

  auto cond_f = [&] (const uintE& x) { return true; };
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
  double sim = u == v ? 1 : 0;
  while (iter++ < max_iters) {
    timer t; t.start();
    // SpMV
    timer tt; tt.start();
    EM_v.template edgeMapReduce_dense<double, double>(Frontier_v, cond_f, map_f_v, reduce_f, apply_f_v, 0.0, no_output);
    EM_u.template edgeMapReduce_dense<double, double>(Frontier_u, cond_f, map_f_u, reduce_f, apply_f_u, 0.0, no_output);
    tt.stop(); tt.reportTotal("em time");

    double dot = 0;
    for (size_t o = 0; o < n; o++) {
    dot += p_next_u[o] * p_next_v[o];
    }
    sim += ((double) pow(c, iter)) * dot;

    sim += ((double) pow(c, iter) * std::inner_product(p_next_u.begin(), p_next_u.end(), p_next_v.begin(), 0));

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
  Frontier_u.del(); Frontier_v.del();
  auto max_pr_v = pbbslib::reduce_max(p_next_v);
  auto max_pr_u = pbbslib::reduce_max(p_next_u);
  cout << "max_pr = " << max_pr_v << ", " << max_pr_u << endl;
  cout << "sim = " << sim << endl;
}
