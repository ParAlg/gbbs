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

#include "edge_map_reduce.h"
#include "ligra.h"
#include "math.h"

template <class G>
struct PR_F {
  using W = typename G::weight_type;
  G& GA;
  double* p_curr, *p_next;
  PR_F(double* _p_curr, double* _p_next, G& GA) :
    p_curr(_p_curr), p_next(_p_next), GA(GA) {}
  inline bool update(const uintE& s, const uintE& d, const W& wgh){ //update function applies PageRank equation
    p_next[d] += p_curr[s]/GA.get_vertex(s).getOutDegree();
    return 1;
  }
  inline bool updateAtomic (const uintE& s, const uintE& d, const W& wgh) { //atomic Update
    pbbs::fetch_and_add(&p_next[d],p_curr[s]/GA.get_vertex(s).getOutDegree());
    return 1;
  }
  inline bool cond (intT d) { return cond_true(d); }};

//vertex map function to update its p value according to PageRank equation
struct PR_Vertex_F {
  double* p_curr;
  double* p_next;
  double damping;
  double addedConstant;
  PR_Vertex_F(double* _p_curr, double* _p_next, double _damping, intE n) :
    p_curr(_p_curr), p_next(_p_next),
    damping(_damping), addedConstant((1-_damping)*(1/(double)n)){}
  inline bool operator () (uintE i) {
    p_next[i] = damping*p_next[i] + addedConstant;
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


template <class G>
void PageRank_edgeMap(G& GA, double eps = 0.000001, size_t max_iters = 100) {
  const uintE n = GA.n;
  const double damping = 0.85;

  double one_over_n = 1/(double)n;
  auto p_curr = pbbs::sequence<double>(n, one_over_n);
  auto p_next = pbbs::sequence<double>(n, static_cast<double>(0));
  auto frontier = pbbs::sequence<bool>(n, true);

  // read from special array of just degrees

  auto degrees = pbbs::sequence<uintE>(n, [&] (size_t i) { return GA.get_vertex(i).getOutDegree(); });

  vertexSubset Frontier(n,n,frontier.to_array());

  size_t iter = 0;
  while (iter++ < max_iters) {
    debug(timer t; t.start(););
    // SpMV
    edgeMap(GA,Frontier,PR_F<G>(p_curr.begin(),p_next.begin(),GA), 0, no_output);
    vertexMap(Frontier,PR_Vertex_F(p_curr.begin(),p_next.begin(),damping,n));

    // Check convergence: compute L1-norm between p_curr and p_next
    auto differences = pbbs::delayed_seq<double>(n, [&] (size_t i) {
      return fabs(p_curr[i]-p_next[i]);
    });
    double L1_norm = pbbs::reduce(differences, pbbs::addm<double>());
    if(L1_norm < eps) break;

    debug(cout << "L1_norm = " << L1_norm << endl;);
    // Reset p_curr
    parallel_for(0, n, [&] (size_t i) { p_curr[i] = static_cast<double>(0); });
    std::swap(p_curr,p_next);

    debug(t.stop(); t.reportTotal("iteration time"););
  }
  Frontier.del();
  auto max_pr = pbbslib::reduce_max(p_next);
  cout << "max_pr = " << max_pr << endl;
  for (size_t i=0; i<100; i++) {
    cout << p_next[i] << endl;
  }
}

template <class G>
void PageRank(G& GA, double eps = 0.000001, size_t max_iters = 100) {
  using W = typename G::weight_type;
  const uintE n = GA.n;
  const double damping = 0.85;
  const double addedConstant = (1 -damping)*(1/static_cast<double>(n));

  double one_over_n = 1/(double)n;
  auto p_curr = pbbs::sequence<double>(n, one_over_n);
  auto p_next = pbbs::sequence<double>(n, static_cast<double>(0));
  auto frontier = pbbs::sequence<bool>(n, true);
  auto p_div = pbbs::sequence<double>(n, [&] (size_t i) -> double {
    return one_over_n / static_cast<double>(GA.get_vertex(i).getOutDegree());
  });

  // read from special array of just degrees

  auto degrees = pbbs::sequence<uintE>(n, [&] (size_t i) { return GA.get_vertex(i).getOutDegree(); });

  vertexSubset Frontier(n,n,frontier.to_array());
  auto EM = EdgeMap<double, G>(GA, std::make_tuple(UINT_E_MAX, static_cast<double>(0)), (size_t)GA.m/1000);

  auto cond_f = [&] (const uintE& v) { return true; };
  auto map_f = [&] (const uintE& d, const uintE& s, const W& wgh) -> double {
    return p_div[s];
//    return p_curr[s] / degrees[s];
  };
  auto reduce_f = [&] (double l, double r) { return l + r; };
  auto apply_f = [&] (std::tuple<uintE, double> k) {
    const uintE& u = std::get<0>(k);
    const double& contribution = std::get<1>(k);
    p_next[u] = damping*contribution + addedConstant;
    p_div[u] = p_next[u]/static_cast<double>(degrees[u]);
    return Maybe<std::tuple<uintE, double>>();
  };

  size_t iter = 0;
  while (iter++ < max_iters) {
    timer t; t.start();
    // SpMV
    timer tt; tt.start();
    EM.template edgeMapReduce_dense<double, double>(Frontier, cond_f, map_f, reduce_f, apply_f, 0.0, no_output);
    tt.stop(); tt.reportTotal("em time");

    // Check convergence: compute L1-norm between p_curr and p_next
    auto differences = pbbs::delayed_seq<double>(n, [&] (size_t i) {
      auto d = p_curr[i];
      p_curr[i] = 0;
      return fabs(d-p_next[i]);
    });
    double L1_norm = pbbs::reduce(differences, pbbs::addm<double>());
    if(L1_norm < eps) break;
    debug(cout << "L1_norm = " << L1_norm << endl;);

    // Reset p_curr
    std::swap(p_curr,p_next);
    t.stop(); t.reportTotal("iteration time");
  }
  Frontier.del();
  auto max_pr = pbbslib::reduce_max(p_next);
  cout << "max_pr = " << max_pr << endl;
  for (size_t i=0; i<100; i++) {
    cout << p_next[i] << endl;
  }
}

namespace delta {

struct delta_and_degree {
  double delta;
  double delta_over_degree;
};

template <class G>
struct PR_Delta_F {
  using W = typename G::weight_type;
  delta_and_degree* Delta;
  double* nghSum;
  PR_Delta_F(delta_and_degree* _Delta, double* _nghSum) :
    Delta(_Delta), nghSum(_nghSum) {}
  inline bool update(const uintE& s, const uintE& d, const W& wgh){
    double oldVal = nghSum[d];
    nghSum[d] += Delta[s].delta_over_degree; // Delta[s].delta/Delta[s].degree; // V[s].getOutDegree();
    return oldVal == 0;
  }
  inline bool updateAtomic (const uintE& s, const uintE&  d, const W& wgh) {
    volatile double oldV, newV;
    do { //basically a fetch-and-add
      oldV = nghSum[d]; newV = oldV + Delta[s].delta_over_degree; // Delta[s]/V[s].getOutDegree();
    } while(!pbbs::atomic_compare_and_swap(&nghSum[d],oldV,newV));
    return oldV == 0.0;
  }
  inline bool cond (uintE d) { return cond_true(d); }
};

template <class G, class E>
void sparse_or_dense(G& GA, E& EM, vertexSubset& Frontier, delta_and_degree* Delta, double* nghSum, const flags fl) {
  using W = typename G::weight_type;

  if (Frontier.size() > GA.n/5) {
    Frontier.toDense();

    auto cond_f = [&] (size_t i) {
      return true;
    };
    auto map_f = [&] (const uintE& s, const uintE& d, const W& wgh) -> double {
      if (Frontier.d[d]) {
        return Delta[d].delta_over_degree; // Delta[d]/GA.V[d].getOutDegree();
      } else {
        return static_cast<double>(0);
      }
    };
    auto reduce_f = [&] (double l, double r) { return l + r; };
    auto apply_f = [&] (std::tuple<uintE, double> k) {
      const uintE& u = std::get<0>(k);
      const double& contribution = std::get<1>(k);
      nghSum[u] = contribution;
      return Maybe<std::tuple<uintE, pbbs::empty>>();
    };
    double id = 0.0;

    flags dense_fl = fl;
    dense_fl ^= in_edges; // todo: check
    timer dt; dt.start();
    EM.template edgeMapReduce_dense<pbbs::empty, double>(Frontier, cond_f, map_f, reduce_f, apply_f, id, dense_fl | no_output);

    dt.stop(); dt.reportTotal("dense time");
  } else {
    edgeMap(GA,Frontier,PR_Delta_F<G>(Delta,nghSum),GA.m/2, no_output);
  }
}

template <class G>
struct PR_Vertex_F_FirstRound {
  double damping, addedConstant, one_over_n, epsilon2;
  double* p;
  delta_and_degree* Delta;
  double* nghSum;
  G& get_degree;
  PR_Vertex_F_FirstRound(double* _p, delta_and_degree* _Delta, double* _nghSum, double _damping, double _one_over_n,double _epsilon2, G& get_degree) :
    damping(_damping),
    addedConstant((1-_damping)*_one_over_n),
    one_over_n(_one_over_n),
    epsilon2(_epsilon2),
    p(_p),
    Delta(_Delta),
    nghSum(_nghSum),
    get_degree(get_degree) {}
  inline bool operator () (uintE i) {
    double pre_init = damping*nghSum[i]+addedConstant;
    p[i] += pre_init;
    double new_delta = pre_init - one_over_n; //subtract off delta from initialization
    Delta[i].delta = new_delta;
    Delta[i].delta_over_degree = new_delta/get_degree(i);
    return (new_delta > epsilon2 * p[i]);
  }
};


template <class G>
auto make_PR_Vertex_F_FirstRound(double* p, delta_and_degree* delta, double* nghSum, double damping, double one_over_n, double epsilon2, G& get_degree) {
  return PR_Vertex_F_FirstRound<G>(p, delta, nghSum, damping, one_over_n, epsilon2, get_degree);
}

template <class G>
struct PR_Vertex_F {
  double* p;
  delta_and_degree* Delta;
  double* nghSum;
  double damping, epsilon2;
  G& get_degree;
  PR_Vertex_F(double* _p, delta_and_degree* _Delta, double* _nghSum, double _damping, double _epsilon2, G& get_degree) :
    p(_p),
    Delta(_Delta), nghSum(_nghSum),
    damping(_damping), epsilon2(_epsilon2), get_degree(get_degree) {}
  inline bool operator () (uintE i) {
    double new_delta = nghSum[i]*damping;
    Delta[i].delta = new_delta;
    Delta[i].delta_over_degree = new_delta/get_degree(i);

    if (fabs(Delta[i].delta) > epsilon2*p[i]) {
      p[i]+=new_delta;
      return 1;
    }
    else return 0;
  }
};

template <class G>
auto make_PR_Vertex_F(double* p, delta_and_degree* delta, double* nghSum, double damping, double epsilon2, G& get_degree) {
  return PR_Vertex_F<G>(p, delta, nghSum, damping, epsilon2, get_degree);
}

template <class G>
void PageRankDelta(G& GA, double eps=0.000001, double local_eps=0.01, size_t max_iters=100) {
  using W = typename G::weight_type;
  const long n = GA.n;
  const double damping = 0.85;

  double one_over_n = 1/(double)n;
  auto p = pbbs::sequence<double>(n);
  auto Delta = pbbs::sequence<delta_and_degree>(n);
  auto nghSum = pbbs::sequence<double>(n);
  auto frontier = pbbs::sequence<bool>(n);
  parallel_for(0, n, [&] (size_t i) {
    uintE degree = GA.get_vertex(i).getOutDegree();
    p[i] = 0.0;//one_over_n;
    Delta[i].delta = one_over_n; //initial delta propagation from each vertex
    Delta[i].delta_over_degree = one_over_n/degree;
    nghSum[i] = 0.0;
    frontier[i] = 1;
  });

  auto get_degree = [&] (size_t i) { return GA.get_vertex(i).getOutDegree(); };
  auto EM = EdgeMap<double, G>(GA, std::make_tuple(UINT_E_MAX, (double)0.0), (size_t)GA.m/1000);
  vertexSubset Frontier(n,n,frontier.to_array());
  auto all = pbbs::sequence<bool>(n, true);
  vertexSubset All(n,n,all.to_array()); //all vertices

  size_t round = 0;
  while(round++ < max_iters) {
    timer t; t.start();
    sparse_or_dense(GA, EM, Frontier, Delta.begin(), nghSum.begin(), no_output);
    vertexSubset active
      = (round == 1) ?
      vertexFilter2(All,delta::make_PR_Vertex_F_FirstRound(p.begin(),Delta.begin(),nghSum.begin(),damping,one_over_n,local_eps,get_degree)) :
      vertexFilter2(All,delta::make_PR_Vertex_F(p.begin(),Delta.begin(),nghSum.begin(),damping,local_eps,get_degree));

    // Check convergence: compute L1-norm between p_curr and p_next
    auto differences = pbbs::delayed_seq<double>(n, [&] (size_t i) {
      return fabs(Delta[i].delta);
    });
    double L1_norm = pbbs::reduce(differences, pbbs::addm<double>());
    if(L1_norm < eps) break;
    debug(cout << "L1_norm = " << L1_norm << endl;);

    // Reset
    parallel_for(0, n, [&] (size_t i) { nghSum[i] = static_cast<double>(0); });

    Frontier.del();
    Frontier = active;
    debug(t.stop(); t.reportTotal("iteration time"););
  }
  auto max_pr = pbbslib::reduce_max(p);
  cout << "max_pr = " << max_pr << endl;

  cout << "Num rounds = " << round << endl;
  Frontier.del(); All.del();
}

}
