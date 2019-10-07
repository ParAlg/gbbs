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

#include "WorkEfficientSDB14/Connectivity.h"
#include "UnionFind/Connectivity.h"

static timer bt;
using uchar = unsigned char;

#define time(_var,_body)    \
  bt.start();               \
  _body;		    \
  double _var = bt.stop();

/* ************************* Connectivity wrappers *************************** */
template <class S1, class S2>
inline void cc_check(S1& correct, S2& check);

template <class Graph>
double t_gbbs_cc(Graph& G, commandLine P, pbbs::sequence<uintE>& correct) {
  double beta = P.getOptionDoubleValue("-beta", 0.2);
  double permute = P.getOptionDoubleValue("-permute", false);
  time(t, auto CC = workefficient_cc::CC(G, beta, false, permute));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}

//template <class Graph,
//          class Find,
//          template <class F> class Unite,
//          template <class F, class U, class G> class UFTemplate>
//double t_union_find_variant(Graph& G, commandLine P, pbbs::sequence<uintE>& correct) {
//  int sampling_rounds = P.getOptionLongValue("-sample_rounds", /*default rounds=*/2);
//  auto find = Find;
//  auto unite = Unite<decltype(find)>(G.n, find);
//  auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds);
//  time(t, auto CC = q.components());
//  if (P.getOptionValue("-check")) {
//    cc_check(correct, CC);
//  }
//  return t;
//}

/* ************************* Benchmark Utils *************************** */

template<typename F>
double reduce(std::vector<double> V, F f) {
  double x = V[0];
  for (size_t i=1; i < V.size(); i++) x = f(x,V[i]);
  return x;
}

double median(std::vector<double> V) {
  std::sort(V.begin(),V.end());
  if (V.size()%2 == 1)
    return V[V.size()/2];
  else
    return (V[V.size()/2] + V[V.size()/2 - 1])/2.0;
}

double sumf(double a, double b) {return a+ b;};
double minf(double a, double b) {return (a < b) ? a : b;};
double maxf(double a, double b) {return (a > b) ? a : b;};

template<typename Graph, typename F>
std::vector<double> repeat(Graph& G, size_t rounds, pbbs::sequence<uintE>& correct, F test, commandLine P) {
  std::vector<double> R;
  for (size_t i=0; i < rounds; i++) {
    R.push_back(test(G, P, correct));
  }
  return R;
}

template<typename Graph, typename F>
bool run_multiple(Graph& G, size_t rounds, pbbs::sequence<uintE>& correct,
		  std::string name, commandLine P, F test) {
  std::vector<double> t = repeat(G, rounds, correct, test, P);

  double mint = reduce(t, minf);
  double maxt = reduce(t, maxf);
  double med = median(t);

  cout << name << std::setprecision(5)
       << ": r=" << rounds
       << ", med=" << med
       << " (" << mint << "," << maxt << "), "
       << endl;
  return 1;
}


template <class Graph>
double pick_test(Graph& G, size_t id, size_t rounds, commandLine P, pbbs::sequence<uintE>& correct) {
  switch (id) {
  case 0:
    return run_multiple(G, rounds, correct, "gbbs_cc", P, t_gbbs_cc<Graph>);
  case 1:
    return run_multiple(G, rounds, correct, "gbbs_hybridcc", P, t_gbbs_hybridcc<G>);
  case 2:
    return run_multiple(G, rounds, correct, "jayanti_rank_cc", P, t_jayanti_rank_cc<G>);
  case 3:
    return run_multiple(G, rounds, correct, "", P, t_union_find_cc<G, Unite, Find>);

  default:
    assert(false);
    std::cout << "Unknown test" << std::endl;
    exit(-1);
    return 0.0 ;
  }
}

/* ************************* Utils *************************** */

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
  std::cout << "n_cc = " << flags[n] << "\n";
  return flags[n];
}

template <class Seq>
inline size_t largest_cc(Seq& labels) {
  size_t n = labels.size();
  // could histogram to do this in parallel.
  auto flags = sequence<uintE>(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    flags[labels[i]] += 1;
  }
  size_t sz = pbbslib::reduce_max(flags);
  std::cout << "largest_cc has size: " << sz << "\n";
  return sz;
}

template <class Seq>
inline size_t RelabelDet(Seq& ids) {
  using T = typename Seq::value_type;
  size_t n = ids.size();
  auto component_map = pbbs::sequence<T>(n + 1, (T)0);
  T cur_comp = 1;
  for (size_t i=0; i<n; i++) {
    T comp = ids[i];
    T new_comp = cur_comp++;
    if (component_map[comp] == 0) {
      component_map[comp] = new_comp;
    }
    ids[i] = new_comp;
  }
  return cur_comp;
}

template <class S1, class S2>
inline void cc_check(S1& correct, S2& check) {
  RelabelDet(check);

  bool is_correct = true;
  uintE max_cor = 0;
  uintE max_chk = 0;
  for (size_t i=0; i<correct.size(); i++) {
    assert(correct[i] == check[i]);
    if ((correct[i] != check[i])) {
      is_correct = false;
      cout << "at i = " << i << " cor = " << correct[i] << " got: " << check[i] << endl;
    }
    if (correct[i] > max_cor) {
      pbbs::write_max(&max_cor, correct[i], std::less<uintE>());
    }
    if (check[i] > max_chk) {
      pbbs::write_max(&max_chk, check[i], std::less<uintE>());
    }
  }//);
  cout << "correctness check: " << is_correct << endl;
  cout << "max_cor = " << max_cor << " max_chk = " << max_chk << endl;
}

/* ************************* ***** *************************** */

template <class Graph>
double Benchmark_runner(Graph& G, commandLine P) {
  int test_num = P.getOptionIntValue("-t", -1);
  int rounds = P.getOptionIntValue("-r", 5);
  bool symmetric = P.getOptionValue("-s");
  int num_tests = 6; // update if new algorithm is added

  cout << "rounds = " << rounds << endl;
  cout << "num threads = " << num_workers() << endl;

  auto correct = pbbs::sequence<uintE>();
  if (P.getOptionValue("-check")) {
    correct = workefficient_cc::CC(G, 0.2, false, true);
    RelabelDet(correct);
  }
  if (test_num == -1) {
    for (int i=0; i < num_tests; i++) {
      pick_test(G, i, rounds, P, correct);
    }
  } else {
    pick_test(G, test_num, rounds, P, correct);
  }
  return 1.0;
}

generate_symmetric_once_main(Benchmark_runner, false);
