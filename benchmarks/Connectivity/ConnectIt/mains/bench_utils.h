#pragma once

#include <iomanip>

#include "check.h"

/* ************************* Benchmark Utils *************************** */

namespace gbbs {
static timer bt;
using uchar = unsigned char;

#define time(_var, _body) \
  bt.start();             \
  _body;                  \
  double _var = bt.stop();

template <typename F>
double reduce(std::vector<double> V, F f) {
  double x = V[0];
  for (size_t i = 1; i < V.size(); i++) x = f(x, V[i]);
  return x;
}

double median(std::vector<double> V) {
  std::sort(V.begin(), V.end());
  if (V.size() % 2 == 1)
    return V[V.size() / 2];
  else
    return (V[V.size() / 2] + V[V.size() / 2 - 1]) / 2.0;
}

double sumf(double a, double b) { return a + b; };
double minf(double a, double b) { return (a < b) ? a : b; };
double maxf(double a, double b) { return (a > b) ? a : b; };

template <typename Graph, typename F>
std::vector<double> repeat(Graph& G, size_t rounds, sequence<parent>& correct,
                           F test, commandLine& P) {
  std::vector<double> R;
  for (size_t i = 0; i < rounds; i++) {
    auto t = test(G, P, correct);
    std::cout << "### t = " << t << std::endl;
    R.push_back(t);
  }
  return R;
}

void print_cpu_stats(std::string& name, size_t rounds, double medt, double mint,
                     double maxt, commandLine& P) {
  std::cout << "{" << std::endl;
  std::cout << "  \"test_type\": \"static_connectivity_result\"," << std::endl;
  std::cout << "  \"test_name\" : \"" << name << "\"," << std::endl;
  std::cout << "  \"graph\" : \"" << P.getArgument(0) << "\"," << std::endl;
  std::cout << "  \"rounds\" : " << rounds << "," << std::endl;
  std::cout << "  \"medt\" : " << std::setprecision(5) << medt << ","
            << std::endl;
  std::cout << "  \"mint\" : " << mint << "," << std::endl;
  std::cout << "  \"maxt\" : " << maxt << std::endl;
  std::cout << "}" << std::endl;
}

template <typename Graph, typename F>
bool run_multiple(Graph& G, size_t rounds, sequence<parent>& correct,
                  std::string name, commandLine& P, F test) {
  std::vector<double> t = repeat(G, rounds, correct, test, P);

  double mint = reduce(t, minf);
  double maxt = reduce(t, maxf);
  double med = median(t);

  print_cpu_stats(name, rounds, med, mint, maxt, P);
  return 1;
}

template <class F, class Graph>
void run_tests(Graph& G, int rounds, commandLine& P, sequence<parent>& correct,
               F test_type, std::initializer_list<F> tests) {
  for (auto test : tests) {
    test(G, rounds, P, correct);
  }
}

/* ************************* Connectivity wrappers ***************************
 */
template <class S1, class S2>
inline void cc_check(S1& correct, S2& check);

template <class Graph>
double t_gbbs_cc(Graph& G, commandLine P, sequence<parent>& correct) {
  double beta = P.getOptionDoubleValue("-beta", 0.2);
  double permute = P.getOptionDoubleValue("-permute", false);
  time(t, auto CC = workefficient_cc::CC(G, beta, false, permute));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}
}  // namespace gbbs
