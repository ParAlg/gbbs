#pragma once

#include "check.h"

/* ************************* Benchmark Utils *************************** */

static timer bt;
using uchar = unsigned char;

#define time(_var,_body)    \
  bt.start();               \
  _body;		    \
  double _var = bt.stop();

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
std::vector<double> repeat(Graph& G, size_t rounds, pbbs::sequence<parent>& correct, F test, commandLine& P) {
  std::vector<double> R;
  for (size_t i=0; i < rounds; i++) {
#ifdef REPORT_PATH_LENGTHS
    max_pathlen.reset();
    total_pathlen.reset();
#endif
#ifdef REPORT_MAX_TRIES
    max_uf_tries.reset();
    total_uf_tries.reset();
#endif
    auto t = test(G, P, correct);
    std::cout << "### t = " << t << std::endl;
    R.push_back(t);
  }
  return R;
}

template<typename Graph, typename F>
bool run_multiple(Graph& G, size_t rounds, pbbs::sequence<parent>& correct,
		  std::string name, commandLine& P, F test) {
  std::vector<double> t = repeat(G, rounds, correct, test, P);

  double mint = reduce(t, minf);
  double maxt = reduce(t, maxf);
  double med = median(t);

  cout << "Test = {"
       << "\"name\": \"" << name << "\"" << std::setprecision(5)
       << ", \"rounds\":" << rounds
       << ", \"med\":" << med
       << ", \"min\":" << mint
       << ", \"max\":" << maxt << "}"
       << endl;
  return 1;
}

/* ************************* ***** *************************** */

template <class CPUStats>
void print_cpu_stats(CPUStats& stats, commandLine& P) {
  std::cout <<
    "Stats = { \"ipc\":" + std::to_string(stats.get_ipc())
    + " ,\"total_cycles\":" + std::to_string(stats.get_total_cycles())
    + " ,\"l2_hit_ratio\":" + std::to_string(stats.get_l2_hit_ratio())
    + " ,\"l3_hit_ratio\":" + std::to_string(stats.get_l3_hit_ratio())
    + " ,\"l2_misses\":" + std::to_string(stats.get_l2_misses())
    + " ,\"l2_hits\":" + std::to_string(stats.get_l2_hits())
    + " ,\"l3_misses\":" + std::to_string(stats.get_l3_misses())
    + " ,\"l3_hits\":" + std::to_string(stats.get_l3_hits())
    + " ,\"throughput\":" + std::to_string(stats.get_throughput())
    + " ,\"max_path_len\":" + std::to_string(max_pathlen.get_value())
    + " ,\"total_path_len\":" + std::to_string(total_pathlen.get_value())
    + " ,\"max_uf_tries\":" + std::to_string(max_uf_tries.get_value())
    + " ,\"total_uf_tries\":" + std::to_string(total_uf_tries.get_value())
    + "}" << std::endl;
}


/* ************************* Connectivity wrappers *************************** */
template <class S1, class S2>
inline void cc_check(S1& correct, S2& check);


template <class Graph>
double t_gbbs_cc(Graph& G, commandLine P, pbbs::sequence<parent>& correct) {
  double beta = P.getOptionDoubleValue("-beta", 0.2);
  double permute = P.getOptionDoubleValue("-permute", false);
  time(t, auto CC = workefficient_cc::CC(G, beta, false, permute));
  if (P.getOptionValue("-check")) {
    cc_check(correct, CC);
  }
  return t;
}
