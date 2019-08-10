// Usage (example):
// numactl -i all ./Benchmark -alg async -rounds 3 -s -m twitter_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -alg: specify which connectivity/spanning forest algorithm to run
//
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc

#include <algorithm>
#include <vector>
#include "Benchmark.h"

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

template<typename G, typename F>
std::vector<double> repeat(G& GA, size_t rounds, pbbs::sequence<uintE>& correct, F test, commandLine P) {
  std::vector<double> R;
  for (size_t i=0; i < rounds; i++) {
    R.push_back(test(GA, P, correct));
  }
  return R;
}

template<typename G, typename F>
bool run_multiple(G& GA, size_t rounds, pbbs::sequence<uintE>& correct,
		  std::string name, commandLine P, F test) {
  std::vector<double> t = repeat(GA, rounds, correct, test, P);

  double mint = reduce(t, minf);
  double maxt = reduce(t, maxf);
  double med = median(t);
  double tt;

  cout << name << std::setprecision(5)
       << ": r=" << rounds
       << ", med=" << med
       << " (" << mint << "," << maxt << "), "
       << endl;
  return 1;
}

template <class G>
double pick_test(G& GA, size_t id, size_t rounds, commandLine P, pbbs::sequence<uintE>& correct) {
  switch (id) {
  case 0:
    return run_multiple(GA, rounds, correct, "gbbs_cc", P, t_gbbs_cc<G>);
  case 1:
    return run_multiple(GA, rounds, correct, "gbbs_sf", P, t_gbbs_sf<G>);
  case 2:
    return run_multiple(GA, rounds, correct, "async_cc", P, t_async_cc<G>);
  case 3:
    return run_multiple(GA, rounds, correct, "async_sf", P, t_async_cc<G>);
  case 4:
    return run_multiple(GA, rounds, correct, "rem_cc", P, t_rem_cc<G>);
  case 5:
    return run_multiple(GA, rounds, correct, "rem_sf", P, t_rem_sf<G>);
  case 6:
    return run_multiple(GA, rounds, correct, "ndopt_cc", P, t_ndopt_cc<G>);
  case 7:
    return run_multiple(GA, rounds, correct, "ndopt_sf", P, t_ndopt_sf<G>);
  case 8:
    return run_multiple(GA, rounds, correct, "gbbs_hybridcc", P, t_gbbs_hybridcc<G>);

  default:
    assert(false);
    exit(-1);
    return 0.0 ;
  }
}


template <class G>
double Benchmark_runner(G& GA, commandLine P) {
  int test_num = P.getOptionIntValue("-t", -1);
  int rounds = P.getOptionIntValue("-r", 5);
  int num_tests = 9; // update if new algorithm is added
  cout << "rounds = " << rounds << endl;
  cout << "num threads = " << num_workers() << endl;

  auto correct = pbbs::sequence<uintE>();
  if (P.getOptionValue("-check")); {
    correct = gbbs_cc::CC(GA, 0.2, false, true);
    RelabelDet(correct);
  }
  if (test_num == -1) {
    for (int i=0; i < num_tests; i++) {
      pick_test(GA, i, rounds, P, correct);
    }
  } else {
    pick_test(GA, test_num, rounds, P, correct);
  }
  return 1.0;
}
// "[-s (symmetrized)] [-c (compressed)] [-r <rounds>] [-t <testid>] infile"
generate_main(Benchmark_runner, false);
