#pragma once

/* ************************* Benchmark Utils *************************** */

namespace gbbs {
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
auto repeat(Graph& G, size_t rounds, F test, commandLine& P) {
  std::vector<double> total;
  std::vector<double> average_batch;
  std::vector<double> thput;
  for (size_t i = 0; i < rounds; i++) {
    double tot, avg, thp;
    std::tie(tot, avg, thp) = test(G, P);
    total.push_back(tot);
    average_batch.push_back(avg);
    thput.push_back(thp);
  }
  return std::make_tuple(total, average_batch, thput);
}

template <typename Graph, typename F>
bool run_multiple(Graph& G, size_t rounds, std::string name, commandLine& P,
                  F test) {
  std::vector<double> t; /* total */
  std::vector<double> a;
  std::vector<double> tp;
  size_t cc_before, cc_after;
  std::tie(t, a, tp, cc_before, cc_after) = repeat(G, rounds, test, P);

  double mint = reduce(t, minf);
  double maxt = reduce(t, maxf);
  double medt = median(t);

  double mina = reduce(a, minf);
  double maxa = reduce(a, maxf);
  double meda = median(a);

  double mintp = reduce(tp, minf);
  double maxtp = reduce(tp, maxf);
  double medtp = median(tp);

  std::cout << "Test = {"
            << "\"name\": \"" << name << "\"" << std::setprecision(5)
            << ", \"rounds\":" << rounds << ", \"med_time\":" << medt
            << ", \"min_time\":" << mint << ", \"max_time\":" << maxt
            << ", \"med_batch_time\":" << meda
            << ", \"min_batch_time\":" << mina
            << ", \"max_batch_time\":" << maxa
            << ", \"med_throughput\":" << medtp
            << ", \"min_throughput\":" << mintp
            << ", \"max_throughput\":" << maxtp << "}" << std::endl;
  return 1;
}

template <class Stats>
void print_cpu_stats(Stats& stats, commandLine& P) {
  std::cout
      << "Stats = { \"ipc\":" + std::to_string(stats.get_ipc()) +
             " ,\"total_cycles\":" + std::to_string(stats.total_time_cycles()) +
             " ,\"l2_hit_ratio\":" + std::to_string(stats.get_l2_hit_ratio()) +
             " ,\"l3_hit_ratio\":" + std::to_string(stats.get_l3_hit_ratio()) +
             " ,\"l2_misses\":" + std::to_string(stats.get_l2_misses()) +
             " ,\"l2_hits\":" + std::to_string(stats.get_l2_hits()) +
             " ,\"l3_misses\":" + std::to_string(stats.get_l3_misses()) +
             " ,\"l3_hits\":" + std::to_string(stats.get_l3_hits()) +
             " ,\"throughput\":" + std::to_string(stats.get_throughput()) + "}"
      << std::endl;
}

}  // namespace gbbs
