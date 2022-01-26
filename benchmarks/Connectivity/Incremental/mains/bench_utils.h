#pragma once

#include <iomanip>

#include "benchmarks/Connectivity/ConnectIt/framework.h"
#include "benchmarks/Connectivity/ConnectIt/mains/check.h"
#include "gbbs/gbbs.h"

/* ************************* Benchmark Utils *************************** */

namespace gbbs {
size_t global_batch_size = 0;
double global_update_pct = 0.0;
double global_insert_to_query = 0.0;

extern bool print_batch_time;

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
  size_t cc_before, cc_after;
  for (size_t i = 0; i < rounds; i++) {
    double tot, avg, thp;
    std::tie(tot, avg, thp, cc_before, cc_after) = test(G, P);
    total.push_back(tot);
    average_batch.push_back(avg);
    thput.push_back(thp);
  }
  return std::make_tuple(total, average_batch, thput, cc_before, cc_after);
}

inline void print_cpu_stats(std::string& name, size_t rounds, size_t cc_before,
                            size_t cc_after, double medt, double mint,
                            double maxt, double med_batch, double min_batch,
                            double max_batch, double med_throughput,
                            double min_throughput, double max_throughput,
                            commandLine& P) {
  std::cout << "{" << std::endl;
  std::cout << "  \"test_type\": \"static_connectivity_result\"," << std::endl;
  std::cout << "  \"test_name\" : \"" << name << "\"," << std::endl;
  std::cout << "  \"graph\" : \"" << P.getArgument(0) << "\"," << std::endl;
  std::cout << "  \"rounds\" : " << rounds << "," << std::endl;
  std::cout << "  \"insert_to_query\" : \"" << global_insert_to_query << "\","
            << std::endl;
  std::cout << "  \"update_pct\" : \"" << global_update_pct << "\","
            << std::endl;
  std::cout << "  \"batch_size\" : \"" << global_batch_size << "\","
            << std::endl;
  std::cout << "  \"cc_before\" : \"" << cc_before << "\"," << std::endl;
  std::cout << "  \"cc_after\" : \"" << cc_after << "\"," << std::endl;
  std::cout << "  \"medt\" : " << std::setprecision(5) << medt << ","
            << std::endl;
  std::cout << "  \"mint\" : " << mint << "," << std::endl;
  std::cout << "  \"maxt\" : " << maxt << "," << std::endl;
  std::cout << "  \"med_batch_time\" : " << med_batch << "," << std::endl;
  std::cout << "  \"min_batch_time\" : " << min_batch << "," << std::endl;
  std::cout << "  \"max_batch_time\" : " << max_batch << "," << std::endl;
  std::cout << "  \"med_throughput\" : " << med_throughput << "," << std::endl;
  std::cout << "  \"min_throughput\" : " << min_throughput << "," << std::endl;
  std::cout << "  \"max_throughput\" : " << max_throughput << std::endl;
  std::cout << "}" << std::endl;
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

  print_cpu_stats(name, rounds, cc_before, cc_after, medt, mint, maxt, meda,
                  mina, maxa, medtp, mintp, maxtp, P);
  return 1;
}

template <class Graph, bool provides_initial_graph, class F>
void run_tests(Graph& G, size_t n, sequence<incremental_update>& updates,
               size_t batch_size, size_t insert_to_query, size_t rounds,
               commandLine P, F test_type, std::initializer_list<F> tests) {
  for (auto test : tests) {
    test(G, n, updates, batch_size, insert_to_query, rounds, P);
  }
}

static timer bt;
using uchar = unsigned char;

#define time(_var, _body) \
  bt.start();             \
  _body;                  \
  double _var = bt.stop();

/* ************************* Connectivity wrappers ***************************
 */
template <class S1, class S2>
inline void cc_check(S1& correct, S2& check);

namespace connectit {

template <class S>
void baseline_process_batch(sequence<uintE>& parents, S& batch) {
  auto F = connectit::get_find_function<find_compress>();
  auto U = connectit::get_unite_function<unite, decltype(F), find_compress>(
      parents.size(), F);
  for (size_t i = 0; i < batch.size(); i++) {
    auto[u, v, utype] = batch[i];
    if (utype == insertion_type) {
      U(u, v, parents);
    }
  }
}

void check_shortcut(sequence<uintE>& parents, uintE u) {
  uintE p_u = parents[u];
  while (p_u != parents[p_u]) {
    p_u = parents[p_u];
  }
  parents[u] = p_u;
}

template <class Graph, class Alg, bool provides_initial_graph,
          bool reorder_batch>
auto run_abstract_alg(Graph& G, size_t n, sequence<incremental_update>& updates,
                      size_t batch_size, size_t insert_to_query, bool check,
                      Alg& alg) {
  auto parents = sequence<uintE>::from_function(n, [&](size_t i) { return i; });

  sequence<uintE> correct_parents;

  /* compute initial components if nec. */
  if
    constexpr(provides_initial_graph) {
      auto find_fn = get_find_function<find_compress>();
      auto unite_fn =
          get_unite_function<unite, decltype(find_fn), find_compress>(n,
                                                                      find_fn);
      using UF =
          union_find::UFAlgorithm<decltype(find_fn), decltype(unite_fn), Graph>;
      auto uf_alg = UF(G, unite_fn, find_fn);

      uf_alg.template compute_components<no_sampling>(parents);
    }

  size_t ncc_before = num_cc(parents);

  if (check) {
    correct_parents = parents; /* copy */
  }

  timer init_t;
  init_t.start();
  alg.initialize(parents);
  init_t.stop();
  init_t.next("#alg initialization time");

  /* process batches */
  timer tt;
  tt.start();
  size_t m = updates.size();
  size_t n_batches = (m + batch_size - 1) / batch_size;

  std::vector<double> batch_times;
  std::cout << "## Total number of updates (all batches): " << m << std::endl;
  std::cout << "## Num batches. " << n_batches << std::endl;
  std::cout << "## Batch size. " << batch_size << std::endl;
  size_t updates_processed = 0;
  for (size_t i = 0; i < n_batches; i++) {
    size_t start = i * batch_size;
    size_t end = std::min((i + 1) * batch_size, m);
    auto update = updates.cut(start, end);
    if (i % 10000 == 0) {
      std::cout << "# starting batch" << std::endl;
    }
    timer batch_tt;
    batch_tt.start();
    alg.template process_batch<reorder_batch>(parents, update);
    double batch_time = batch_tt.stop();
    //      std::cout << "# finished batch" << std::endl;
    //      if (i % 10000 == 0) {
    //        std::cout << "## Finished : " << i << " out of " << n_batches << "
    //        batches." << std::endl;
    //      }
    batch_times.emplace_back(batch_time);
    updates_processed += update.size();
    if (print_batch_time) {
      std::cout << "batch-time: " << batch_time << std::endl;
    }

    if (check) {
      /* 1. run check algorithm */
      baseline_process_batch(correct_parents, update);

      /* 2. check that labels are consistent */
      auto parents_copy = parents;
      auto correct_parents_copy = correct_parents;
      parallel_for(0, n, [&](size_t j) {
        check_shortcut(parents_copy, j);
        check_shortcut(correct_parents_copy, j);
      });

      RelabelDet(correct_parents_copy);
      cc_check(correct_parents_copy, parents_copy);
    }
  }
  double t = tt.stop();
  double med_batch_time = median(batch_times);
  parallel_for(0, n, [&](size_t i) { check_shortcut(parents, i); });

  size_t ncc_after = num_cc(parents);

  double throughput = static_cast<double>(updates_processed) / t;
  return std::make_tuple(t, med_batch_time, throughput, ncc_before, ncc_after);
}

}  // namespace connectit
}  // namespace gbbs
