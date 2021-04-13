#include "benchmarks/Connectivity/common.h"

#include "framework.h"
#include "sampling.h"

#include "utils/benchmark.h"

namespace gbbs {
template <class CpuStats>
void print_result(commandLine& P, std::string method, size_t rounds, double time, CpuStats& stats) {
  std::cout << "{" << std::endl;
  std::cout << "  \"test_type\": \"edge_gather_result\"," << std::endl;
  std::cout << "  \"graph\" : \"" << P.getArgument(0) << "\"," << std::endl;
  std::cout << "  \"method\" : " << "\"" << method << "\"," << std::endl;
  std::cout << "  \"rounds\" : " << rounds << "," << std::endl;
  std::cout << "  \"time\" : " << time << "," << std::endl;
  std::cout << "  \"ipc\" : " << std::to_string(stats.get_ipc()) << "," << std::endl;
  std::cout << "  \"total_cycles\" : " << std::to_string(stats.get_total_cycles()) << "," << std::endl;
  std::cout << "  \"l2_hit_ratio\" : " << std::to_string(stats.get_l2_hit_ratio()) << "," << std::endl;
  std::cout << "  \"l3_hit_ratio\" : " << std::to_string(stats.get_l3_hit_ratio()) << "," << std::endl;
  std::cout << "  \"l2_misses\" : " << std::to_string(stats.get_l2_misses()) << "," << std::endl;
  std::cout << "  \"l2_hits\" : " << std::to_string(stats.get_l2_hits()) << "," << std::endl;
  std::cout << "  \"l3_misses\" : " << std::to_string(stats.get_l3_misses()) << "," << std::endl;
  std::cout << "  \"l3_hits\" : " << std::to_string(stats.get_l3_hits()) << "," << std::endl;
  std::cout << "  \"throughput\" : " << std::to_string(stats.get_throughput()) << "," << std::endl;
  std::cout << "  \"max_path_len\" : " << std::to_string(max_pathlen.get_value()) << "," << std::endl;
  std::cout << "  \"total_path_len\" : " << std::to_string(total_pathlen.get_value()) << std::endl;
  std::cout << "}" << std::endl;
}

template <class Graph>
double TestEdgeGather(Graph& G, commandLine& P) {
  using W = typename Graph::weight_type;
  long rounds = P.getOptionLongValue("-r", 5);
  size_t n = G.n;
  auto test = [&] () {
    auto arr = sequence<size_t>(n, (size_t)1);
    auto out = sequence<size_t>(n);

    timer t;
    t.start();
    parallel_for(0, n, [&] (size_t u) {
      auto u_vtx = G.get_vertex(u);
      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        return arr[v];
      };
      auto reduce_m = pbbslib::addm<size_t>();
      u_vtx.out_neighbors().reduce(map_f, reduce_m);
    }, 512);
    double edge_gather_time = t.stop();

    return edge_gather_time;
  };

#ifdef USE_PCM_LIB
  auto before_state = get_pcm_state();
  timer ot; ot.start();
#endif
  auto [mint, maxt, medt] = benchmark::run_multiple(rounds, test);
#ifdef USE_PCM_LIB
  double elapsed = ot.stop();
  auto after_state = get_pcm_state();
  cpu_stats stats = get_pcm_stats(before_state, after_state, elapsed, rounds);
#else
  cpu_stats stats;
#endif
  print_result(P, "edge-gather", rounds, medt, stats);
}

template <class Graph>
double TestEdgeMap(Graph& G, commandLine& P) {
  using W = typename Graph::weight_type;
  long rounds = P.getOptionLongValue("-r", 5);
  size_t n = G.n;
  auto test = [&] () {
    auto out = sequence<size_t>(n);

    timer t;
    t.start();
    parallel_for(0, n, [&] (size_t u) {
      auto u_vtx = G.get_vertex(u);
      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        return static_cast<size_t>(1);
      };
      auto reduce_m = pbbslib::addm<size_t>();
      out[u] = u_vtx.out_neighbors().reduce(map_f, reduce_m);
    }, 512);
    double edge_gather_time = t.stop();

    return edge_gather_time;
  };

#ifdef USE_PCM_LIB
  auto before_state = get_pcm_state();
  timer ot; ot.start();
#endif
  auto [mint, maxt, medt] = benchmark::run_multiple(rounds, test);
#ifdef USE_PCM_LIB
  double elapsed = ot.stop();
  auto after_state = get_pcm_state();
  cpu_stats stats = get_pcm_stats(before_state, after_state, elapsed, rounds);
#else
  cpu_stats stats;
#endif

  print_result(P, "edge-map", rounds, medt, stats);
}

template <class Graph>
double Sampler(Graph& G, commandLine& P) {
  std::cout << "[" << std::endl;
  TestEdgeGather<Graph>(G, P);
  std::cout << "," << std::endl;
  TestEdgeMap<Graph>(G, P);

  std::cout << "]" << std::endl;
}
}  // namespace gbbs

generate_symmetric_once_main(gbbs::Sampler, false);
