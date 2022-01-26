#include "benchmarks/Connectivity/common.h"

#include "framework.h"
#include "sampling.h"

#include "utils/benchmark.h"

namespace gbbs {

void print_result(commandLine& P, std::string method, size_t rounds,
                  double time) {
  std::cout << "{" << std::endl;
  std::cout << "  \"test_type\": \"edge_gather_result\"," << std::endl;
  std::cout << "  \"graph\" : \"" << P.getArgument(0) << "\"," << std::endl;
  std::cout << "  \"method\" : "
            << "\"" << method << "\"," << std::endl;
  std::cout << "  \"rounds\" : " << rounds << "," << std::endl;
  std::cout << "  \"time\" : " << time << "," << std::endl;
  std::cout << "}" << std::endl;
}

template <class Graph>
double TestEdgeGather(Graph& G, commandLine& P) {
  using W = typename Graph::weight_type;
  long rounds = P.getOptionLongValue("-r", 5);
  size_t n = G.n;
  auto test = [&]() {
    auto arr = sequence<size_t>(n, (size_t)1);
    auto out = sequence<size_t>(n);

    timer t;
    t.start();
    parallel_for(0, n,
                 [&](size_t u) {
                   auto u_vtx = G.get_vertex(u);
                   auto map_f = [&](const uintE& u, const uintE& v,
                                    const W& wgh) { return arr[v]; };
                   auto reduce_m = parlay::addm<size_t>();
                   u_vtx.out_neighbors().reduce(map_f, reduce_m);
                 },
                 512);
    double edge_gather_time = t.stop();

    return edge_gather_time;
  };

  auto[mint, maxt, medt] = benchmark::run_multiple(rounds, test);
  print_result(P, "edge-gather", rounds, medt);
}

template <class Graph>
double TestEdgeMap(Graph& G, commandLine& P) {
  using W = typename Graph::weight_type;
  long rounds = P.getOptionLongValue("-r", 5);
  size_t n = G.n;
  auto test = [&]() {
    auto out = sequence<size_t>(n);

    timer t;
    t.start();
    parallel_for(0, n,
                 [&](size_t u) {
                   auto u_vtx = G.get_vertex(u);
                   auto map_f = [&](const uintE& u, const uintE& v,
                                    const W& wgh) {
                     return static_cast<size_t>(1);
                   };
                   auto reduce_m = parlay::addm<size_t>();
                   out[u] = u_vtx.out_neighbors().reduce(map_f, reduce_m);
                 },
                 512);
    double edge_gather_time = t.stop();

    return edge_gather_time;
  };

  auto[mint, maxt, medt] = benchmark::run_multiple(rounds, test);
  print_result(P, "edge-map", rounds, medt);
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
