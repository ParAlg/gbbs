#include "benchmarks/Connectivity/common.h"

#include "framework.h"
#include "sampling.h"

#include "utils/benchmark.h"

void print_result(commandLine& P, std::string sampling_method, size_t rounds, double sampling_time, double pct_covered, double pct_ic_edges, size_t k=0) {
  std::cout << "{" << std::endl;
  std::cout << "  \"test_type\": \"sampling_result\"," << std::endl;
  std::cout << "  \"graph\" : \"" << P.getArgument(0) << "\"," << std::endl;
  std::cout << "  \"sampling_method\" : " << "\"" << sampling_method << "\"," << std::endl;
  std::cout << "  \"rounds\" : " << rounds << "," << std::endl;
  std::cout << "  \"sampling_time\" : " << sampling_time << "," << std::endl;
  std::cout << "  \"pct_covered\" : " << pct_covered << "," << std::endl;
  std::cout << "  \"pct_ic_edges\" : " << pct_ic_edges << "," << std::endl;
  std::cout << "  \"k\" : " << k << std::endl;
  std::cout << "}" << std::endl;
}

template <class Graph>
size_t intercomponent_edges(Graph& G, pbbs::sequence<parent>& parents) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  auto ic_edges = pbbs::sequence<size_t>(n, (size_t)0);
  parallel_for(0, n, [&] (size_t u) {
    auto p_u = parents[u];
    auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      if (p_u != parents[v]) {
        return static_cast<size_t>(1);
      }
      return static_cast<size_t>(0);
    };
    auto red_m = pbbs::addm<size_t>();
    ic_edges[u] = G.get_vertex(u).template reduceOutNgh<size_t>(u, map_f, red_m);
  });
  return pbbslib::reduce_add(ic_edges.slice());
}

template <class Graph>
void TestBFSSampling(Graph& G, commandLine& P) {
  long rounds = P.getOptionLongValue("-r", 5);
  using BFS = BFSSamplingTemplate<Graph>;
  auto sampler = BFS(G, P);
  std::vector<double> pcts;
  std::vector<double> ic_edge_pcts;

  auto test = [&] () {
    timer t; t.start();
    auto parents = sampler.initial_components();
    double sampling_time = t.stop();

    auto [frequent_comp, pct] = connectit::sample_frequent_element(parents);
    pcts.emplace_back(pct);

    size_t ic_edges = intercomponent_edges(G, parents);
    double ic_edge_pct = ((double)ic_edges) / ((double)G.m);
    ic_edge_pcts.emplace_back(ic_edge_pct);

    return sampling_time;
  };
  auto [mint, maxt, medt] = benchmark::run_multiple(rounds, test);
  auto medpct = benchmark::median(pcts);
  auto med_icpct = benchmark::median(ic_edge_pcts);
  print_result(P, "bfs", rounds, medt, medpct, med_icpct);
}

template <class Graph>
void TestLDDSampling(Graph& G, commandLine& P) {
  long rounds = P.getOptionLongValue("-r", 5);
  using LDD = LDDSamplingTemplate<Graph>;
  auto sampler = LDD(G, P);
  std::vector<double> pcts;
  std::vector<double> ic_edge_pcts;
  auto test = [&] () {
    timer t; t.start();
    auto parents = sampler.initial_components();
    double sampling_time = t.stop();

    auto [frequent_comp, pct] = connectit::sample_frequent_element(parents);
    pcts.emplace_back(pct);

    size_t ic_edges = intercomponent_edges(G, parents);
    double ic_edge_pct = ((double)ic_edges) / ((double)G.m);
    ic_edge_pcts.emplace_back(ic_edge_pct);

    return sampling_time;
  };

  auto [mint, maxt, medt] = benchmark::run_multiple(rounds, test);
  auto medpct = benchmark::median(pcts);
  auto med_icpct = benchmark::median(ic_edge_pcts);
  print_result(P, "ldd", rounds, medt, medpct, med_icpct);
}

template <class Graph>
double TestKOutSampling(Graph& G, commandLine& P, int neighbor_rounds) {
  long rounds = P.getOptionLongValue("-r", 5);
  size_t n = G.n;
  std::vector<double> pcts;
  std::vector<double> ic_edge_pcts;
  auto test = [&] () {
    using KOut = KOutSamplingTemplate<Graph>;
    auto sampler = KOut(G, P, neighbor_rounds);

    timer t; t.start();
    auto parents = sampler.initial_components();
    double sampling_time = t.stop();

    auto [frequent_comp, pct] = connectit::sample_frequent_element(parents);
    pcts.emplace_back(pct);

    size_t ic_edges = intercomponent_edges(G, parents);
    double ic_edge_pct = ((double)ic_edges) / ((double)G.m);
    ic_edge_pcts.emplace_back(ic_edge_pct);

    return sampling_time;
  };

  auto [mint, maxt, medt] = benchmark::run_multiple(rounds, test);
  auto medpct = benchmark::median(pcts);
  auto med_icpct = benchmark::median(ic_edge_pcts);
  print_result(P, "kout", rounds, medt, medpct, med_icpct, neighbor_rounds);
}

template <class Graph>
double TestKOutAfforestSampling(Graph& G, commandLine& P, int neighbor_rounds) {
  long rounds = P.getOptionLongValue("-r", 5);
  size_t n = G.n;
  std::vector<double> pcts;
  std::vector<double> ic_edge_pcts;
  auto test = [&] () {
    using KOut = KOutSamplingTemplate<Graph>;
    auto sampler = KOut(G, P, neighbor_rounds);

    timer t; t.start();
    auto parents = sampler.initial_components_afforest();
    double sampling_time = t.stop();

    auto [frequent_comp, pct] = connectit::sample_frequent_element(parents);
    pcts.emplace_back(pct);

    size_t ic_edges = intercomponent_edges(G, parents);
    double ic_edge_pct = ((double)ic_edges) / ((double)G.m);
    ic_edge_pcts.emplace_back(ic_edge_pct);

    return sampling_time;
  };

  auto [mint, maxt, medt] = benchmark::run_multiple(rounds, test);
  auto medpct = benchmark::median(pcts);
  auto med_icpct = benchmark::median(ic_edge_pcts);
  print_result(P, "kout-afforest", rounds, medt, medpct, med_icpct, neighbor_rounds);
}

template <class Graph>
double TestKOutPureSampling(Graph& G, commandLine& P, int neighbor_rounds) {
  long rounds = P.getOptionLongValue("-r", 5);
  size_t n = G.n;
  std::vector<double> pcts;
  std::vector<double> ic_edge_pcts;
  auto test = [&] () {
    using KOut = KOutSamplingTemplate<Graph>;
    auto sampler = KOut(G, P, neighbor_rounds);

    timer t; t.start();
    auto parents = sampler.initial_components_pure();
    double sampling_time = t.stop();

    auto [frequent_comp, pct] = connectit::sample_frequent_element(parents);
    pcts.emplace_back(pct);

    size_t ic_edges = intercomponent_edges(G, parents);
    double ic_edge_pct = ((double)ic_edges) / ((double)G.m);
    ic_edge_pcts.emplace_back(ic_edge_pct);

    return sampling_time;
  };

  auto [mint, maxt, medt] = benchmark::run_multiple(rounds, test);
  auto medpct = benchmark::median(pcts);
  auto med_icpct = benchmark::median(ic_edge_pcts);
  print_result(P, "kout-pure", rounds, medt, medpct, med_icpct, neighbor_rounds);
}

template <class Graph>
double TestKOutMaxDegreeSampling(Graph& G, commandLine& P, int neighbor_rounds) {
  long rounds = P.getOptionLongValue("-r", 5);
  size_t n = G.n;
  std::vector<double> pcts;
  std::vector<double> ic_edge_pcts;
  auto test = [&] () {
    using KOut = KOutSamplingTemplate<Graph>;
    auto sampler = KOut(G, P, neighbor_rounds);

    timer t; t.start();
    auto parents = sampler.initial_components_max_degree();
    double sampling_time = t.stop();

    auto [frequent_comp, pct] = connectit::sample_frequent_element(parents);
    pcts.emplace_back(pct);

    size_t ic_edges = intercomponent_edges(G, parents);
    double ic_edge_pct = ((double)ic_edges) / ((double)G.m);
    ic_edge_pcts.emplace_back(ic_edge_pct);

    return sampling_time;
  };

  auto [mint, maxt, medt] = benchmark::run_multiple(rounds, test);
  auto medpct = benchmark::median(pcts);
  auto med_icpct = benchmark::median(ic_edge_pcts);
  print_result(P, "kout-maxdeg", rounds, medt, medpct, med_icpct, neighbor_rounds);
}

template <class Graph>
double Sampler(Graph& G, commandLine& P) {
  std::cout << "[" << std::endl;
  TestBFSSampling<Graph>(G, P);
  std::cout << "," << endl;
  TestLDDSampling<Graph>(G, P);
  std::cout << "," << endl;

  TestKOutSampling<Graph>(G, P, 1);
  std::cout << "," << endl;
  TestKOutSampling<Graph>(G, P, 2);
  std::cout << "," << endl;
  TestKOutSampling<Graph>(G, P, 3);
  std::cout << "," << endl;
  TestKOutSampling<Graph>(G, P, 4);
  std::cout << "," << endl;

  TestKOutAfforestSampling<Graph>(G, P, 1);
  std::cout << "," << endl;
  TestKOutAfforestSampling<Graph>(G, P, 2);
  std::cout << "," << endl;
  TestKOutAfforestSampling<Graph>(G, P, 3);
  std::cout << "," << endl;
  TestKOutAfforestSampling<Graph>(G, P, 4);
  std::cout << "," << endl;

  TestKOutPureSampling<Graph>(G, P, 1);
  std::cout << "," << std::endl;
  TestKOutPureSampling<Graph>(G, P, 2);
  std::cout << "," << std::endl;
  TestKOutPureSampling<Graph>(G, P, 3);
  std::cout << "," << std::endl;
  TestKOutPureSampling<Graph>(G, P, 4);

  TestKOutMaxDegreeSampling<Graph>(G, P, 1);
  std::cout << "," << std::endl;
  TestKOutMaxDegreeSampling<Graph>(G, P, 2);
  std::cout << "," << std::endl;
  TestKOutMaxDegreeSampling<Graph>(G, P, 3);
  std::cout << "," << std::endl;
  TestKOutMaxDegreeSampling<Graph>(G, P, 4);

  std::cout << "]" << std::endl;
}

generate_symmetric_once_main(Sampler, false);
