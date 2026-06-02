// HaMIS driver. Mirrors the HaMM driver: no algorithm-toggle flags, just
// runs the asynchronous priority-by-id maximal independent set.

#include "HaMIS.h"
#include "benchmarks/MaximalMatching/shared/decompress_to_csr.h"

#include "gbbs/gbbs.h"

#include <iostream>

namespace gbbs {

template <template <class W> class vertex, class W>
double run_mis(symmetric_graph<vertex, W>& G, commandLine P) {
  timer t;
  t.start();
  auto state = hamis::MaximalIndependentSet<vertex, W>(G);
  double tt = t.stop();
  if (P.getOption("-stats")) {
    size_t mis_size = parlay::reduce(parlay::delayed_seq<size_t>(
        G.n,
        [&](size_t i) -> size_t { return state[i] == hamis::kInMis ? 1 : 0; }));
    std::cout << "MaximalIndependentSet size: " << mis_size << "\n";
  }
  if (P.getOption("-verify")) {
    hamis::verify_mis(G, state);
  }
  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

template <template <class W> class vertex, class W>
double HaMIS_runner(symmetric_graph<vertex, W>& G, commandLine P) {
  std::cout << "### Application: HaMIS" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  bool decompress = P.getOption("-decompress");
  std::cout << "### Params: decompress=" << (decompress ? "true" : "false")
            << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  assert(P.getOption("-s"));

  if (decompress) {
    timer dec_t;
    dec_t.start();
    auto G_csr = decompress_to_csr(G);
    double dec_time = dec_t.stop();
    std::cout << "### Decompression Time: " << dec_time << std::endl;
    return run_mis(G_csr, P);
  } else {
    return run_mis(G, P);
  }
}

}  // namespace gbbs

generate_symmetric_main(gbbs::HaMIS_runner, false);
