// Usage:
// numactl -i all ./SCC -beta 1.5 -rounds 2 -s -m twitter_J
// flags:
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -beta <value> : the base of the exponent to use (controls how quickly
//                     vertices are added)
//     -rounds : the number of times to run the algorithm
//     -stats : print the #sccs, and the #vertices in the largest scc

#include "SCC.h"
#include "ligra.h"

template <class vertex>
void SCC_runner(graph<vertex>& GA, commandLine P) {
  double beta = P.getOptionDoubleValue("-beta", 1.1);
  timer scc_t; scc_t.start();
  auto labels = SCC(GA, beta);
  scc_t.stop(); scc_t.reportTotal("SCC time");
  if (P.getOption("-stats")) {
    num_scc(labels);
    scc_stats(labels);
  }
}

generate_main(SCC_runner, false);
