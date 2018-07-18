// Usage:
// numactl -i all ./CC -rounds 1 -pack -s -m twitter_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -pack : delete edges from G
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc

#include "CC.h"
#include "ligra.h"

template <class vertex>
void CC_runner(graph<vertex>& GA, commandLine P) {
  auto beta = P.getOptionDoubleValue("-beta", 0.2);
  auto pack = P.getOption("-pack");
  auto sym = P.getOption("-s");
  assert(sym);
  auto components = cc::CC(GA, beta, pack);
  if (P.getOption("-stats")) {
    auto cc_im = make_in_imap<uintE>(GA.n, [&] (size_t i) { return components[i];});
    cc::num_cc(cc_im);
    cc::largest_cc(cc_im);
  }
  free(components);
  if (pack) {
    // packing mutates the graph, packing out all intra-cluster edges, and can
    // only be run once unless the input graph is copied.
    exit(0);
  }
}

generate_main(CC_runner, true);
