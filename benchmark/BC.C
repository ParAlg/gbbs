// Usage:
// numactl -i all ./BC -src 10012 -s -m -rounds 3 twitter_SJ
// flags:
//   required:
//     -src: the source to compute centrality contributions from
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#include "BC.h"
#include "ligra.h"

template <class vertex>
void BC_runner(graph<vertex>& GA, commandLine P) {
  uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
  auto scores = bc::BC(GA, src);
}

generate_main(BC_runner, false);
