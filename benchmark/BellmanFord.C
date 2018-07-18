// Usage:
// numactl -i all ./BellmanFord -src 10012 -s -m -rounds 3 twitter_wgh_SJ
// flags:
//   required:
//     -src: the source to compute shortest path distances from
//     -w: indicate that the graph is weighted
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#include "BellmanFord.h"

template <class vertex>
void BellmanFord_runner(graph<vertex>& GA, commandLine P) {
  uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
  auto distances = BellmanFord(GA, src);
}

generate_main(BellmanFord_runner, false);
