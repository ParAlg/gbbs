// Usage:
// numactl -i all ./BFS -src 10012 -s -m -rounds 3 twitter_SJ
// flags:
//   required:
//     -src: the source to compute the BFS from
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#include "BFS.h"

template <class vertex>
void BFS_runner(graph<vertex>& GA, commandLine P) {
  uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
  auto parents = BFS(GA, src);
  parents.del();
}

generate_main(BFS_runner, false);
