// Usage:
// numactl -i all ./wBFS -src 10012 -s -m -rounds 3 twitter_wgh_SJ
// flags:
//   required:
//     -src: the source to compute shortest path distances from
//     -w: indicate that the graph is weighted
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#define WEIGHTED 1

#include "wBFS.h"

template <class vertex>
void wBFS_runner(graph<vertex>& GA, commandLine P) {
  uintE src = P.getOptionLongValue("-src", 0);
  size_t num_buckets = P.getOptionLongValue("-nb", 32);
  bool no_blocked = P.getOptionValue("-noblocked");
  bool largemem = P.getOptionValue("-largemem");
  if (num_buckets != (1 << pbbs::log2_up(num_buckets))) {
    cout << "Please specify a number of buckets that is a power of two" << endl;
    exit(-1);
  }
  wBFS(GA, src, num_buckets, largemem, no_blocked);
}

generate_main(wBFS_runner, false);
