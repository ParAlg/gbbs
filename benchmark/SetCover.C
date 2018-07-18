// Usage:
// numactl -i all ./SetCover -nb 512 -rounds 1 -s -c clueweb_sym.bytepda
// flags:
//   optional:
//     -s : indicate that the graph is symmetric
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -nb : the number of buckets to use in the bucketing implementation

#include "SetCover.h"
#include "ligra.h"

template <class vertex>
void SetCover_runner(graph<vertex>& GA, commandLine P) {
  timer t;
  t.start();
  size_t num_buckets = P.getOptionLongValue("-nb", 128);
  cout << "here, nb = " << num_buckets << endl;
  auto cover = SetCover(GA, num_buckets);
  cover.del();
  t.stop();
  t.reportTotal();

  // Set-cover mutates the underlying graph (unless it is copied, which
  // we don't do to prevent memory issues), so we make sure the algorithm is run
  // exactly once.
  exit(0);
}

generate_main(SetCover_runner, true)
