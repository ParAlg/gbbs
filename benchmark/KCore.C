// Usage:
// numactl -i all ./KCore -rounds 3 -s -m com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -fa : run the fetch-and-add implementation of k-core
//     -nb : the number of buckets to use in the bucketing implementation

#include "KCore.h"
#include "ligra.h"

template <class vertex>
void KCore_runner(graph<vertex>& GA, commandLine P) {
  size_t num_buckets = P.getOptionLongValue("-nb", 16);
  if (num_buckets != (1 << pbbs::log2_up(num_buckets))) {
    cout << "Number of buckets must be a power of two." << endl;
    exit(-1);
  }
  assert(P.getOption("-s"));

  // runs the fetch-and-add based implementation if set.
  bool fa = P.getOption("-fa");
  auto cores = (fa) ? KCore_FA(GA, num_buckets) : KCore(GA, num_buckets);
}

generate_main(KCore_runner, false);
