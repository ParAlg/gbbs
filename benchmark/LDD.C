// Usage:
// numactl -i all ./LDD -beta 0.2 -rounds 3 -s -m com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -fa : run the fetch-and-add implementation of k-core
//     -nb : the number of buckets to use in the bucketing implementation

#include "LDD.h"
#include "ligra.h"

template <class vertex>
void LDD_runner(graph<vertex>& GA, commandLine P) {
  double beta = P.getOptionDoubleValue("-beta", 0.2);
  bool permute = P.getOption("-permute");
  assert(P.getOption("-s"));
  auto ldd = LDD(GA, beta, permute, false);
  if (P.getOption("-stats")) {
    ldd_utils::num_clusters(ldd);
    ldd_utils::num_intercluster_edges(GA, ldd);
  }
}

generate_main(LDD_runner, false);
