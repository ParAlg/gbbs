// Usage:
// numactl -i all ./BellmanFord -src 10012 -s -m -rounds 3 twitter_wgh_SJ
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//     -w: indicate that the graph is weighted
//   optional:
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -specfor : run the speculative_for (union-find) based algorithm from pbbs
//     -largemem : set the sampling thresholds to utilize less memory
//
// Note: in our experiments we set -largemem when running MST on the weighted
// hyperlink2012 graph.

#include "MST.h"
#include "ligra.h"

template <class vertex>
void MST_runner(graph<vertex>& GA, commandLine P) {
  bool spec_for = P.getOption("-specfor");
  bool largemem = P.getOption("-largemem");
  timer mst_t; mst_t.start();
  if (spec_for) {
    MST_spec_for::MST(GA);
  } else {
    MST_boruvka::MST(GA, largemem);
  }
  mst_t.stop(); mst_t.reportTotal("MST time");
  // MST mutates the underlying graph (unless it is copied, which we don't do to
  // prevent memory issues), so we make sure the algorithm is run exactly once.
  exit(0);
}

generate_main(MST_runner, true);
