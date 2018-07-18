// Usage:
// > numactl -i all ./Coloring -s -m clueweb_sym.bytepda
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//   optional:
//     -c : indicate that the graph should be mmap'd
//     -m : indicate that the graph is compressed
//     -lf : use the LF (largest degree first) herustic
//     -stats : output statistics on the resulting coloring
//     -verify : verify that the algorithm produced a valid coloring
//
// The default heuristic used is LLF, which provably achieves polynomial
// parallelism (see "Ordering Heuristics for Parallel Graph Coloring" by
// Hasenplaugh et al.)

#include "Coloring.h"

#include <fstream>
#include <iostream>

template <class vertex>
void Coloring_runner(graph<vertex>& GA, commandLine P) {
  bool runLF = P.getOption("-lf");
  auto colors = Coloring(GA, runLF);
  if (P.getOption("-stats")) {
    cout << "num_colors = " << pbbs::reduce_max(colors) << endl;
  }
  if (P.getOption("-verify)")) {
    verify_coloring(GA, colors);
  }
}

generate_main(Coloring_runner, false);
