#pragma once

#include "immutable_graph.h"
#include "traversable_graph.h"

#include "gbbs/benchmark.h"

namespace aspen {

template <class Graph>
auto symmetric_graph_from_static_graph(Graph& GA) {
  using W = typename Graph::weight_type;
  using inner_graph = symmetric_graph<W>;
  using outer_graph = traversable_graph<inner_graph>;
  timer build_t;
  build_t.start();
  auto G = outer_graph(GA);
  build_t.stop();
  build_t.reportTotal("Aspen: build time");
  return G;
}

}  // namespace aspen

/* Macro to generate binary for unweighted graph applications that can ingest
 * only
 * symmetric graph inputs */
#define generate_symmetric_aspen_main(APP, mutates)                           \
  int main(int argc, char* argv[]) {                                          \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                        \
    char* iFile = P.getArgument(0);                                           \
    bool symmetric = P.getOptionValue("-s");                                  \
    bool compressed = P.getOptionValue("-c");                                 \
    bool mmap = P.getOptionValue("-m");                                       \
    bool mmapcopy = mutates;                                                  \
    if (!symmetric) {                                                         \
      std::cout                                                               \
          << "# The application expects the input graph to be symmetric (-s " \
             "flag)."                                                         \
          << std::endl;                                                       \
      std::cout << "# Please run on a symmetric input." << std::endl;         \
    }                                                                         \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                       \
    gbbs::pcm_init();                                                         \
    if (compressed) {                                                         \
      auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(   \
          iFile, mmap, mmapcopy);                                             \
      auto AG = aspen::symmetric_graph_from_static_graph(G);                  \
      run_app(AG, APP, rounds)                                                \
    } else {                                                                  \
      auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);   \
      auto AG = aspen::symmetric_graph_from_static_graph(G);                  \
      run_app(AG, APP, rounds)                                                \
    }                                                                         \
  }
