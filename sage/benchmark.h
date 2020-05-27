#pragma once

#include "sage_io.h"

/* Macro to generate binary for unweighted graph applications that can ingest
 * either symmetric or asymmetric graph inputs */
#define generate_sage_main(APP)                                                \
  int main(int argc, char* argv[]) {                                           \
    commandLine P(argc, argv, " [-s] <inFile>");                               \
    char* f1 = P.getOptionValue("-f1");                                        \
    char* f2 = P.getOptionValue("-f2");                                        \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    pcm_init();                                                                \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G =                                                               \
            sage_io::read_compressed_symmetric_graph<pbbslib::empty>(f1, f2);  \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      } else {                                                                 \
        auto G =                                                               \
            sage_io::read_compressed_asymmetric_graph<pbbslib::empty>(f1, f2); \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = sage_io::read_symmetric_binary_graph<pbbslib::empty>(f1, f2); \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      } else {                                                                 \
        auto G =                                                               \
            sage_io::read_asymmetric_binary_graph<pbbslib::empty>(f1, f2);     \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      }                                                                        \
    }                                                                          \
    alloc_finish();                                                            \
  }
