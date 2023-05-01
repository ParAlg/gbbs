// Utilities for creating main files that read graphs.
#pragma once

#include "assert.h"
#include "graph_io.h"

#define run_app(G, APP, mutates, rounds)    \
  double total_time = 0.0;                  \
  for (size_t r = 0; r < rounds; r++) {     \
    if (mutates) {                          \
      auto G_copy = G;                      \
      total_time += APP(G_copy, P);         \
    } else {                                \
      total_time += APP(G, P);              \
    }                                       \
  }                                         \
  auto time_per_iter = total_time / rounds; \
  std::cout << "# time per iter: " << time_per_iter << "\n";

/* Macro to generate binary for graph applications that read a graph (either
 * asymmetric or symmetric) and transform it into a COO (edge-array)
 * representation for the algorithm. This is currently only used to measure
 * the performance of CSR vs. COO in the graph connectivity benchmark. */
#define generate_coo_main(APP, mutates)                                        \
  int main(int argc, char* argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(  \
            iFile, mmap);                                                      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, rounds)                                   \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>( \
            iFile, mmap);                                                      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, rounds)                                   \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap,   \
                                                                binary);       \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, rounds)                                   \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap,  \
                                                                 binary);      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, rounds)                                   \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for graph applications that read a graph (either
 * asymmetric or symmetric) and transform it into a COO (edge-array)
 * representation for the algorithm. This is currently only used to measure
 * the performance of CSR vs. COO in the graph connectivity benchmark. */
#define generate_coo_once_main(APP, mutates)                                   \
  int main(int argc, char* argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(  \
            iFile, mmap);                                                      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, 1)                                        \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>( \
            iFile, mmap);                                                      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, 1)                                        \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap,   \
                                                                binary);       \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, 1)                                        \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap,  \
                                                                 binary);      \
        auto G_coo = to_edge_array<gbbs::empty>(G);                            \
        run_app(G_coo, APP, mutates, 1)                                        \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for unweighted graph applications that can ingest
 * only
 * either symmetric or asymmetric graph inputs */
#define generate_main(APP, mutates)                                            \
  int main(int argc, char* argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool binary = P.getOptionValue("-b");                                      \
    bool mmap = P.getOptionValue("-m");                                        \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(  \
            iFile, mmap);                                                      \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>( \
            iFile, mmap);                                                      \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap,   \
                                                                binary);       \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap,  \
                                                                 binary);      \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for unweighted graph applications that can ingest
 * only asymmetric graph inputs */
#define generate_asymmetric_main(APP, mutates)                               \
  int main(int argc, char* argv[]) {                                         \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                       \
    char* iFile = P.getArgument(0);                                          \
    bool compressed = P.getOptionValue("-c");                                \
    bool mmap = P.getOptionValue("-m");                                      \
    bool binary = P.getOptionValue("-b");                                    \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                      \
    if (compressed) {                                                        \
      auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>( \
          iFile, mmap);                                                      \
      run_app(G, APP, mutates, rounds)                                       \
    } else {                                                                 \
      auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap,  \
                                                               binary);      \
      run_app(G, APP, mutates, rounds)                                       \
    }                                                                        \
  }

/* Macro to generate binary for unweighted graph applications that can ingest
 * only
 * symmetric graph inputs */
#define generate_symmetric_main(APP, mutates)                                  \
  int main(int argc, char* argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    if (!symmetric) {                                                          \
      std::cout                                                                \
          << "# The application expects the input graph to be symmetric (-s "  \
             "flag)."                                                          \
          << std::endl;                                                        \
      std::cout << "# Please run on a symmetric input." << std::endl;          \
    }                                                                          \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(    \
          iFile, mmap);                                                        \
      run_app(G, APP, mutates, rounds)                                         \
    } else {                                                                   \
      auto G =                                                                 \
          gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap, binary); \
      run_app(G, APP, mutates, rounds)                                         \
    }                                                                          \
  }

/* Macro to generate binary for unweighted graph applications that can ingest
 * only
 * symmetric graph inputs */
#define generate_symmetric_once_main(APP, mutates)                             \
  int main(int argc, char* argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    if (!symmetric) {                                                          \
      std::cout                                                                \
          << "# The application expects the input graph to be symmetric (-s "  \
             "flag)."                                                          \
          << std::endl;                                                        \
      std::cout << "# Please run on a symmetric input." << std::endl;          \
    }                                                                          \
    if (compressed) {                                                          \
      auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(    \
          iFile, mmap);                                                        \
      run_app(G, APP, mutates, 1)                                              \
    } else {                                                                   \
      auto G =                                                                 \
          gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap, binary); \
      run_app(G, APP, mutates, 1)                                              \
    }                                                                          \
  }

/* Macro to generate binary for weighted graph applications that can ingest
 * either symmetric or asymmetric graph inputs */
#define generate_weighted_main(APP, mutates)                                  \
  int main(int argc, char* argv[]) {                                          \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                        \
    char* iFile = P.getArgument(0);                                           \
    bool symmetric = P.getOptionValue("-s");                                  \
    bool compressed = P.getOptionValue("-c");                                 \
    bool mmap = P.getOptionValue("-m");                                       \
    bool binary = P.getOptionValue("-b");                                     \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                       \
    if (compressed) {                                                         \
      if (symmetric) {                                                        \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::intE>(  \
            iFile, mmap);                                                     \
        run_app(G, APP, mutates, rounds)                                      \
      } else {                                                                \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::intE>( \
            iFile, mmap);                                                     \
        run_app(G, APP, mutates, rounds)                                      \
      }                                                                       \
    } else {                                                                  \
      if (symmetric) {                                                        \
        auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<gbbs::intE>(    \
            iFile, mmap, binary);                                             \
        run_app(G, APP, mutates, rounds)                                      \
      } else {                                                                \
        auto G = gbbs::gbbs_io::read_weighted_asymmetric_graph<gbbs::intE>(   \
            iFile, mmap, binary);                                             \
        run_app(G, APP, mutates, rounds)                                      \
      }                                                                       \
    }                                                                         \
  }

/* Macro to generate binary for weighted graph applications that can ingest
 * either symmetric or asymmetric graph inputs */
#define generate_float_main(APP, mutates)                                      \
  int main(int argc, char* argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool binary = P.getOptionValue("-b");                                      \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<float>(iFile,  \
                                                                       mmap);  \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<float>(iFile, \
                                                                        mmap); \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<float>(          \
            iFile, mmap, binary);                                              \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_weighted_asymmetric_graph<float>(         \
            iFile, mmap, binary);                                              \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for weighted graph applications that can ingest
 * only symmetric graph inputs */
#define generate_symmetric_weighted_main(APP, mutates)                     \
  int main(int argc, char* argv[]) {                                       \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                     \
    char* iFile = P.getArgument(0);                                        \
    gbbs_debug(bool symmetric = P.getOptionValue("-s"); assert(symmetric););    \
    bool compressed = P.getOptionValue("-c");                              \
    bool mmap = P.getOptionValue("-m");                                    \
    bool binary = P.getOptionValue("-b");                                  \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                    \
    if (compressed) {                                                      \
      auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::intE>( \
          iFile, mmap);                                                    \
      run_app(G, APP, mutates, rounds)                                     \
    } else {                                                               \
      auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<gbbs::intE>(   \
          iFile, mmap, binary);                                            \
      run_app(G, APP, mutates, rounds)                                     \
    }                                                                      \
  }

/* Macro to generate binary for floating-point weighted graph applications that
 * can ingest only symmetric graph inputs */
#define generate_symmetric_float_weighted_main(APP)                     \
  int main(int argc, char* argv[]) {                                    \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                  \
    char* iFile = P.getArgument(0);                                     \
    gbbs_debug(bool symmetric = P.getOptionValue("-s"); assert(symmetric);); \
    bool compressed = P.getOptionValue("-c");                           \
    bool mmap = P.getOptionValue("-m");                                 \
    bool binary = P.getOptionValue("-b");                               \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                 \
    if (compressed) {                                                   \
      ABORT("Graph compression not yet implemented for float weights"); \
    } else {                                                            \
      auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<float>(     \
          iFile, mmap, binary);                                         \
      run_app(G, APP, mutates, rounds)                                  \
    }                                                                   \
  }
