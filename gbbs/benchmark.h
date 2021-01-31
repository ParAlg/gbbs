
// Utilities for creating main files that read graphs, and other benchmark tools
// (e.g., profiling)
#pragma once

#include "assert.h"
#include "graph_io.h"

namespace gbbs {
/* Aggregate metrics for a repeated experiment, repeated num_rounds times. */
struct cpu_stats {
  double ipc;
  size_t total_cycles;
  double l2_hit_ratio;
  double l3_hit_ratio;
  size_t l2_misses;
  size_t l2_hits;
  size_t l3_misses;
  size_t l3_hits;
  size_t bytes_read;
  size_t bytes_written;
  double total_time;
  size_t num_rounds;

  cpu_stats();

  cpu_stats(double ipc, size_t total_cycles, double l2_hit_ratio,
            double l3_hit_ratio, size_t l2_misses, size_t l2_hits,
            size_t l3_misses, size_t l3_hits, size_t bytes_read,
            size_t bytes_written, double total_time, size_t num_rounds);

  double get_ipc();
  size_t get_total_cycles();
  double get_l2_hit_ratio();
  double get_l3_hit_ratio();
  size_t get_l2_misses();
  size_t get_l2_hits();
  size_t get_l3_misses();
  size_t get_l3_hits();
  double get_throughput();
};

void print_pcm_stats(size_t before, size_t after, size_t rounds,
                     double elapsed);
void pcm_init();
#ifdef USE_PCM_LIB

#include "cpucounters.h"

inline auto get_pcm_state() { return getSystemCounterState(); }
#else
inline auto get_pcm_state() { return (size_t)1; }
#endif
}  // namespace gbbs

#define run_app(G, APP, rounds)                                            \
  auto before_state = gbbs::get_pcm_state();                               \
  gbbs::timer st;                                                          \
  double total_time = 0.0;                                                 \
  for (size_t r = 0; r < rounds; r++) {                                    \
    total_time += APP(G, P);                                               \
  }                                                                        \
  auto time_per_iter = total_time / rounds;                                \
  std::cout << "# time per iter: " << time_per_iter << "\n";               \
  auto after_state = gbbs::get_pcm_state();                                \
  gbbs::print_pcm_stats(before_state, after_state, rounds, time_per_iter);

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
    bool mmapcopy = mutates;                                                   \
    debug(std::cout << "# mmapcopy = " << mmapcopy << "\n";);                  \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    gbbs::pcm_init();                                                          \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G =                                                               \
            gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(    \
                iFile, mmap, mmapcopy);                                        \
        gbbs::alloc_init(G);                                                   \
        auto G_coo = to_edge_array<gbbs::empty>(G);                         \
        run_app(G_coo, APP, rounds)                                            \
      } else {                                                                 \
        auto G =                                                               \
            gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>(   \
                iFile, mmap, mmapcopy);                                        \
        gbbs::alloc_init(G);                                                   \
        auto G_coo = to_edge_array<gbbs::empty>(G);                         \
        run_app(G_coo, APP, rounds)                                            \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);  \
        gbbs::alloc_init(G);                                                   \
        auto G_coo = to_edge_array<gbbs::empty>(G);                         \
        run_app(G_coo, APP, rounds)                                            \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap); \
        gbbs::alloc_init(G);                                                   \
        auto G_coo = to_edge_array<gbbs::empty>(G);                         \
        run_app(G_coo, APP, rounds)                                            \
      }                                                                        \
    }                                                                          \
    gbbs::alloc_finish();                                                      \
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
    bool mmapcopy = mutates;                                                   \
    debug(std::cout << "# mmapcopy = " << mmapcopy << "\n";);                  \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    gbbs::pcm_init();                                                          \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G =                                                               \
            gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(    \
                iFile, mmap, mmapcopy);                                        \
        gbbs::alloc_init(G);                                                   \
        auto G_coo = to_edge_array<gbbs::empty>(G);                         \
        run_app(G_coo, APP, 1)                                                 \
      } else {                                                                 \
        auto G =                                                               \
            gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>(   \
                iFile, mmap, mmapcopy);                                        \
        gbbs::alloc_init(G);                                                   \
        auto G_coo = to_edge_array<gbbs::empty>(G);                         \
        run_app(G_coo, APP, 1)                                                 \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);  \
        gbbs::alloc_init(G);                                                   \
        auto G_coo = to_edge_array<gbbs::empty>(G);                         \
        run_app(G_coo, APP, 1)                                                 \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap); \
        gbbs::alloc_init(G);                                                   \
        auto G_coo = to_edge_array<gbbs::empty>(G);                         \
        run_app(G_coo, APP, 1)                                                 \
      }                                                                        \
    }                                                                          \
    gbbs::alloc_finish();                                                      \
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
    bool mmap = P.getOptionValue("-m");                                        \
    bool mmapcopy = mutates;                                                   \
    debug(std::cout << "# mmapcopy = " << mmapcopy << "\n";);                  \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    gbbs::pcm_init();                                                          \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G =                                                               \
            gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(    \
                iFile, mmap, mmapcopy);                                        \
        gbbs::alloc_init(G);                                                   \
        run_app(G, APP, rounds)                                                \
      } else {                                                                 \
        auto G =                                                               \
            gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>(   \
                iFile, mmap, mmapcopy);                                        \
        gbbs::alloc_init(G);                                                   \
        run_app(G, APP, rounds)                                                \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);  \
        gbbs::alloc_init(G);                                                   \
        run_app(G, APP, rounds)                                                \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap); \
        gbbs::alloc_init(G);                                                   \
        run_app(G, APP, rounds)                                                \
      }                                                                        \
    }                                                                          \
    gbbs::alloc_finish();                                                      \
  }

/* Macro to generate binary for unweighted graph applications that can ingest
 * only asymmetric graph inputs */
#define generate_asymmetric_main(APP, mutates)                               \
  int main(int argc, char* argv[]) {                                         \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                       \
    char* iFile = P.getArgument(0);                                          \
    bool compressed = P.getOptionValue("-c");                                \
    bool mmap = P.getOptionValue("-m");                                      \
    bool mmapcopy = mutates;                                                 \
    assert(!symmetric);                                                      \
    debug(std::cout << "# mmapcopy = " << mmapcopy << "\n";);                \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                      \
    gbbs::pcm_init();                                                        \
    if (compressed) {                                                        \
      auto G =                                                               \
          gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>(   \
              iFile, mmap, mmapcopy);                                        \
      gbbs::alloc_init(G);                                                   \
      run_app(G, APP, rounds)                                                \
    } else {                                                                 \
      auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap); \
      gbbs::alloc_init(G);                                                   \
      run_app(G, APP, rounds)                                                \
    }                                                                        \
    gbbs::alloc_finish();                                                    \
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
    bool mmapcopy = mutates;                                                   \
    if (!symmetric) {                                                          \
      std::cout                                                                \
          << "# The application expects the input graph to be symmetric (-s "  \
             "flag)."                                                          \
          << std::endl;                                                        \
      std::cout << "# Please run on a symmetric input." << std::endl;          \
    }                                                                          \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    gbbs::pcm_init();                                                          \
    if (compressed) {                                                          \
      auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>( \
          iFile, mmap, mmapcopy);                                              \
      gbbs::alloc_init(G);                                                     \
      run_app(G, APP, rounds)                                                  \
    } else {                                                                   \
      auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);    \
      gbbs::alloc_init(G);                                                     \
      run_app(G, APP, rounds)                                                  \
    }                                                                          \
    gbbs::alloc_finish();                                                      \
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
    bool mmapcopy = mutates;                                                   \
    if (!symmetric) {                                                          \
      std::cout                                                                \
          << "# The application expects the input graph to be symmetric (-s "  \
             "flag)."                                                          \
          << std::endl;                                                        \
      std::cout << "# Please run on a symmetric input." << std::endl;          \
    }                                                                          \
    gbbs::pcm_init();                                                          \
    if (compressed) {                                                          \
      auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>( \
          iFile, mmap, mmapcopy);                                              \
      gbbs::alloc_init(G);                                                     \
      run_app(G, APP, 1)                                                       \
    } else {                                                                   \
      auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);    \
      gbbs::alloc_init(G);                                                     \
      run_app(G, APP, 1)                                                       \
    }                                                                          \
    gbbs::alloc_finish();                                                      \
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
    bool mmapcopy = mutates;                                                  \
    debug(std::cout << "# mmapcopy = " << mmapcopy << "\n";);                 \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                       \
    gbbs::pcm_init();                                                         \
    if (compressed) {                                                         \
      if (symmetric) {                                                        \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<float>(  \
            iFile, mmap, mmapcopy);                                           \
        gbbs::alloc_init(G);                                                  \
        run_app(G, APP, rounds)                                               \
      } else {                                                                \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<float>( \
            iFile, mmap, mmapcopy);                                           \
        gbbs::alloc_init(G);                                                  \
        run_app(G, APP, rounds)                                               \
      }                                                                       \
    } else {                                                                  \
      if (symmetric) {                                                        \
        auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<float>(    \
            iFile, mmap);                                                     \
        gbbs::alloc_init(G);                                                  \
        run_app(G, APP, rounds)                                               \
      } else {                                                                \
        auto G = gbbs::gbbs_io::read_weighted_asymmetric_graph<float>(   \
            iFile, mmap);                                                     \
        gbbs::alloc_init(G);                                                  \
        run_app(G, APP, rounds)                                               \
      }                                                                       \
    }                                                                         \
    gbbs::alloc_finish();                                                     \
  }

/* Macro to generate binary for weighted graph applications that can ingest
 * only symmetric graph inputs */
#define generate_symmetric_weighted_main(APP, mutates)                         \
  int main(int argc, char* argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char* iFile = P.getArgument(0);                                            \
    debug(bool symmetric = P.getOptionValue("-s"); assert(symmetric););        \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool mmapcopy = mutates;                                                   \
    debug(std::cout << "# mmapcopy = " << mmapcopy << "\n";);                  \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    gbbs::pcm_init();                                                          \
    if (compressed) {                                                          \
      auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<float>(     \
          iFile, mmap, mmapcopy);                                              \
      gbbs::alloc_init(G);                                                     \
      run_app(G, APP, rounds)                                                  \
    } else {                                                                   \
      auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<float>(iFile, \
                                                                        mmap); \
      gbbs::alloc_init(G);                                                     \
      run_app(G, APP, rounds)                                                  \
    }                                                                          \
    gbbs::alloc_finish();                                                      \
  }

/* Macro to generate binary for floating-point weighted graph applications that
 * can ingest only symmetric graph inputs */
#define generate_symmetric_float_weighted_main(APP)                         \
  int main(int argc, char* argv[]) {                                        \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                      \
    char* iFile = P.getArgument(0);                                         \
    debug(bool symmetric = P.getOptionValue("-s"); assert(symmetric););     \
    bool compressed = P.getOptionValue("-c");                               \
    bool mmap = P.getOptionValue("-m");                                     \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                     \
    gbbs::pcm_init();                                                       \
    if (compressed) {                                                       \
      ABORT("Graph compression not yet implemented for float weights");     \
    } else {                                                                \
      auto G =                                                              \
          gbbs::gbbs_io::read_weighted_symmetric_graph<float>(iFile, mmap); \
      gbbs::alloc_init(G);                                                  \
      run_app(G, APP, rounds)                                               \
    }                                                                       \
    gbbs::alloc_finish();                                                   \
  }
