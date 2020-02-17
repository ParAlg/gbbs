// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#pragma once

#include <limits.h>
#include <stdlib.h>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "bridge.h"
#include "compressed_vertex.h"
#include "edge_map_utils.h"
#include "edge_map_blocked.h"
#include "flags.h"
#include "graph.h"
#include "graph_io.h"
#include "graph_mutation.h"
#include "parse_command_line.h"
#include "vertex.h"
#include "vertex_subset.h"

template <class Data  /* per-vertex data in the emitted vertex_subset */,
          class Graph /* graph type */,
          class VS    /* vertex_subset type */,
          class F     /* edgeMap struct */>
inline vertexSubsetData<Data> edgeMapDense(Graph& GA, VS& vertexSubset, F& f,
                                           const flags fl) {
  using D = std::tuple<bool, Data>;
  size_t n = GA.n;
  auto dense_par = fl & dense_parallel;
  if (should_output(fl)) {
    D* next = pbbslib::new_array_no_init<D>(n);
    auto g = get_emdense_gen<Data>(next);
    parallel_for(
        0, n,
        [&](size_t v) {
          std::get<0>(next[v]) = 0;
          if (f.cond(v)) {
            auto vert = GA.get_vertex(v);
            (fl & in_edges)
                ? vert.decodeOutNghBreakEarly(v, vertexSubset, f, g, dense_par)
                : vert.decodeInNghBreakEarly(v, vertexSubset, f, g, dense_par);
          }
        },
        (fl & fine_parallel) ? 1 : 2048);
    return vertexSubsetData<Data>(n, next);
  } else {
    auto g = get_emdense_nooutput_gen<Data>();
    parallel_for(0, n,
                 [&](size_t v) {
                   if (f.cond(v)) {
                     (fl & in_edges) ? GA.get_vertex(v).decodeOutNghBreakEarly(
                                           v, vertexSubset, f, g, dense_par)
                                     : GA.get_vertex(v).decodeInNghBreakEarly(
                                           v, vertexSubset, f, g, dense_par);
                   }
                 },
                 (fl & fine_parallel) ? 1 : 2048);
    return vertexSubsetData<Data>(n);
  }
}

template <class Data  /* per-vertex data in the emitted vertex_subset */,
          class Graph /* graph type */,
          class VS    /* vertex_subset type */,
          class F     /* edgeMap struct */>
inline vertexSubsetData<Data> edgeMapDenseForward(Graph& GA, VS& vertexSubset, F& f,
                                                  const flags fl) {
  debug(std::cout << "# dense forward" << std::endl;);
  using D = std::tuple<bool, Data>;
  size_t n = GA.n;
  if (should_output(fl)) {
    D* next = pbbslib::new_array_no_init<D>(n);
    auto g = get_emdense_forward_gen<Data>(next);
    par_for(0, n, pbbslib::kSequentialForThreshold,
            [&](size_t i) { std::get<0>(next[i]) = 0; });
    par_for(0, n, 1, [&](size_t i) {
      if (vertexSubset.isIn(i)) {
        (fl & in_edges) ? GA.get_vertex(i).decodeInNgh(i, f, g)
                        : GA.get_vertex(i).decodeOutNgh(i, f, g);
      }
    });
    return vertexSubsetData<Data>(n, next);
  } else {
    auto g = get_emdense_forward_nooutput_gen<Data>();
    par_for(0, n, 1, [&](size_t i) {
      if (vertexSubset.isIn(i)) {
        (fl & in_edges) ? GA.get_vertex(i).decodeInNgh(i, f, g)
                        : GA.get_vertex(i).decodeOutNgh(i, f, g);
      }
    });
    return vertexSubsetData<Data>(n);
  }
}

// Decides on sparse or dense base on number of nonzeros in the active vertices.
template <
    class Data /* data associated with vertices in the output vertex_subset */,
    class G /* graph type */, class VS /* vertex_subset type */,
    class F /* edgeMap struct */>
inline vertexSubsetData<Data> edgeMapData(G& GA, VS& vs, F f,
                                          intT threshold = -1,
                                          const flags& fl = 0) {
  size_t numVertices = GA.n, numEdges = GA.m, m = vs.numNonzeros();
  size_t dense_threshold = threshold;
  if (threshold == -1) dense_threshold = numEdges / 20;
  if (vs.size() == 0) return vertexSubsetData<Data>(numVertices);

  if (vs.isDense && vs.size() > numVertices / 10) {
    return (fl & dense_forward)
               ? edgeMapDenseForward<Data, G, VS, F>(GA, vs, f, fl)
               : edgeMapDense<Data, G, VS, F>(GA, vs, f, fl);
  }

  size_t out_degrees = 0;
  if (vs.out_degrees_set()) {
    out_degrees = vs.get_out_degrees();
  } else {
    vs.toSparse();
    auto degree_f = [&](size_t i) {
      return (fl & in_edges) ? GA.get_vertex(vs.vtx(i)).getInDegree()
                             : GA.get_vertex(vs.vtx(i)).getOutDegree();
    };
    auto degree_im = pbbslib::make_sequence<size_t>(vs.size(), degree_f);
    out_degrees = pbbslib::reduce_add(degree_im);
    vs.set_out_degrees(out_degrees);
  }

  if (out_degrees == 0) return vertexSubsetData<Data>(numVertices);
  if (m + out_degrees > dense_threshold && !(fl & no_dense)) {
    vs.toDense();
    return (fl & dense_forward)
               ? edgeMapDenseForward<Data, G, VS, F>(GA, vs, f, fl)
               : edgeMapDense<Data, G, VS, F>(GA, vs, f, fl);
  } else {
    auto vs_out = edgeMapChunked<Data, G, VS, F>(GA, vs, f, fl);
//    auto vs_out = edgeMapBlocked<Data, G, VS, F>(GA, vs, f, fl);
//    auto vs_out = edgeMapSparse<Data, G, VS, F>(GA, vs, f, fl);
    return vs_out;
  }
}

// Regular edgeMap, where no extra data is stored per vertex.
template <class G /* graph type */, class VS /* vertex_subset type */,
          class F /* edgeMap struct */>
inline vertexSubset edgeMap(G& GA, VS& vs, F f, intT threshold = -1,
                            const flags& fl = 0) {
  return edgeMapData<pbbslib::empty>(GA, vs, f, threshold, fl);
}


/*========================== Vertex Functions =============================*/

template <class F, class VS,
          typename std::enable_if<!std::is_same<VS, vertexSubset>::value,
                                  int>::type = 0>
inline void vertexMap(VS& V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (V.dense()) {
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      if (V.isIn(i)) {
        f(i, V.ithData(i));
      }
    });
  } else {
    par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { f(V.vtx(i), V.vtxData(i)); });
  }
}

template <class VS, class F,
          typename std::enable_if<std::is_same<VS, vertexSubset>::value,
                                  int>::type = 0>
inline void vertexMap(VS& V, F f, size_t granularity=pbbslib::kSequentialForThreshold) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (V.dense()) {
    par_for(0, n, granularity, [&] (size_t i) {
      if (V.isIn(i)) {
        f(i);
      }
    }, granularity);
  } else {
    par_for(0, m, granularity, [&] (size_t i)
                    { f(V.vtx(i)); });
  }
}

// Note: this is the version of vertexMap in which only a subset of the
// input vertexSubset is returned
template <class F>
inline vertexSubset vertexFilter(vertexSubset V, F filter) {
  size_t n = V.numRows();
  V.toDense();
  bool* d_out = pbbslib::new_array_no_init<bool>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { d_out[i] = 0; });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    if (V.d[i]) d_out[i] = filter(i);
  });
  return vertexSubset(n, d_out);
}

template <class F>
inline vertexSubset vertexFilter2(vertexSubset V, F filter) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  bool* bits = pbbslib::new_array_no_init<bool>(m);
  V.toSparse();
  par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintE v = V.vtx(i);
    bits[i] = filter(v);
  });
  auto v_imap = pbbslib::make_sequence<uintE>(m, [&](size_t i) { return V.vtx(i); });
  auto bits_m = pbbslib::make_sequence<bool>(m, [&](size_t i) { return bits[i]; });
  auto out = pbbslib::pack(v_imap, bits_m);
  pbbslib::free_array(bits);
  size_t out_size = out.size();
  return vertexSubset(n, out_size, out.to_array());
}

template <class Data, class F>
inline vertexSubset vertexFilter2(vertexSubsetData<Data> V, F filter) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  bool* bits = pbbslib::new_array_no_init<bool>(m);
  V.toSparse();
  par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    auto t = V.vtxAndData(i);
    bits[i] = filter(std::get<0>(t), std::get<1>(t));
  });
  auto v_imap = pbbslib::make_sequence<uintE>(m, [&](size_t i) { return V.vtx(i); });
  auto bits_m = pbbslib::make_sequence<bool>(m, [&](size_t i) { return bits[i]; });
  auto out = pbbslib::pack(v_imap, bits_m);
  size_t out_size = out.size();
  pbbslib::free_array(bits);
  return vertexSubset(n, out_size, out.to_array());
}

// Adds vertices to a vertexSubset vs.
// Caller must ensure that every v in new_verts is not already in vs
// Note: Mutates the given vertexSubset.
void add_to_vsubset(vertexSubset& vs, uintE* new_verts, uintE num_new_verts);

// cond function that always returns true
inline bool cond_true(intT d) { return 1; }

// Sugar to pass in a single f and get a struct suitable for edgeMap.
template <class W, class F>
struct EdgeMap_F {
  F f;
  EdgeMap_F(F _f) : f(_f) {}
  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    return f(s, d, wgh);
  }

  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    return f(s, d, wgh);
  }

  inline bool cond(const uintE& d) const { return true; }
};

template <class W, class F>
inline EdgeMap_F<W, F> make_em_f(F f) {
  return EdgeMap_F<W, F>(f);
}


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

  cpu_stats(
      double ipc,
      size_t total_cycles,
      double l2_hit_ratio,
      double l3_hit_ratio,
      size_t l2_misses,
      size_t l2_hits,
      size_t l3_misses,
      size_t l3_hits,
      size_t bytes_read,
      size_t bytes_written,
      double total_time,
      size_t num_rounds);

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

void print_pcm_stats(size_t before, size_t after, size_t rounds, double elapsed);
void pcm_init();
#ifdef USE_PCM_LIB

#include "cpucounters.h"

inline auto get_pcm_state() { return getSystemCounterState(); }
#else
inline auto get_pcm_state() { return (size_t)1; }
#endif

#define run_app(G, APP, rounds)                                      \
  auto before_state = get_pcm_state();                               \
  timer st;                                                          \
  double total_time = 0.0;                                           \
  for (size_t r = 0; r < rounds; r++) {                              \
    total_time += APP(G, P);                                         \
  }                                                                  \
  auto time_per_iter = total_time / rounds;                          \
  std::cout << "# time per iter: " << time_per_iter << "\n";           \
  auto after_state = get_pcm_state();                                \
  print_pcm_stats(before_state, after_state, rounds, time_per_iter); \
  G.del();


/* Macro to generate binary for graph applications that read a graph (either
 * asymmetric or symmetric) and transform it into a COO (edge-array)
 * representation for the algorithm. This is currently only used to measure
 * the performance of CSR vs. COO in the graph connectivity benchmark. */
#define generate_coo_main(APP, mutates)                                        \
  int main(int argc, char* argv[]) {                                           \
    commandLine P(argc, argv, " [-s] <inFile>");                               \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool mmapcopy = mutates;                                                   \
    debug(std::cout << "# mmapcopy = " << mmapcopy << "\n";);                    \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    pcm_init();                                                                \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs_io::read_compressed_symmetric_graph<pbbslib::empty>(     \
            iFile, mmap, mmapcopy);                                            \
        alloc_init(G);                                                         \
        auto G_coo = to_edge_array<pbbslib::empty>(G);                         \
        run_app(G_coo, APP, rounds)                                            \
      } else {                                                                 \
        auto G = gbbs_io::read_compressed_asymmetric_graph<pbbslib::empty>(    \
            iFile, mmap, mmapcopy);                                            \
        alloc_init(G);                                                         \
        auto G_coo = to_edge_array<pbbslib::empty>(G);                         \
        run_app(G_coo, APP, rounds)                                            \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G =                                                               \
            gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);             \
        alloc_init(G);                                                         \
        auto G_coo = to_edge_array<pbbslib::empty>(G);                         \
        run_app(G_coo, APP, rounds)                                            \
      } else {                                                                 \
        auto G =                                                               \
            gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap);            \
        alloc_init(G);                                                         \
        auto G_coo = to_edge_array<pbbslib::empty>(G);                         \
        run_app(G_coo, APP, rounds)                                            \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for graph applications that read a graph (either
 * asymmetric or symmetric) and transform it into a COO (edge-array)
 * representation for the algorithm. This is currently only used to measure
 * the performance of CSR vs. COO in the graph connectivity benchmark. */
#define generate_coo_once_main(APP, mutates)                                   \
  int main(int argc, char* argv[]) {                                           \
    commandLine P(argc, argv, " [-s] <inFile>");                               \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool mmapcopy = mutates;                                                   \
    debug(std::cout << "# mmapcopy = " << mmapcopy << "\n";);                    \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    pcm_init();                                                                \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs_io::read_compressed_symmetric_graph<pbbslib::empty>(     \
            iFile, mmap, mmapcopy);                                            \
        alloc_init(G);                                                         \
        auto G_coo = to_edge_array<pbbslib::empty>(G);                         \
        run_app(G_coo, APP, 1)                                                 \
      } else {                                                                 \
        auto G = gbbs_io::read_compressed_asymmetric_graph<pbbslib::empty>(    \
            iFile, mmap, mmapcopy);                                            \
        alloc_init(G);                                                         \
        auto G_coo = to_edge_array<pbbslib::empty>(G);                         \
        run_app(G_coo, APP, 1)                                                 \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G =                                                               \
            gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);             \
        alloc_init(G);                                                         \
        auto G_coo = to_edge_array<pbbslib::empty>(G);                         \
        run_app(G_coo, APP, 1)                                                 \
      } else {                                                                 \
        auto G =                                                               \
            gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap);            \
        alloc_init(G);                                                         \
        auto G_coo = to_edge_array<pbbslib::empty>(G);                         \
        run_app(G_coo, APP, 1)                                                 \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for unweighted graph applications that can ingest only
 * either symmetric or asymmetric graph inputs */
#define generate_main(APP, mutates)                                            \
  int main(int argc, char* argv[]) {                                           \
    commandLine P(argc, argv, " [-s] <inFile>");                               \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool mmapcopy = mutates;                                                   \
    debug(std::cout << "# mmapcopy = " << mmapcopy << "\n";);                    \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    pcm_init();                                                                \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs_io::read_compressed_symmetric_graph<pbbslib::empty>(     \
            iFile, mmap, mmapcopy);                                            \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      } else {                                                                 \
        auto G = gbbs_io::read_compressed_asymmetric_graph<pbbslib::empty>(    \
            iFile, mmap, mmapcopy);                                 \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G =                                                               \
            gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);  \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      } else {                                                                 \
        auto G =                                                               \
            gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap); \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for unweighted graph applications that can ingest only
 * symmetric graph inputs */
#define generate_symmetric_main(APP, mutates)                                  \
  int main(int argc, char* argv[]) {                                           \
    commandLine P(argc, argv, " [-s] <inFile>");                               \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool mmapcopy = mutates;                                                   \
    if (!symmetric) { \
      cout << "# The application expects the input graph to be symmetric (-s flag)." << endl; \
      cout << "# Please run on a symmetric input." << endl; \
    } \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    pcm_init();                                                                \
    if (compressed) {                                                          \
      auto G = gbbs_io::read_compressed_symmetric_graph<pbbslib::empty>(       \
          iFile, mmap, mmapcopy);                                              \
      alloc_init(G);                                                           \
      run_app(G, APP, rounds)                                                  \
    } else {                                                                   \
        auto G =                                                               \
            gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);             \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
    }                                                                          \
  }

/* Macro to generate binary for unweighted graph applications that can ingest only
 * symmetric graph inputs */
#define generate_symmetric_once_main(APP, mutates)                                  \
  int main(int argc, char* argv[]) {                                           \
    commandLine P(argc, argv, " [-s] <inFile>");                               \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool mmapcopy = mutates;                                                   \
    if (!symmetric) { \
      cout << "# The application expects the input graph to be symmetric (-s flag)." << endl; \
      cout << "# Please run on a symmetric input." << endl; \
    } \
    pcm_init();                                                                \
    if (compressed) {                                                          \
      auto G = gbbs_io::read_compressed_symmetric_graph<pbbslib::empty>(       \
          iFile, mmap, mmapcopy);                                              \
      alloc_init(G);                                                           \
      run_app(G, APP, 1)                                                       \
    } else {                                                                   \
        auto G =                                                               \
            gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);             \
        alloc_init(G);                                                         \
        run_app(G, APP, 1)                                                     \
    }                                                                          \
  }

/* Macro to generate binary for weighted graph applications that can ingest
 * either symmetric or asymmetric graph inputs */
#define generate_weighted_main(APP, mutates)                                   \
  int main(int argc, char* argv[]) {                                           \
    commandLine P(argc, argv, " [-s] <inFile>");                               \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool mmapcopy = mutates;                                                   \
    debug(std::cout << "# mmapcopy = " << mmapcopy << "\n";);                    \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    pcm_init();                                                                \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs_io::read_compressed_symmetric_graph<intE>(               \
            iFile, mmap, mmapcopy);                                            \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      } else {                                                                 \
        auto G = gbbs_io::read_compressed_asymmetric_graph<intE>(              \
            iFile, mmap, mmapcopy);                                            \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G =                                                               \
            gbbs_io::read_weighted_symmetric_graph(iFile, mmap);               \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      } else {                                                                 \
        auto G =                                                               \
            gbbs_io::read_weighted_asymmetric_graph(iFile, mmap);              \
        alloc_init(G);                                                         \
        run_app(G, APP, rounds)                                                \
      }                                                                        \
    }                                                                          \
  }

/* Macro to generate binary for weighted graph applications that can ingest
 * only symmetric graph inputs */
#define generate_symmetric_weighted_main(APP, mutates)                         \
  int main(int argc, char* argv[]) {                                           \
    commandLine P(argc, argv, " [-s] <inFile>");                               \
    char* iFile = P.getArgument(0);                                            \
    debug(bool symmetric = P.getOptionValue("-s");                             \
    assert(symmetric););                                                       \
    bool compressed = P.getOptionValue("-c");                                  \
    bool mmap = P.getOptionValue("-m");                                        \
    bool mmapcopy = mutates;                                                   \
    debug(std::cout << "# mmapcopy = " << mmapcopy << "\n";);                    \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    pcm_init();                                                                \
    if (compressed) {                                                          \
      auto G = gbbs_io::read_compressed_symmetric_graph<intE>(                 \
          iFile, mmap, mmapcopy);                                              \
      alloc_init(G);                                                           \
      run_app(G, APP, rounds)                                                  \
    } else {                                                                   \
      auto G =                                                                 \
          gbbs_io::read_weighted_symmetric_graph(iFile, mmap);                 \
      alloc_init(G);                                                           \
      run_app(G, APP, rounds)                                                  \
    }                                                                          \
  }
