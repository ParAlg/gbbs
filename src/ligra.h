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

#include "IO.h"
#include "bridge.h"
#include "compressed_vertex.h"
#include "edge_map_utils.h"
#include "edge_map_blocked.h"
#include "flags.h"
#include "graph.h"
#include "packed_graph.h"
#include "parse_command_line.h"
#include "vertex.h"
#include "vertex_subset.h"

template <
    class data /* data associated with vertices in the output vertex_subset */,
    class G /* graph type */, class VS /* vertex_subset type */,
    class F /* edgeMap struct */>
inline vertexSubsetData<data> edgeMapDense(G& GA, VS& vertexSubset, F& f,
                                           const flags fl) {
  using D = std::tuple<bool, data>;
  size_t n = GA.n;
  auto dense_par = fl & dense_parallel;
  if (should_output(fl)) {
    D* next = pbbslib::new_array_no_init<D>(n);
    auto g = get_emdense_gen<data>(next);
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
    return vertexSubsetData<data>(n, next);
  } else {
    auto g = get_emdense_nooutput_gen<data>();
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
    return vertexSubsetData<data>(n);
  }
}

template <
    class data /* data associated with vertices in the output vertex_subset */,
    class G /* graph type */, class VS /* vertex_subset type */,
    class F /* edgeMap struct */>
inline vertexSubsetData<data> edgeMapDenseForward(G& GA, VS& vertexSubset, F& f,
                                                  const flags fl) {
  debug(std::cout << "dense forward" << std::endl;);
  using D = std::tuple<bool, data>;
  size_t n = GA.n;
  if (should_output(fl)) {
    D* next = pbbslib::new_array_no_init<D>(n);
    auto g = get_emdense_forward_gen<data>(next);
    par_for(0, n, pbbslib::kSequentialForThreshold,
            [&](size_t i) { std::get<0>(next[i]) = 0; });
    par_for(0, n, 1, [&](size_t i) {
      if (vertexSubset.isIn(i)) {
        (fl & in_edges) ? GA.get_vertex(i).decodeInNgh(i, f, g)
                        : GA.get_vertex(i).decodeOutNgh(i, f, g);
      }
    });
    return vertexSubsetData<data>(n, next);
  } else {
    auto g = get_emdense_forward_nooutput_gen<data>();
    par_for(0, n, 1, [&](size_t i) {
      if (vertexSubset.isIn(i)) {
        (fl & in_edges) ? GA.get_vertex(i).decodeInNgh(i, f, g)
                        : GA.get_vertex(i).decodeOutNgh(i, f, g);
      }
    });
    return vertexSubsetData<data>(n);
  }
}


// Decides on sparse or dense base on number of nonzeros in the active vertices.
template <
    class data /* data associated with vertices in the output vertex_subset */,
    class G /* graph type */, class VS /* vertex_subset type */,
    class F /* edgeMap struct */>
inline vertexSubsetData<data> edgeMapData(G& GA, VS& vs, F f,
                                          intT threshold = -1,
                                          const flags& fl = 0) {
  size_t numVertices = GA.n, numEdges = GA.m, m = vs.numNonzeros();
  size_t dense_threshold = threshold;
  if (threshold == -1) dense_threshold = numEdges / 20;
  if (vs.size() == 0) return vertexSubsetData<data>(numVertices);

  if (vs.isDense && vs.size() > numVertices / 10) {
    return (fl & dense_forward)
               ? edgeMapDenseForward<data, G, VS, F>(GA, vs, f, fl)
               : edgeMapDense<data, G, VS, F>(GA, vs, f, fl);
  }

  vs.toSparse();
  auto degree_f = [&](size_t i) {
    return (fl & in_edges) ? GA.get_vertex(vs.vtx(i)).getInDegree()
                           : GA.get_vertex(vs.vtx(i)).getOutDegree();
  };
  auto degree_im = pbbslib::make_sequence<size_t>(vs.size(), degree_f);
  size_t out_degrees = pbbslib::reduce_add(degree_im);
  cout << "out_degree = " << out_degrees << endl;

  if (out_degrees == 0) return vertexSubsetData<data>(numVertices);
  if (m + out_degrees > dense_threshold && !(fl & no_dense)) {
    vs.toDense();
    return (fl & dense_forward)
               ? edgeMapDenseForward<data, G, VS, F>(GA, vs, f, fl)
               : edgeMapDense<data, G, VS, F>(GA, vs, f, fl);
  } else {
    auto vs_out = edgeMapBlocked<data, G, VS, F>(GA, vs, f, fl);
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



// Packs out the adjacency lists of all vertex in vs. A neighbor, ngh, is kept
// in the new adjacency list if p(ngh) is true.
template <class G /* packable graph type */, class P /* predicate function */>
inline vertexSubsetData<uintE> edgeMapPack(G& GA, vertexSubset& vs, P& p,
                                           const flags& fl = 0) {
  using S = std::tuple<uintE, uintE>;
  vs.toSparse();
  size_t m = vs.numNonzeros();
  size_t n = vs.numRows();
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  auto space = sequence<uintT>(m);
  par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintE v = vs.vtx(i);
    space[i] = GA.get_vertex(v).calculateOutTemporarySpaceBytes();
  });

//  size_t total_space = pbbslib::scan_add_inplace(space.slice());
//  std::cout << "packNghs: total space allocated = " << total_space << "\n";
////  sequence<uint8_t> tmp(total_space);
//  uint8_t* tmp = nullptr;
//  if (total_space > 0) {
//    tmp = pbbs::new_array_no_init<uint8_t>(total_space);
//  }
  S* outV;
  if (should_output(fl)) {
    outV = pbbslib::new_array_no_init<S>(vs.size());
    {
      auto for_inner = [&](size_t i) {
        uintE v = vs.vtx(i);
//        uint8_t* tmp_v = tmp + space[i];
        size_t ct = GA.get_vertex(v).packOutNgh(v, p, nullptr);
        outV[i] = std::make_tuple(v, ct);
      };
      par_for(0, m, 1, [&](size_t i) { for_inner(i); });
    }
//    if (total_space > 0) {
//      pbbs::free_array(tmp);
//    }
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    {
      auto for_inner = [&](size_t i) {
        uintE v = vs.vtx(i);
//        uint8_t* tmp_v = tmp + space[i];
        GA.get_vertex(v).packOutNgh(v, p, nullptr);
      };
      par_for(0, m, 1, [&](size_t i) { for_inner(i); });
    }
//    if (total_space > 0) {
//      pbbs::free_array(tmp);
//    }
    return vertexSubsetData<uintE>(n);
  }
  // TODO: update degrees
}


template <class G /* graph type */, class P /* predicate function */>
inline vertexSubsetData<uintE> edgeMapFilter(G& GA, vertexSubset& vs, P& p,
                                             const flags& fl = 0) {
  vs.toSparse();
  size_t m = vs.numNonzeros();
  size_t n = vs.numRows();
  using S = std::tuple<uintE, uintE>;
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  S* outV;
  if (should_output(fl)) {
    outV = pbbslib::new_array_no_init<S>(vs.size());
    parallel_for(0, m,
                 [&](size_t i) {
                   uintE v = vs.vtx(i);
                   size_t ct = GA.get_vertex(v).countOutNgh(v, p);
                   outV[i] = std::make_tuple(v, ct);
                 },
                 1);
  } else {
    parallel_for(0, m,
                 [&](size_t i) {
                   uintE v = vs.vtx(i);
                   GA.get_vertex(v).countOutNgh(v, p);
                 },
                 1);
  }
  if (should_output(fl)) {
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    return vertexSubsetData<uintE>(n);
  }
}

//*****VERTEX FUNCTIONS*****

template <class F, class VS,
          typename std::enable_if<!std::is_same<VS, vertexSubset>::value,
                                  int>::type = 0>
inline void vertexMap(VS& V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (V.dense()) {
    par_for(0, n, pbbslib::kSequentialForThreshold, [&](size_t i) {
      if (V.isIn(i)) {
        f(i, V.ithData(i));
      }
    });
  } else {
    par_for(0, m, pbbslib::kSequentialForThreshold,
            [&](size_t i) { f(V.vtx(i), V.vtxData(i)); });
  }
}

template <class VS, class F,
          typename std::enable_if<std::is_same<VS, vertexSubset>::value,
                                  int>::type = 0>
inline void vertexMap(VS& V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (V.dense()) {
    par_for(0, n, pbbslib::kSequentialForThreshold, [&](size_t i) {
      if (V.isIn(i)) {
        f(i);
      }
    });
  } else {
    par_for(0, m, pbbslib::kSequentialForThreshold,
            [&](size_t i) { f(V.vtx(i)); });
  }
}

// Note: this is the version of vertexMap in which only a subset of the
// input vertexSubset is returned
template <class F>
inline vertexSubset vertexFilter(vertexSubset V, F filter) {
  size_t n = V.numRows();
  V.toDense();
  bool* d_out = pbbslib::new_array_no_init<bool>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold,
          [&](size_t i) { d_out[i] = 0; });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&](size_t i) {
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
  par_for(0, m, pbbslib::kSequentialForThreshold, [&](size_t i) {
    uintE v = V.vtx(i);
    bits[i] = filter(v);
  });
  auto v_imap =
      pbbslib::make_sequence<uintE>(m, [&](size_t i) { return V.vtx(i); });
  auto bits_m =
      pbbslib::make_sequence<bool>(m, [&](size_t i) { return bits[i]; });
  auto out = pbbslib::pack(v_imap, bits_m);
  pbbslib::free_array(bits);
  size_t out_size = out.size();
  return vertexSubset(n, out_size, out.to_array());
}

template <class data, class F>
inline vertexSubset vertexFilter2(vertexSubsetData<data> V, F filter) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  bool* bits = pbbslib::new_array_no_init<bool>(m);
  V.toSparse();
  par_for(0, m, pbbslib::kSequentialForThreshold, [&](size_t i) {
    auto t = V.vtxAndData(i);
    bits[i] = filter(std::get<0>(t), std::get<1>(t));
  });
  auto v_imap =
      pbbslib::make_sequence<uintE>(m, [&](size_t i) { return V.vtx(i); });
  auto bits_m =
      pbbslib::make_sequence<bool>(m, [&](size_t i) { return bits[i]; });
  auto out = pbbslib::pack(v_imap, bits_m);
  size_t out_size = out.size();
  pbbslib::free_array(bits);
  return vertexSubset(n, out_size, out.to_array());
}

// Adds vertices to a vertexSubset vs.
// Caller must ensure that every v in new_verts is not already in vs
// Note: Mutates the given vertexSubset.
inline void add_to_vsubset(vertexSubset& vs, uintE* new_verts,
                           uintE num_new_verts) {
  if (vs.isDense) {
    par_for(0, num_new_verts, pbbslib::kSequentialForThreshold,
            [&](size_t i) { vs.d[new_verts[i]] = true; });
    vs.m += num_new_verts;
  } else {
    const size_t vs_size = vs.numNonzeros();
    const size_t new_size = num_new_verts + vs_size;
    uintE* all_verts = pbbslib::new_array_no_init<uintE>(new_size);
    par_for(0, new_size, pbbslib::kSequentialForThreshold, [&](size_t i) {
      if (i < vs_size) {
        all_verts[i] = vs.s[i];
      } else {
        all_verts[i] = new_verts[i - vs_size];
      }
    });
    uintE* old_s = vs.s;
    vs.s = all_verts;
    vs.m = new_size;
    if (old_s) {
      pbbslib::free_array(old_s);
    }
  }
}

#ifdef USE_PCM_LIB
inline void print_pcm_stats(SystemCounterState& before_sstate,
                            SystemCounterState& after_sstate, size_t rounds,
                            double elapsed) {
  std::cout << "Instructions per clock:        "
            << (getIPC(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "Total Cycles:                  "
            << (getCycles(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "========= Cache misses/hits ========="
            << "\n";
  std::cout << "L2 Hit ratio:                  "
            << (getL2CacheHitRatio(before_sstate, after_sstate) / rounds)
            << "\n";
  std::cout << "L3 Hit ratio:                  "
            << (getL3CacheHitRatio(before_sstate, after_sstate) / rounds)
            << "\n";
  std::cout << "L2 Misses:                     "
            << (getL2CacheMisses(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "L2 Hits:                       "
            << (getL2CacheHits(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "L3 Misses:                     "
            << (getL3CacheMisses(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "L3 Hits:                       "
            << (getL3CacheHits(before_sstate, after_sstate) / rounds) << "\n";
  std::cout << "========= Bytes read/written ========="
            << "\n";
  auto bytes_read = getBytesReadFromMC(before_sstate, after_sstate) / rounds;
  auto bytes_written =
      getBytesWrittenToMC(before_sstate, after_sstate) / rounds;
  size_t GB = 1024 * 1024 * 1024;
  auto throughput = ((bytes_read + bytes_written) / elapsed) / GB;
  std::cout << "Bytes read:                    " << bytes_read << "\n";
  std::cout << "Bytes written:                 " << bytes_written << "\n";
  std::cout << "Throughput: " << throughput << " GB/s"
            << "\n";
  std::cout << "========= Other statistics ========="
            << "\n";
  std::cout << "Average relative frequency:    "
            << (getActiveRelativeFrequency(before_sstate, after_sstate) /
                rounds)
            << "\n";
}
inline void pcm_init() {
  auto* m = PCM::getInstance();
  if (m->program() != PCM::Success) {
    std::cout << "Could not enable program counters"
              << "\n";
    exit(0);
  }
}
inline size_t get_pcm_state() { return getSystemCounterState(); }
#else
inline void print_pcm_stats(size_t before, size_t after, size_t rounds,
                            double elapsed) {}
inline void pcm_init() {}
inline size_t get_pcm_state() { return (size_t)1; }
#endif

#define run_app(G, APP, rounds)                                      \
  auto before_state = get_pcm_state();                               \
  timer st;                                                          \
  double total_time = 0.0;                                           \
  for (size_t r = 0; r < rounds; r++) {                              \
    total_time += APP(G, P);                                         \
  }                                                                  \
  auto time_per_iter = total_time / rounds;                          \
  std::cout << "time per iter: " << time_per_iter << "\n";           \
  auto after_state = get_pcm_state();                                \
  print_pcm_stats(before_state, after_state, rounds, time_per_iter); \
  G.del();

#define generate_main(APP, mutates)                                        \
  int main(int argc, char* argv[]) {                                       \
    commandLine P(argc, argv, " [-s] <inFile>");                           \
    char* iFile = P.getArgument(0);                                        \
    bool symmetric = P.getOptionValue("-s");                               \
    bool compressed = P.getOptionValue("-c");                              \
    bool weighted = P.getOptionValue("-w");                                \
    assert(weighted == false);                                             \
    bool mmap = P.getOptionValue("-m");                                    \
    bool mmapcopy = mutates || P.getOptionValue("-mc");                    \
    debug(std::cout << "mmapcopy = " << mmapcopy << "\n";);                \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                    \
    if (P.getOptionValue("-b")) { \
      auto G = readUncompressedBinaryGraph<symmetricVertex, pbbslib::empty>(iFile, symmetric, mmap, mmapcopy); \
      alloc_init(G); \
      run_app(G, APP, rounds)                                                \
    } else if (P.getOptionValue("-adj")) { \
      auto G = readUnweightedGraph<symmetricVertex>(iFile, symmetric, mmap); \
      alloc_init(G); \
      run_app(G, APP, rounds)                                                \
    } else { \
      auto G = readCompressedGraph<csv_bytepd_amortized, pbbs::empty>(iFile, symmetric, mmap, mmapcopy); \
      alloc_init(G); \
      run_app(G, APP, rounds)                                                \
    } \
  }

//        auto GA = packed_graph<csv_bytepd_amortized, pbbs::empty>(G);                 \








//      auto G = readCompressedGraph<csv_bytepd_amortized, pbbslib::empty>(       \
//      auto GA = packed_graph<csv_bytepd_amortized, pbbs::empty>(G);           \

//    auto G = readCompressedGraph<csv_bytepd_amortized, intE>(iFile, symmetric, mmap, mmapcopy); \
//   auto G = readWeightedGraph<symmetricVertex>(iFile, symmetric, mmap); \

#define generate_weighted_main(APP, mutates)                             \
  int main(int argc, char* argv[]) {                                     \
    commandLine P(argc, argv, " [-s] <inFile>");                         \
    char* iFile = P.getArgument(0);                                      \
    bool symmetric = P.getOptionValue("-s");                             \
    bool compressed = P.getOptionValue("-c");                            \
    bool weighted = P.getOptionValue("-w");                              \
    assert(weighted == true);                                            \
    bool mmap = P.getOptionValue("-m");                                  \
    bool mmapcopy = mutates;                                             \
    debug(std::cout << "mmapcopy = " << mmapcopy << "\n";);              \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                  \
    if (P.getOptionValue("-b")) { \
      auto G = readUncompressedBinaryGraph<symmetricVertex, intE>(iFile, symmetric, mmap, mmapcopy); \
      alloc_init(G); \
      run_app(G, APP, rounds)                                              \
    } else { \
      auto G = readCompressedGraph<csv_bytepd_amortized, intE>(iFile, symmetric, mmap, mmapcopy); \
      alloc_init(G); \
      run_app(G, APP, rounds)                                              \
    } \
    data_block_allocator::print_stats(); \
  }
