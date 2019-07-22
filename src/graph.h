// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <string>
#include <stdlib.h>
#include <fstream>
#include <functional>
#include <iostream>

#include "bridge.h"
#include "compressed_vertex.h"
#include "flags.h"
#include "vertex.h"
#include "pbbslib/parallel.h"

#ifdef NVM
#include <utmpx.h>
#include <numa.h>
#endif

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

template <class vertex>
struct graph {
  vertex* V;

#ifdef NVM
  vertex* V0;
  vertex* V1;
#endif

  size_t n;
  size_t m;
  bool transposed;
  uintE* flags;
  std::function<void()> deletion_fn;
  std::function<graph<vertex>()> copy_fn;

#ifndef NVM
  vertex get_vertex(size_t i) {
      	return V[i];
  }
#else
  vertex get_vertex(size_t i) {
//		int cpu = sched_getcpu();
//		int node = numa_node_of_cpu(cpu);
    if (numanode() == 0) {
      return V0[i];
    } else {
      return V1[i];
    }
  }
#endif

  graph(vertex* _V, long _n, long _m, std::function<void()> _d,
        uintE* _flags = NULL)
#ifdef NVM
      : V(_V), n(_n), m(_m), transposed(0), flags(_flags), deletion_fn(_d), V0(_V), V1(_V) {}
#else
      : V(_V), n(_n), m(_m), transposed(0), flags(_flags), deletion_fn(_d) {}
#endif

  graph(vertex* _V, long _n, long _m, std::function<void()> _d,
        std::function<graph<vertex>()> _c, uintE* _flags = NULL)
#ifdef NVM
      : V(_V), n(_n), m(_m), transposed(0), flags(_flags), deletion_fn(_d), copy_fn(_c), V0(_V), V1(_V) {}
#else
      : V(_V), n(_n), m(_m), transposed(0), flags(_flags), deletion_fn(_d), copy_fn(_c) {}
#endif

  auto copy() -> graph<vertex> { return copy_fn(); }

  void del() {
    if (flags != NULL) pbbslib::free_array(flags);
    deletion_fn();
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map=true) {
    par_for(0, n, 1, [&] (size_t i) {
      V[i].mapOutNgh(i, f, parallel_inner_map);
    });
  }
};

inline auto get_deletion_fn(void* V, void* edges) -> std::function<void()> {
  auto df = [&](void* V, void* edges) {
    pbbslib::free_array(V);
    pbbslib::free_array(edges);
  };
  return std::bind(df, V, edges);
}

inline auto get_deletion_fn(void* V, void* in_edges, void* out_edges)
    -> std::function<void()> {
  auto df = [&](void* V, void* in_edges, void* out_edges) {
    pbbslib::free_array(V);
    pbbslib::free_array(in_edges);
    pbbslib::free_array(out_edges);
  };
  return std::bind(df, V, in_edges, out_edges);
}

template <class vertex, class E>
inline std::function<graph<vertex>()> get_copy_fn(vertex* V, E* edges, size_t n,
                                                  size_t m, size_t sizeofe) {
  auto df = [&](vertex* V, E* edges, size_t n, size_t m, size_t sizeofe) {
    auto NV = pbbslib::new_array_no_init<vertex>(n);
    auto NE = pbbslib::new_array_no_init<E>(sizeofe);
    par_for(0, sizeofe, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { NE[i] = edges[i]; });
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      NV[i].setOutDegree(V[i].getOutDegree());
      size_t off = (V[i].getOutNeighbors() - edges);
      NV[i].setOutNeighbors(NE + off);
    });
    graph<vertex> G = graph<vertex>(NV, n, m, get_deletion_fn(NV, NE));
    G.copy_fn = get_copy_fn(NV, NE, n, m, sizeofe);
    //    std::cout << "Returning copied graph" << "\n";
    return G;
  };
  return std::bind(df, V, edges, n, m, sizeofe);
}

template <class vertex, class E>
inline std::function<graph<vertex>()> get_copy_fn(vertex* V, E* in_edges,
                                                  E* out_edges, size_t n,
                                                  size_t m, size_t m_in,
                                                  size_t m_out) {
  auto df = [&](vertex* V, E* in_edges, E* out_edges, size_t n, size_t m,
                size_t m_in, size_t m_out) {
    auto NV = pbbslib::new_array_no_init<vertex>(n);
    auto Nin = pbbslib::new_array_no_init<E>(m_in);
    auto Nout = pbbslib::new_array_no_init<E>(m_out);
    par_for(0, m_in, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { Nin[i] = in_edges[i]; });
    par_for(0, m_in, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { Nout[i] = out_edges[i]; });
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      NV[i].setOutDegree(V[i].getOutDegree());
      NV[i].setInDegree(V[i].getInDegree());
      size_t out_off = (V[i].getOutNeighbors() - out_edges);
      NV[i].setOutNeighbors(Nout + out_off);
      size_t in_off = (V[i].getInNeighbors() - in_edges);
      NV[i].setOutNeighbors(Nin + in_off);
    });
    graph<vertex> G = graph<vertex>(
        V, n, m, get_deletion_fn((void*)NV, (void*)Nin, (void*)Nout));
    G.copy_fn = get_copy_fn(NV, Nin, Nout, n, m, m_in, m_out);
    return G;
  };
  return std::bind(df, V, in_edges, out_edges, n, m, m_in, m_out);
}
