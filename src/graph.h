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

#include <stdlib.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

#include "bridge.h"
#include "compressed_vertex.h"
#include "flags.h"
#include "pbbslib/parallel.h"
#include "vertex.h"

#ifdef NVM
#include <numa.h>
#include <utmpx.h>
#endif

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

template <template <class W> class vertex, class W>
struct graph {
  using w_vertex = vertex<W>;
  using weight_type = W;
  w_vertex* V;

#ifdef NVM
  w_vertex* V0;
  w_vertex* V1;
#endif

  size_t n;
  size_t m;
  bool transposed;
  uintE* flags;
  std::function<void()> deletion_fn;

#ifndef NVM
  w_vertex get_vertex(size_t i) { return V[i]; }
#else
  w_vertex get_vertex(size_t i) {
    //		int cpu = sched_getcpu();
    //		int node = numa_node_of_cpu(cpu);
    if (numanode() == 0) {
      return V0[i];
    } else {
      return V1[i];
    }
  }
#endif

  graph(w_vertex* _V, long _n, long _m, std::function<void()> _d,
        uintE* _flags = NULL)
#ifdef NVM
      : V(_V),
        n(_n),
        m(_m),
        transposed(0),
        flags(_flags),
        deletion_fn(_d),
        V0(_V),
        V1(_V) {
  }
#else
      : V(_V), n(_n), m(_m), transposed(0), flags(_flags), deletion_fn(_d) {
  }
#endif

  void del() {
    if (flags != NULL) pbbslib::free_array(flags);
    deletion_fn();
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    par_for(0, n, 1,
            [&](size_t i) { V[i].mapOutNgh(i, f, parallel_inner_map); });
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
