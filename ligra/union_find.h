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

#include <tuple>

#include "speculative_for.h"
#include "macros.h"

struct UnionFind {
  size_t n;
  intT* parents;
  UnionFind(size_t _n) : n(_n) {
    parents = pbbslib::new_array_no_init<intT>(n);
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { parents[i] = -1; });
  }

  intT find(int32_t i) {
    if (parents[i] < 0) return i;
    intT j = parents[i];
    if (parents[j] < 0) return j;
    do
      j = parents[j];
    while (parents[j] >= 0);
    intT tmp;
    while ((tmp = parents[i]) != j) {
      parents[i] = j;
      i = tmp;
    }
    return j;
  }

  void link(intT u, intT v) { parents[u] = v; }

  void clear() { pbbslib::free_array(parents); }
};

// edges: <uintE, uintE, W>
template <class intT, class Edges, class ST, class UF>
struct UnionFindStep {
  using res = reservation<intT>;
  using storage = std::tuple<intT, intT>;
  Edges& E;
  res* R;
  ST& inST;
  UF& uf;
  storage* indices;

  size_t n;

  void clear() { pbbslib::free_array(indices); }

  UnionFindStep(Edges& _E, res* _R, ST& ist, UF& _uf)
      : E(_E), R(_R), inST(ist), uf(_uf) {
    n = uf.n;
    indices = pbbslib::new_array_no_init<storage>(E.non_zeros);
  }

  bool reserve(intT i) {
    assert(i < E.non_zeros);
    auto e = E.E[i];
    intT u = uf.find(std::get<0>(e));
    intT v = uf.find(std::get<1>(e));
    if (u != v) {
      indices[i] = std::make_tuple(u, v);
      R[v].reserve(i);
      R[u].reserve(i);
      assert(u < n);
      assert(v < n);
      return 1;  // active
    } else
      return 0;  // done
  }

  bool commit(intT i) {
    assert(i < E.non_zeros);
    // read back u and v from 'reserve'
    auto st = indices[i];
    intT u = std::get<0>(st), v = std::get<1>(st);
    if (u >= n || v >= n) {
      std::cout << "u = " << u << " v = " << v << " i = " << i << "\n";
      exit(0);
    }
    //    assert(u < n); assert(v < n);
    if (R[v].checkReset(i)) {
      R[u].checkReset(i);
      uf.link(v, u);
      inST[i] = 1;
      return 1;
    } else if (R[u].checkReset(i)) {
      uf.link(u, v);
      inST[i] = 1;
      return 1;
    } else
      return 0;
  }
};

template <class intT, class Edges, class R, class ST, class UF>
inline UnionFindStep<intT, Edges, ST, UF> make_uf_step(Edges& e, R r, ST& ist,
                                                       UF& uf) {
  return UnionFindStep<intT, Edges, ST, UF>(e, r, ist, uf);
}
