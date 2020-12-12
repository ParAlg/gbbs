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

#include "gbbs/gbbs.h"
#include "benchmarks/Connectivity/common.h"

namespace gbbs {
namespace simple_union_find {

inline uintE find_compress(uintE i, pbbs::sequence<parent>& parents) {
  uintE pathlen = 1;
  parent j = i;
  if (parents[j] == j) return j;
  do {
    j = parents[j];
    pathlen++;
  } while (parents[j] != j);
  parent tmp;
  while ((tmp=parents[i])>j) {
    parents[i] = j; i=tmp;
  }
  report_pathlen(pathlen);
  return j;
}

void unite_impl(uintE u_orig, uintE v_orig, pbbs::sequence<parent>& parents) {
  parent u = u_orig;
  parent v = v_orig;
  while(1) {
    u = find_compress(u,parents);
    v = find_compress(v,parents);
    if(u == v) break;
    else if (u > v && parents[u] == u && pbbs::atomic_compare_and_swap(&parents[u],u,v)) {
      break;
    }
    else if (v > u && parents[v] == v && pbbs::atomic_compare_and_swap(&parents[v],v,u)) {
      break;
    }
  }
}

struct SimpleUnionAsyncStruct {
  size_t n;
  sequence<parent> parents;
  SimpleUnionAsyncStruct(size_t n) : n(n) {
    parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return (uintE)i; });
  }
  void unite(uintE u, uintE v) {
    unite_impl(u, v, parents);
  }
  sequence<parent> finish() {
    parallel_for(0, n, [&] (size_t i) {
      find_compress(i, parents);
    });
    return std::move(parents);
  }
};

// Outputs a sequence `S` of length `G.n` such that the `i`-th vertex is in
// connected component `S[i]`. The component IDs will be in the range `[0, G.n)`
// but are not necessarily contiguous.
template <class Graph>
inline sequence<parent> SimpleUnionAsync(Graph& G) {
  using W = typename Graph::weight_type;

  size_t n = G.n;
  auto uf = SimpleUnionAsyncStruct(n);

  parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      uf.unite(u, v);
    };
    G.get_vertex(i).mapOutNgh(i, map_f);
  });
  return uf.finish();
}


template <class Seq>
inline size_t num_cc(Seq& labels) {
  size_t n = labels.size();
  auto flags = sequence<uintE>(n + 1, [&](size_t i) { return 0; });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    if (!flags[labels[i]]) {
      flags[labels[i]] = 1;
    }
  });
  pbbslib::scan_add_inplace(flags);
  std::cout << "# n_cc = " << flags[n] << "\n";
  return flags[n];
}

template <class Seq>
inline size_t largest_cc(Seq& labels) {
  size_t n = labels.size();
  // could histogram to do this in parallel.
  auto flags = sequence<uintE>(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    flags[labels[i]] += 1;
  }
  size_t sz = pbbslib::reduce_max(flags);
  std::cout << "# largest_cc has size: " << sz << "\n";
  return sz;
}

}  // namespace simple_union_find
}  // namespace gbbs
