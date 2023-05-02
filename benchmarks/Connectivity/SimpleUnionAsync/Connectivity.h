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

#include "benchmarks/Connectivity/common.h"
#include "gbbs/gbbs.h"

namespace gbbs {
namespace simple_union_find {

template <class IntT>
inline IntT find_compress(uintE i, IntT* parents) {
  IntT j = i;
  if (parents[j] == j) return j;
  do {
    j = parents[j];
  } while (parents[j] != j);
  IntT tmp;
  while ((tmp = parents[i]) > j) {
    parents[i] = j;
    i = tmp;
  }
  return j;
}

template <class IntT>
inline bool unite_impl(uintE u_orig, uintE v_orig, IntT* parents) {
  IntT u = u_orig;
  IntT v = v_orig;
  while (u != v) {
    u = find_compress(u, parents);
    v = find_compress(v, parents);
    if (u > v && parents[u] == u &&
             gbbs::atomic_compare_and_swap(&parents[u], u, v)) {
      return true;
    } else if (v > u && parents[v] == v &&
               gbbs::atomic_compare_and_swap(&parents[v], v, u)) {
      return true;
    }
  }
  return false;
}

template <class IntT>
inline IntT find_compress_atomic(uintE i, IntT* parents) {
  IntT j = i;
  if (gbbs::atomic_load(&parents[j]) == j) return j;
  do {
    j = gbbs::atomic_load(&parents[j]);
  } while (gbbs::atomic_load(&parents[j]) != j);
  IntT tmp;
  while ((tmp = gbbs::atomic_load(&parents[i])) > j) {
    if (!gbbs::atomic_compare_and_swap(&parents[i], tmp, j)) {
      return j;
    }
    i = tmp;
  }
  return j;
}

template <class IntT>
inline bool unite_impl_atomic(uintE u_orig, uintE v_orig,
                              IntT* parents) {
  IntT u = u_orig;
  IntT v = v_orig;
  while (u != v) {
    u = find_compress_atomic(u, parents);
    v = find_compress_atomic(v, parents);
    if (u > v && parents[u] == u &&
             gbbs::atomic_compare_and_swap(&parents[u], u, v)) {
      return true;
    } else if (v > u && parents[v] == v &&
               gbbs::atomic_compare_and_swap(&parents[v], v, u)) {
      return true;
    }
  }
  return false;
}

struct SimpleUnionAsyncStruct {
  size_t n;
  sequence<parent> parents;
  SimpleUnionAsyncStruct(size_t n) : n(n) {
    parents =
        sequence<uintE>::from_function(n, [&](size_t i) { return (uintE)i; });
  }
  void unite(uintE u, uintE v) { unite_impl(u, v, parents.begin()); }
  sequence<parent> finish() {
    parallel_for(0, n, [&](size_t i) { find_compress(i, parents.begin()); });
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

  parallel_for(0, n, [&](size_t i) {
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      uf.unite(u, v);
    };
    G.get_vertex(i).out_neighbors().map(map_f);
  });
  return uf.finish();
}

// Same as SimpleUnionAsync above, but takes a predicate function and
// only unites edges s.t. pred(u, v) = true.
template <class Graph, class Predicate>
inline sequence<parent> CC_predicate(Graph& G, Predicate pred) {
  using W = typename Graph::weight_type;

  size_t n = G.n;
  auto uf = SimpleUnionAsyncStruct(n);

  parallel_for(0, n, [&](size_t i) {
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (pred(u,v,wgh)) {
        uf.unite(u, v);
      }
    };
    G.get_vertex(i).out_neighbors().map(map_f);
  });
  return uf.finish();
}

template <class Seq>
inline size_t num_cc(Seq& labels) {
  size_t n = labels.size();
  auto flags =
      sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    if (!flags[labels[i]]) {
      flags[labels[i]] = 1;
    }
  });
  parlay::scan_inplace(flags);
  std::cout << "# n_cc = " << flags[n] << "\n";
  return flags[n];
}

template <class Seq>
inline size_t largest_cc(Seq& labels) {
  size_t n = labels.size();
  // could histogram to do this in parallel.
  auto flags =
      sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    flags[labels[i]] += 1;
  }
  size_t sz = parlay::reduce_max(flags);
  std::cout << "# largest_cc has size: " << sz << "\n";
  return sz;
}

}  // namespace simple_union_find
}  // namespace gbbs
