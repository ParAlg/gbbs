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

/*
#pragma once

#include "benchmarks/Connectivity/common.h"
#include "gbbs/gbbs.h"

namespace gbbs {
namespace simple_union_find {

inline bool is_parent(uintE i, sequence<parent>& parents){
  auto max_val = parents[i];
  std::size_t max_bit = sizeof(parent) * 8;
  parent check_bit = (max_val >> (max_bit - 1)) & 1U;
  return check_bit != 0 || max_val == i;
}

inline bool is_fake_parent(uintE i, sequence<parent>& parents){
  auto max_val = parents[i];
  std::size_t max_bit = sizeof(parent) * 8;
  parent check_bit = (max_val >> (max_bit - 1)) & 1U;
  return check_bit != 0;
}

inline bool is_parent_atomic (uintE i, sequence<parent>& parents){
  auto max_val = gbbs::atomic_load(&parents[i]);
  std::size_t max_bit = sizeof(parent) * 8;
  parent check_bit = (max_val >> (max_bit - 1)) & 1U;
  return check_bit != 0 || max_val == i;
}

inline void set_fake_parent(uintE i, sequence<parent>& parents, parent fake_parent){
  // Set up fake_parent
  std::size_t max_bit = sizeof(parent) * 8;
  parent one = 1;
  fake_parent ^= (one << (max_bit - 1));

  // We'll set parent to fake_parent only if parent is the root and has not
  // been previously set to a fake_parent; that is to say, parents[i] == i
  gbbs::atomic_compare_and_swap(&parents[i], i, fake_parent);
}

inline parent extract_fake_parent(uintE i, sequence<parent>& parents) {
  auto max_val = parents[i];
  std::size_t max_bit = sizeof(parent) * 8;
  parent one = 1;
  parent check_bit = (max_val >> (max_bit - 1)) & 1U;
  if (check_bit != 0) {
    max_val ^= (one << (max_bit - 1));
    return static_cast<parent>(max_val);
  }
  return max_val;
}

inline parent pretend_parent(uintE i, sequence<parent>& parents) {
  if (is_parent(i, parents)) return i;
  return parents[i];
}

inline bool replace_parent(uintE i, sequence<parent>& parents, parent replace){
  // As long as parents[i] is either i or a fake parent, we must try to replace
  while(true) {
    auto max_val = parents[i];
    std::size_t max_bit = sizeof(parent) * 8;
    parent check_bit = (max_val >> (max_bit - 1)) & 1U;

    if (check_bit == 0 && max_val != i) return false;

    if (gbbs::atomic_compare_and_swap(&parents[i], max_val, replace)) return true;
  }
  return false;
}

inline void print_uf(sequence<parent>& parents) {
  for (size_t i = 0; i < parents.size(); i++) {
    if (!is_parent(i, parents) || parents[i] == i) std::cout << i << ": " << parents[i] << std::endl;
    else std::cout << i << " fake: " << extract_fake_parent(i, parents) << std::endl;
  }
}

inline uintE find_compress(uintE i, sequence<parent>& parents) {
  uintE pathlen = 1;
  parent j = i;
  if (is_parent(j, parents)) return j;
  do {
    j = pretend_parent(j, parents);
    pathlen++;
  } while (!is_parent(j, parents));
  parent tmp;
  while ((tmp = pretend_parent(i, parents)) > j) {
    parents[i] = parents[j];
    i = tmp;
  }
  report_pathlen(pathlen);
  return j;
}

inline void unite_impl(uintE u_orig, uintE v_orig, sequence<parent>& parents) {
  parent u = u_orig;
  parent v = v_orig;
  while (1) {
    u = find_compress(u, parents);
    v = find_compress(v, parents);
    if (u == v)
      break;
    else if (u > v && is_parent(u, parents) && replace_parent(u, parents, v)) {
      break;
    } else if (v > u && is_parent(v, parents) && replace_parent(v, parents, u)) {
      break;
    }
  }
}

inline uintE find_compress_atomic(uintE i, sequence<parent>& parents) {
  uintE pathlen = 1;
  parent j = i;
  if (is_parent_atomic(j, parents)) return j;
  do {
    j = gbbs::atomic_load(&parents[j]);
    pathlen++;
  } while (!is_parent_atomic(j, parents));
  parent tmp;
  while ((tmp = gbbs::atomic_load(&parents[i])) > j) {
    if (!gbbs::atomic_compare_and_swap(&parents[i], tmp, j)) {
      return j;
    }
    i = tmp;
  }
  report_pathlen(pathlen);
  return j;
}

inline void unite_impl_atomic(uintE u_orig, uintE v_orig,
                              sequence<parent>& parents) {
  parent u = u_orig;
  parent v = v_orig;
  while (1) {
    u = find_compress_atomic(u, parents);
    v = find_compress_atomic(v, parents);
    if (u == v)
      break;
    else if (u > v && is_parent(u, parents) && replace_parent(u, parents, v)) {
      break;
    } else if (v > u && is_parent(v, parents) && replace_parent(v, parents, u)) {
      break;
    }
  }
}

struct SimpleUnionAsyncStruct {
  size_t n;
  sequence<parent> parents;
  SimpleUnionAsyncStruct(size_t n) : n(n) {
    parents =
        sequence<uintE>::from_function(n, [&](size_t i) { return (uintE)i; });
  }
  void unite(uintE u, uintE v) { unite_impl(u, v, parents); }
  sequence<parent> finish() {
    parallel_for(0, n, [&](size_t i) { find_compress(i, parents); });
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
*/


#pragma once

#include "benchmarks/Connectivity/common.h"
#include "gbbs/gbbs.h"

namespace gbbs {
namespace simple_union_find {

inline uintE find_compress(uintE i, sequence<parent>& parents) {
  uintE pathlen = 1;
  parent j = i;
  if (parents[j] == j) return j;
  do {
    j = parents[j];
    pathlen++;
  } while (parents[j] != j);
  parent tmp;
  while ((tmp = parents[i]) > j) {
    parents[i] = j;
    i = tmp;
  }
  report_pathlen(pathlen);
  return j;
}

inline void unite_impl(uintE u_orig, uintE v_orig, sequence<parent>& parents) {
  parent u = u_orig;
  parent v = v_orig;
  while (1) {
    u = find_compress(u, parents);
    v = find_compress(v, parents);
    if (u == v)
      break;
    else if (u > v && parents[u] == u &&
             gbbs::atomic_compare_and_swap(&parents[u], u, v)) {
      break;
    } else if (v > u && parents[v] == v &&
               gbbs::atomic_compare_and_swap(&parents[v], v, u)) {
      break;
    }
  }
}

inline uintE find_compress_atomic(uintE i, sequence<parent>& parents) {
  uintE pathlen = 1;
  parent j = i;
  if (gbbs::atomic_load(&parents[j]) == j) return j;
  do {
    j = gbbs::atomic_load(&parents[j]);
    pathlen++;
  } while (gbbs::atomic_load(&parents[j]) != j);
  parent tmp;
  while ((tmp = gbbs::atomic_load(&parents[i])) > j) {
    if (!gbbs::atomic_compare_and_swap(&parents[i], tmp, j)) {
      return j;
    }
    i = tmp;
  }
  report_pathlen(pathlen);
  return j;
}

inline void unite_impl_atomic(uintE u_orig, uintE v_orig,
                              sequence<parent>& parents) {
  parent u = u_orig;
  parent v = v_orig;
  while (1) {
    u = find_compress_atomic(u, parents);
    v = find_compress_atomic(v, parents);
    if (u == v)
      break;
    else if (u > v && parents[u] == u &&
             gbbs::atomic_compare_and_swap(&parents[u], u, v)) {
      break;
    } else if (v > u && parents[v] == v &&
               gbbs::atomic_compare_and_swap(&parents[v], v, u)) {
      break;
    }
  }
}

struct SimpleUnionAsyncStruct {
  size_t n;
  sequence<parent> parents;
  SimpleUnionAsyncStruct(size_t n) : n(n) {
    parents =
        sequence<uintE>::from_function(n, [&](size_t i) { return (uintE)i; });
  }
  void unite(uintE u, uintE v) { unite_impl(u, v, parents); }
  sequence<parent> finish() {
    parallel_for(0, n, [&](size_t i) { find_compress(i, parents); });
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