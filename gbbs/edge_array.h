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

#include "bridge.h"
#include "macros.h"

namespace gbbs {

// Edge Array Representation
template <class W>
struct edge_array {
  using weight_type = W;
  using edge = std::tuple<uintE, uintE, W>;
  edge* E;
  // for sq matrices, num_rows == num_cols
  size_t num_rows;  // n TODO: deprecate #rows/#cols
  size_t num_cols;  // TODO deprecate

  size_t n;
  size_t m;

  // non_zeros is the #edges
  size_t non_zeros;  // m TODO rename to "m"
  void del() { pbbslib::free_array(E); }
  edge_array(edge* _E, size_t r, size_t c, size_t nz)
      : E(_E), num_rows(r), num_cols(c), non_zeros(nz) {
    if (r != c) {
      std::cout << "# edge_array format currently expects square matrix"
                << std::endl;
      exit(0);
    }
    n = r;
    m = nz;
  }
  edge_array() {}
  size_t size() { return non_zeros; }

  pbbs::sequence<edge> to_seq() {
    auto ret = pbbs::sequence<edge>(E, non_zeros);
    non_zeros = 0;
    E = nullptr;
    return std::move(ret);
  }

  template <class F>
  void map_edges(F f, bool parallel_inner_map = true) {
    parallel_for(0, m,
                 [&](size_t i) {
                   uintE u, v;
                   W w;
                   std::tie(u, v, w) = E[i];
                   f(u, v, w);
                 },
                 512);
  }
};

template <class W, class Graph>
inline edge_array<W> to_edge_array(Graph& G) {
  using edge = std::tuple<uintE, uintE, W>;
  size_t n = G.n;
  auto sizes = pbbs::sequence<uintT>(n);
  parallel_for(0, n,
               [&](size_t i) { sizes[i] = G.get_vertex(i).getOutDegree(); });
  size_t m = pbbslib::scan_add_inplace(sizes.slice());
  assert(m == G.m);

  edge* arr = pbbs::new_array_no_init<edge>(m);
  parallel_for(0, n, [&](size_t i) {
    size_t idx = 0;
    uintT offset = sizes[i];
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      arr[offset + idx] = std::make_tuple(u, v, wgh);
      idx++;
    };
    G.get_vertex(i).mapOutNgh(i, map_f, /* parallel = */ false);
  });
  return edge_array<W>(arr, n, n, m);
}

}  // namespace gbbs
