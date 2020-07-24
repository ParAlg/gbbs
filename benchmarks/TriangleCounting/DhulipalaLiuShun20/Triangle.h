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

#include <algorithm>
#include <cmath>
#include "pbbslib/sample_sort.h"
#include "pbbslib/monoid.h"
#include "gbbs/gbbs.h"
// #include "gbbs/pbbslib/sparse_table.h"
// #include "two_level_tables.h"
#include "dynamic_graph.h"


namespace gbbs {
using namespace std;




// // struct hash_pair {
// //   inline size_t operator () (const std::tuple<uintE, uintE>& t) {
// //     size_t l = std::min(std::get<0>(t), std::get<1>(t));
// //     size_t r = std::max(std::get<0>(t), std::get<1>(t));
// //     size_t key = (l << 32) + r;
// //     return pbbslib::hash64_2(key);
// //   }
// // };

template <class Graph, class F>
inline size_t Triangle(Graph& G, const F& f, commandLine& P) {
  // auto C0 = P.getOptionIntValue("-c", 0);
  // test(G);
  DBTGraph::DyGraph DG = DBTGraph::DyGraph(100, G);

  // auto state = Initialize<Graph>(G);

  return 0;
}


// template <class Graph>
// inline auto Initialize(Graph& G){
//   using K = uintE;
//   using V = int;
//   using BT = std::tuple<K, V>;
//   // using W = typename Graph::weight_type;
//   size_t n = G.n;
//   BT empty = std::make_tuple(UINT_E_MAX, -1);

//   auto HH = NestHash::nested_table<K, V, hash_vertex>(n, empty, hash_vertex());
//   auto HL = NestHash::nested_table<K, V, hash_vertex>(n, empty, hash_vertex());
//   auto LH = NestHash::nested_table<K, V, hash_vertex>(n, empty, hash_vertex());
//   auto LL = NestHash::nested_table<K, V, hash_vertex>(n, empty, hash_vertex());
//   auto T  = NestHash::nested_table<K, size_t, hash_vertex>(2*n, empty, hash_vertex());
//   // auto D  = sequence<size_t>(n);

//   using E = decltype(HH);
//   using F = decltype(T);
  
//   auto state = BPDTriangleCountState<Graph, E, F>(G, HH, HL, LH, LL, T);

//   return state;
// }

}  // namespace gbbs
