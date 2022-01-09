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

#include "MinimumSpanningForest_lib.h"

#include "benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
//#include "benchmarks/MinimumSpanningForest/Kruskal/MinimumSpanningForest.h"

#include "gbbs/gbbs.h"

namespace gbbs {
namespace compiled {

sequence<edge> MinimumSpanningForest(symmetric_unweighted_graph& G) {
  return gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
}
sequence<uint_edge> MinimumSpanningForest(symmetric_uint32_graph& G) {
  return gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
}
sequence<float_edge> MinimumSpanningForest(symmetric_float_graph& G) {
  return gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
}
sequence<double_edge> MinimumSpanningForest(symmetric_double_graph& G) {
  return gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
}

sequence<edge> MinimumSpanningForest(symmetric_unweighted_cgraph& G) {
  return gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
}
sequence<uint_edge> MinimumSpanningForest(symmetric_uint32_cgraph& G) {
  return gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
}
sequence<float_edge> MinimumSpanningForest(symmetric_float_cgraph& G) {
  return gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
}
sequence<double_edge> MinimumSpanningForest(symmetric_double_cgraph& G) {
  return gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
}

}  // namespace compiled
}  // namespace gbbs
