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
#include "gbbs/gbbs.h"

#include "dynamic_graph.h"

namespace gbbs {

template <class Graph, class F>
inline size_t Dynamic_Triangle(
    Graph& G,
    const vector<gbbs::gbbs_io::Edge<int>>& updates,
    const F& f,
    int batch_num,
    commandLine& P) {
  size_t n = P.getOptionLongValue("-n", 1000000);
  size_t batch_size = 1000;
  size_t num_batches = (updates.size() + batch_size - 1) / batch_size;

  // Just convert to sequence for convenince.
  using Edge = gbbs::gbbs_io::Edge<int>;
  pbbs::sequence<Edge> U(updates.size());
  parallel_for(0, updates.size(), [&] (size_t i) {
    U[i].from = updates[i].from;
    U[i].to = updates[i].to;
    U[i].weight = updates[i].weight;
  });

  timer t;
  t.start();

  auto DG = gbbs::DynamicGraph(n);

  for (size_t i=0; i<num_batches; i++) {
    // process batch i
    size_t batch_start = i*batch_size;
    size_t batch_end = std::min(updates.size(), (i+1)*batch_size);
    auto batch = U.slice(batch_start, batch_end);
    timer bt; bt.start();
    DG.process_batch(batch);
    bt.stop(); bt.reportTotal("batch time");
//    size_t tc = DG.triangle_count();
//    std::cout << "Triangle count = " << tc << std::endl;
    std::cout << std::endl;
  }
  DG.report_stats();

//  t.end();
//  t.reportTotal("total processing time");
}

}  // namespace gbbs
