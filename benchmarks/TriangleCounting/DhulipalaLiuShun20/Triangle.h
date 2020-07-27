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
#include "benchmark.h"
#include "preprocess.h"


namespace gbbs {
using namespace std;


template <class Graph, class F>
inline size_t Triangle(Graph& G, const F& f, commandLine& P) {
  // auto C0 = P.getOptionIntValue("-c", 0);
  using EdgeT = DBTGraph::EdgeT;
  using UpdatesT = pbbs::sequence<pair<EdgeT, bool>>;
  timer t;

  t.start();
  DBTGraph::DyGraph DG = DBTGraph::DyGraph(3, G); t.stop();t.reportTotal("init");

  UpdatesT updates = UTIL::generateEdgeUpdates<EdgeT>(DG.n, 5);
  t.start(); //step 1
  UpdatesT updates_final = Preprocessing(DG, updates); 
  pbbs::sequence<bool> flag = pbbs::sequence<bool>::no_init(updates_final.size());
  par_for(0, updates_final.size(), [&] (size_t i) {
      flag[i] = updates_final[i].second;
  });
  pair<UpdatesT, size_t> tmp = pbbs::split_two(updates_final, flag);
  updates.clear();
  updates = tmp.first;
  size_t m_del = tmp.second;
  size_t m_ins = updates_final.size() - tmp.second;

  
  t.stop();t.reportTotal("preprocess");
  
  UTIL::PrintFunctionItem("1", "m", updates_final.size());


  t.start(); //step 2
  // mark_inserted_edges(updates_final); t.stop();t.reportTotal("preprocess");





  cout << "done" << endl;

  return 0;
}

}  // namespace gbbs
