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

#include "benchmark.h"
#include "dynamic_graph.h"
#include "gbbs/gbbs.h"
#include "pbbslib/sample_sort.h"
#include "shared.h"
#include <algorithm>
#include <cmath>

namespace gbbs {
using namespace std;

template <class Graph, class F>
inline size_t Triangle(Graph& G, const F& f, commandLine& P) {
  // auto C0 = P.getOptionIntValue("-c", 0);
  using EdgeT = DBTGraph::EdgeT;
  using UpdatesT = pbbs::sequence<pair<EdgeT, bool>>;
  timer t;
  size_t m, m_ins;

  t.start();
  DBTGraph::DyGraph DG = DBTGraph::DyGraph(5, G); t.stop();t.reportTotal("init");

  UpdatesT updates = UTIL::generateEdgeUpdates<EdgeT>(DG.num_vertices(), 10);

  t.start(); //step 1
  UpdatesT updates_final = Preprocessing(DG, updates); // Already sorted
  m = updates_final.size();
  updates.clear();
  t.stop();t.reportTotal("1. preprocess");
  UTIL::PrintFunctionItem("1.", "m", m);

  t.start(); //toCSR
  UpdatesT edges = UpdatesT::no_init(2*m);
  auto result = toCSR(DG, updates_final, edges, DG.num_vertices());
  pbbs::sequence<DBTGraph::VtxUpdate> vtxNew = result.first;
  pbbs::sequence<size_t> vtxMap = result.second;
  t.stop();t.reportTotal("process update graph and create structure");

  /////Perform merging of sorted adjacency lists////
  t.start();
  int n = updates_final.size();
  pbbs::sequence<size_t> unique_flags =
      pbbs::sequence<size_t>::no_init(n); // count the set of unique adj lists
  par_for(0, n - 1, [&](size_t i) {
    if (updates[i].first.first != updates[i + 1].first.first) {
      unique_flags[i] = 1;
    } else {
      unique_flags[i] = 0;
    }
  });
  unique_flags[n - 1] = 1;

  auto monoid = pbbslib::addm<size_t>();
  size_t num_distinct_adj_lists =
      pbbs::scan_inplace(unique_flags.slice(), monoid);

  pbbs::sequence<size_t> adj_list_starts =
      pbbs::sequence<size_t>::no_init(num_distinct_adj_lists);
  adj_list_starts[0] = 0;
  par_for(0, num_distinct_adj_lists - 1, [&](size_t i) {
    if (unique_flags[i] != unique_flags[i + 1]) {
      adj_list_starts[unique_flags[i + 1]] = i + 1;
    }
  });

  // Clearing flags
  unique_flags.clear();

  // Merging with DG adjacency lists

  t.end(); t.reportTotal("Merging sorted adjacency lists.");

  cout << "done" << endl;

  return 0;
}

}  // namespace gbbs
