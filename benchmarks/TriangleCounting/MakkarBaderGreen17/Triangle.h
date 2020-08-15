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
#include "sample_sort.h"
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
  t.stop();t.reportTotal("count degrees");

  /////Perform merging of sorted adjacency lists////
  pbbs::sequence<size_t> unique_flags =
      pbbs::sequence<size_t>::no_init(n); // count the set of unique adj lists
  par_for(0, n - 1, [&](size_t i) {
    if (updates[i].first.first != updates[i + 1].first.first) {
      unique_flags[i] = 1;
    } else {
      s unique_flags[i] = 0;
    }
  });
  flag[n - 1] = 1;

  auto monoid = pbbslib::addm<size_t>();
  size_t num_distinct_adj_lists =
      pbbs::scan_inplace(unique_flags.slice(), monoid);

  pbbs::sequence<size_t> adj_list_starts =
      pbbs::sequence<size_t>::no_init(num_distinct_adj_lists);
  adj_list_starts[0] = 0;
  par_for(0, num_distinct_adj_list - 1, [&](size_t i) {
    if (unique_flags[i] != unique_flags[i + 1]) {
      adj_list_starts[unique_flags[i + 1]] = i + 1;
    }
  });

  // Clearing flags
  unique_flags.clear();

  // Merging with DG adjacency lists

  ////////////////Remainder of the file is from yushangdi's branch////////////
  t.start(); //step 6. count triangles // updates final has one copy of each edge
  DBTGraph::TriangleCounts tc = DBTGraph::TriangleCounts();
  par_for(0, updates_final.size(), [&] (size_t i) {
    EdgeT e = updates_final[i].first;
    bool flag = updates_final[i].second;
    DBTGraph::VtxUpdate u = vtxNew[vtxMap[e.first]];
    DBTGraph::VtxUpdate v = vtxNew[vtxMap[e.second]];
    DG.countTriangles(u,v,flag, tc);
  });
  pbbs::sequence<size_t> triCounts = tc.report();
  size_t delta_triangles_pos = triCounts[0] + triCounts[1]/2 + triCounts[2]/3;
  size_t delta_triangles_neg = triCounts[3] + triCounts[4]/2 + triCounts[5]/3;
  UTIL::PrintFunctionItem("6.", "# tri +", delta_triangles_pos);
  UTIL::PrintFunctionItem("6.", "# tri -", delta_triangles_neg);
  tc.clear();
  triCounts.clear();
  t.stop();t.reportTotal("6. count triangles");

  t.start(); //  first remark inserts, then remove
  par_for(0, vtxNew.size(), [&] (size_t i) { // remark inserts
    DG.cleanUpEdgeInsertion(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].insOffset()));
  });
  par_for(0, vtxNew.size(), [&] (size_t i) { // remove deletes
    DG.cleanUpEdgeDeletion(vtxNew[i], edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
  });
  par_for(0, vtxNew.size(), [&] (size_t i) { //cleanup T
    DG.cleanUpTable(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].end()));
  });
  t.stop();t.reportTotal("7. clean up tables");


  auto insertDegrees = pbbs::delayed_sequence<size_t, DBTGraph::VtxUpdateInsDeg>(vtxNew.size(), DBTGraph::VtxUpdateInsDeg(vtxNew));
  auto monoid = pbbslib::addm<size_t>();
  m_ins = pbbs::reduce(insertDegrees, monoid) / 2;
  if(DG.majorRebalance(2 * m_ins - m)){
    cout <<  "major rebalancing not implemented " << endl;
  }else{
  // t.start(); //  minor rebalancing
  // move between tables
  // move table to array

  // move if not in array and change from L->H or H->L
  // change top level first, then bottom level

  // remove empty table

  // resize T since num high changes

  // t.stop();t.reportTotal("8. 9. minor rebalancing");
  }

  // insertDegrees.clear();
  updates_final.clear();


  // resize table (increase)
  // find table using oldD, increase

  // minor rebalancing

  // // resize table (decrease)
  // t.stop();t.reportTotal("resizing");


  //CLEANUP TODO: move to array,  remove empty tables
  //TODO: do not keep block nodes in T,  in minor rebalance add to T

  cout << "done" << endl;

  return 0;
}

}  // namespace gbbs
