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
#include "shared.h"


namespace gbbs {
using namespace std;

template <class EdgeT>
inline uintE getFirst(pbbs::sequence<pair<EdgeT,bool>> edges, size_t i){
  return edges[i].first.first;
}

template <class EdgeT>
inline uintE getSecond(pbbs::sequence<pair<EdgeT,bool>> edges, size_t i){
  return edges[i].first.second;
}

//TODO: not keeping vtxMap if not used later
template <class Graph, class EdgeT>
pair<pbbs::sequence<DBTGraph::VtxUpdate>, pbbs::sequence<size_t>> toCSR(DBTGraph::DyGraph<Graph>& G, pbbs::sequence<pair<EdgeT,bool>> &edgesIn, pbbs::sequence<pair<EdgeT,bool>> edges, size_t n){
  size_t m = edges.size();
  pbbs::sequence<DBTGraph::VtxUpdate> vtxNew;
  pbbs::sequence<size_t> vtxMap = pbbs::sequence<size_t>::no_init(n);
  // pbbs::sequence<pair<EdgeT,bool>> edges = pbbs::sequence<pair<EdgeT,bool>>::no_init(2*m);
  pbbs::sequence<size_t> flag = pbbs::sequence<size_t>::no_init(2*m+1);

  //sort edges
  par_for(0, m, [&] (size_t i) {
    edges[2*i] = edgesIn[i];
    edges[2*i+1] = make_pair(EdgeT(getSecond(edgesIn,i), getFirst(edgesIn,i)), edgesIn[i].second);
  });
  // size_t bits = pbbslib::log2_up(n);
  // pbbslib::integer_sort_inplace(A.slice(), get_u, bits); // which sort?
  pbbs::sample_sort_inplace(edges.slice(), [&](const pair<EdgeT,bool>& i, const pair<EdgeT,bool>& j) {
    if(i.first.first == j.first.first) return i.second && !j.second;
      return i.first.first < j.first.first; 
    });

  //find offsets of vertices
  par_for(0, 2*m-1, [&] (size_t i) {
    if(getFirst(edges,i) != getFirst(edges,i+1)){flag[i] = 1;
    }else{flag[i] = 0;}});
  flag[2*m-1] = 1;
  flag[2*m] = 1;
  auto monoid = pbbslib::addm<size_t>();
  size_t numVtx = pbbs::scan_inplace(flag.slice(), monoid) - 1 ;
  vtxNew =  pbbs::sequence<DBTGraph::VtxUpdate>::no_init(numVtx);

  par_for(0, 2*m, [&] (size_t i) {
  if(flag[i]!=flag[i+1]){
    uintE u = getFirst(edges,i);
    vtxNew[flag[i]] = DBTGraph::VtxUpdate(u,i);
    vtxMap[u] = flag[i];
  }});

  //count D and insert D
  par_for(0, 2*m-1, [&] (size_t i) {
  if(getFirst(edges,i) == getFirst(edges,i+1) && edges[i].second && !edges[i+1].second){
    uintE u = getFirst(edges,i);
    vtxNew[vtxMap[u]].insert_degree = i-vtxNew[vtxMap[u]].offset;
  }else if(getFirst(edges,i) != getFirst(edges,i+1)){
    uintE u = getFirst(edges,i);
    uintE next_v = getFirst(edges,i+1);
    vtxNew[vtxMap[u]].setDeg(vtxNew[vtxMap[next_v]].offset - vtxNew[vtxMap[u]].offset);
    if(edges[i].second)vtxNew[vtxMap[u]].insert_degree = i-vtxNew[vtxMap[u]].offset;
  }
  });
  vtxNew[numVtx-1].setDeg(2*m - vtxNew[numVtx-1].offset);
  if(edges[2*m-1].second) vtxNew[numVtx-1].insert_degree = vtxNew[numVtx-1].degree;

  //count lowD
    par_for(0, 2*m, [&] (size_t i) {
      flag[i] = G.is_low_v(getSecond(edges,i));
    });
    par_for(0, numVtx, [&] (size_t i) {
      size_t s = vtxNew[i].offset;
      size_t s2 = vtxNew[i].offset + vtxNew[i].insert_degree;
      // size_t e = s+vtxNew[i].degree;
      vtxNew[i].insert_low_degree = pbbslib::reduce(flag.slice(s,s2 ), monoid);
      // size_t insert_high_degree = pbbslib::reduce(flag.slice(s2,e), monoid);
      // cilk_sync;
      // vtxNew[i].insert_degree = vtxNew[i].insert_low_degree + insert_high_degree;
    });

    flag.clear();
    return make_pair(vtxNew,vtxMap);
}


template <class Graph, class F>
inline size_t Triangle(Graph& G, const F& f, commandLine& P) {
  // auto C0 = P.getOptionIntValue("-c", 0);
  using EdgeT = DBTGraph::EdgeT;
  using UpdatesT = pbbs::sequence<pair<EdgeT, bool>>;
  timer t;
  size_t m, m_ins;

  t.start();
  DBTGraph::DyGraph DG = DBTGraph::DyGraph(3, G); t.stop();t.reportTotal("init");

  // UpdatesT updates = UTIL::generateEdgeUpdates<EdgeT>(DG.num_vertices(), 10);
  UpdatesT updates = UpdatesT::no_init(8);
  updates[0] = make_pair(EdgeT(1,2), true);
  updates[1] = make_pair(EdgeT(3,2), true);
  updates[2] = make_pair(EdgeT(1,2), true);
  updates[3] = make_pair(EdgeT(1,2), false);
  updates[4] = make_pair(EdgeT(1,2), true);
  updates[5] = make_pair(EdgeT(4,2), false);//remove
  updates[6] = make_pair(EdgeT(0,4), false);
  updates[7] = make_pair(EdgeT(1,4), true);//add

  t.start(); //step 1
  UpdatesT updates_final = Preprocessing(DG, updates);  // mark delete edge as well
  m = updates_final.size();
  updates.clear();
  t.stop();t.reportTotal("1. preprocess + 3. mark deletes");
  UTIL::PrintFunctionItem("1.", "m", m);

  t.start(); //toCSR
  UpdatesT edges = UpdatesT::no_init(2*m);
  auto result = toCSR(DG, updates_final, edges, DG.num_vertices());
  pbbs::sequence<DBTGraph::VtxUpdate> vtxNew = result.first;
  pbbs::sequence<size_t> vtxMap = result.second;
  t.stop();t.reportTotal("count degrees");

  t.start(); //step 2 mark insert ,  some array moves to tables
  par_for(0, vtxNew.size(), [&] (size_t i) {
    DG.markEdgeInsertion(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].offset + vtxNew[i].insert_degree));
  });
  t.stop();t.reportTotal("2. mark insertions");

  t.start(); //step 3 mark deletion 
  par_for(0, vtxNew.size(), [&] (size_t i) {
    DG.markEdgeDeletion(vtxNew[i], edges.slice(vtxNew[i].offset + vtxNew[i].insert_degree, vtxNew[i].offset + vtxNew[i].degree));
  });
  t.stop();t.reportTotal("3. mark deletions");

  t.start(); //step 4 and 5 update insertions  and deletions
  // loop over the low degree vertices, process if the other is high
  // only process each edge once
  par_for(0, vtxNew.size(), [&] (size_t i) {
    DBTGraph::VtxUpdate w = vtxNew[i];
    DG.updateTable(w, edges.slice(vtxNew[i].offset, vtxNew[i].offset + vtxNew[i].degree));
  });
  t.stop();t.reportTotal("4.5. update insertions and deletions");

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
  size_t delta_triangles = triCounts[0] + triCounts[1]/2 + triCounts[2]/3 - triCounts[3] - triCounts[4]/2 - triCounts[5]/3;
  UTIL::PrintFunctionItem("6.", "tri", delta_triangles);
  tc.clear();
  triCounts.clear();


  t.stop();t.reportTotal("6. count triangles");

  t.start(); //  first remark inserts, then remove
  par_for(0, vtxNew.size(), [&] (size_t i) { // remark inserts
    DG.cleanUpEdgeInsertion(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].offset + vtxNew[i].insert_degree));
  });
  par_for(0, vtxNew.size(), [&] (size_t i) { // remove deletes
    DG.cleanUpEdgeDeletion(vtxNew[i], edges.slice(vtxNew[i].offset + vtxNew[i].insert_degree, vtxNew[i].offset + vtxNew[i].degree));
  });
  par_for(0, vtxNew.size(), [&] (size_t i) {
    DBTGraph::VtxUpdate w = vtxNew[i];
    DG.cleanUpTable(w, edges.slice(vtxNew[i].offset, vtxNew[i].offset + vtxNew[i].degree));
  });
  t.stop();t.reportTotal("7. clean up tables");
  

  auto insertDegrees = pbbs::delayed_sequence<size_t, DBTGraph::VtxUpdateInsDeg>(vtxNew.size(), DBTGraph::VtxUpdateInsDeg(vtxNew));
  auto monoid = pbbslib::addm<size_t>();
  m_ins = pbbs::reduce(insertDegrees, monoid);
  if(DG.majorRebalance(2 * m_ins - m)){
    cout <<  "major rebalancing not implemented " << endl;
  }else{
  // t.start(); //  minor rebalancing 
  // move between tables 
  // move table to array

  // move if not in array and change from L->H or H->L
  // change top level first, then bottom level

  // remove empty table
  
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
