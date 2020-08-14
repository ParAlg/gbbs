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

// #define TURNOFFMAJOR

#include "shared.h"
#include "dynamic_graph.h"
#include "benchmark.h"
#include "preprocess.h"
#include "rebalancing.h"



namespace gbbs {
using namespace std;

// gbbs_io::write_graph_to_file

template <class Graph, class F, class UT>
inline size_t Dynamic_Triangle_Helper(DBTGraph::DyGraph<Graph>& DG, std::vector<UT>& updates, size_t C0, commandLine& P) {
  using EdgeT = DBTGraph::EdgeT;
  using UpdatesT = pbbs::sequence<pair<EdgeT, bool>>;
  timer t;
  size_t m, m_ins;
  size_t delta_triangles_pos, delta_triangles_neg;

  t.start(); //step 1
  UpdatesT updates_final = DBTInternal::Preprocessing<Graph, EdgeT, UT>(DG, updates);
  m = updates_final.size();
  updates.clear();
  t.stop();t.reportTotal("1. preprocess");
  DBTInternal::PrintFunctionItem("1.", "valid b", m);

  t.start(); //toCSR
  UpdatesT edges = UpdatesT::no_init(2*m);
  auto result = DBTInternal::toCSR(DG, updates_final, edges, DG.num_vertices());
  pbbs::sequence<DBTGraph::VtxUpdate> vtxNew = result.first;
  pbbs::sequence<size_t> vtxMap = result.second;
  t.stop();t.reportTotal("count degrees");

  auto insertDegrees = pbbs::delayed_sequence<size_t, DBTGraph::VtxUpdateInsDeg>(vtxNew.size(), DBTGraph::VtxUpdateInsDeg(vtxNew));
  auto monoid = pbbslib::addm<size_t>();
  m_ins = pbbs::reduce(insertDegrees, monoid) / 2;
  if(DG.majorRebalance(m_ins, m-m_ins)){
    // size_t new_m = DG.num_edges() + 2*m_ins - m;
    DBTGraph::DyGraph<DBTGraph::SymGraph> DGnew;
    size_t new_ct = DBTGraph::majorRebalancing(DG, DGnew, edges,vtxNew, vtxMap, P);
    // DG.del();
    return new_ct;
  }else{
  // insertion must be before deletion, because when resizing write OLD_EDGE into tables
  t.start(); //step 2 mark insert,  some array moves to tables
  par_for(0, vtxNew.size(), [&] (size_t i) {
    DG.markEdgeInsertion(vtxNew[i], edges.slice(vtxNew[i].offset,      vtxNew[i].insOffset()));
    DG.markEdgeDeletion(vtxNew[i],  edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
  });
  t.stop();t.reportTotal("2. 3. mark insertions + deletions");

  t.start(); //step 4 and 5 update insertions  and deletions
  // loop over the low degree vertices w, process if the other is high
  // only process each edge once
  par_for(0, vtxNew.size(), [&] (size_t i) {
    DG.updateTable(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].end()));
  });
  t.stop();t.reportTotal("4. 5. update insertions and deletions");

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
  delta_triangles_pos = triCounts[0] + triCounts[1]/2 + triCounts[2]/3;
  delta_triangles_neg = triCounts[3] + triCounts[4]/2 + triCounts[5]/3;
  DBTInternal::PrintFunctionItem("6.", "# tri +", delta_triangles_pos);
  DBTInternal::PrintFunctionItem("6.", "# tri -", delta_triangles_neg);
  tc.clear();
  triCounts.clear();
  t.stop();t.reportTotal("6. count triangles");

  t.start(); //  first remark inserts, then remove
  par_for(0, vtxNew.size(), [&] (size_t i) { //cleanup T
    DG.cleanUpTable(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].end()), false);
  });
  par_for(0, vtxNew.size(), [&] (size_t i) { //cleanup T, delete 0 wedges
    DG.cleanUpTable(vtxNew[i], edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()), true); //TODO: change to tombstone
  });
  par_for(0, vtxNew.size(), [&] (size_t i) { // remark inserts
    DG.cleanUpEdgeInsertion(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].insOffset()));
  });
  par_for(0, vtxNew.size(), [&] (size_t i) { // remove deletes
    DG.cleanUpEdgeDeletion(vtxNew[i], edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
  });
  t.stop();t.reportTotal("7. clean up tables");
  
  t.start(); //  minor rebalancing 
  DBTGraph::minorRebalancing(DG, vtxNew, vtxMap);
  t.stop();t.reportTotal("8. 9. update degree + minor rebalancing");
  }

  par_for(0, vtxNew.size(), [&] (size_t i) {
    vtxMap[vtxNew[i].id] = EMPTYVMAP;
  });

  //TODO: do not keep block nodes in T,  in minor rebalance add to T
  updates_final.clear();
  edges.clear();
  vtxNew.clear(); vtxMap.clear();

  return C0 + delta_triangles_pos - delta_triangles_neg;

}

// if es flag is there assume edges is sorted and there is no duplicates
template <class Graph, class F, class UT>
inline size_t Dynamic_Triangle(Graph& G, std::vector<UT>& updates, const F& f, commandLine& P) {
  bool empty_graph = P.getOptionValue("-eg"); 
  size_t block_size = P.getOptionLongValue("-blocksize", 5);        
   

  // DBTInternal::generateEdgeUpdates<EdgeT>(DG.num_vertices(), 10);

  if(!empty_graph){
    auto C0 = P.getOptionIntValue("-trict", 0);
    timer t;t.start();
    DBTGraph::DyGraph<Graph> DG = DBTGraph::DyGraph<Graph>(block_size, G);
    t.stop();t.reportTotal("init");
    return Dynamic_Triangle_Helper<Graph, F, UT>(DG, updates, C0, P); 
  }

  size_t n = P.getOptionLongValue("-n", 0);  
  size_t batch_offset = P.getOptionLongValue("-bo", 0);  
  bool edges_sorted = P.getOptionValue("-es"); 

  if(batch_offset == 0 && edges_sorted){ // all edges are inserts updates
    DBTGraph::DyGraph<DBTGraph::SymGraph> DG = DBTGraph::DyGraph<DBTGraph::SymGraph>(block_size, n); 
    DBTGraph::DyGraph<DBTGraph::SymGraph> DGnew;
    size_t  new_ct = DBTGraph::majorRebalancing(updates, 0, updates.size(), DG.get_block_size(), DGnew, P);
    return new_ct;
  }

  vector<gbbs::gbbs_io::Edge<pbbs::empty>> edges = DBTGraph::getEdgeVec(updates, 0, batch_offset);
  vector<UT> updates2 = DBTGraph::getEdgeVecWeighted(updates, batch_offset, updates.size());
  DBTGraph::SymGraph G2 = gbbs::gbbs_io::edge_list_to_symmetric_graph(edges);
  auto C0 = Triangle(G, f, "degree", P);  //TODO: which ordering?, how to ini commandline object?
  DBTGraph::DyGraph<DBTGraph::SymGraph> DG = DBTGraph::DyGraph<DBTGraph::SymGraph>(block_size, G2);
  G.del();
  return Dynamic_Triangle_Helper<DBTGraph::SymGraph, F, UT>(DG, updates2, C0, P); 

}

}  // namespace gbbs
