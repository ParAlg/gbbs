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

namespace DBTInternal{
using DSymGraph = DBTGraph::DyGraph<DBTGraph::SymGraph>;


// gbbs_io::write_graph_to_file

template <class Graph, class F, class UT>
inline tuple<size_t, bool, DSymGraph *> Dynamic_Triangle_Helper(DBTGraph::DyGraph<Graph>& DG, pbbs::sequence<size_t>& vtxMap, std::vector<UT>& updates, size_t C0, commandLine& P, size_t n, size_t s, size_t e) {
  using EdgeT = DBTGraph::EdgeT;
  using UpdatesT = pbbs::sequence<pair<EdgeT, bool>>;
  timer t;
  size_t m, m_ins;
  size_t delta_triangles_pos, delta_triangles_neg;
  DSymGraph *DGnew;
    
  if(DG.num_edges() == 0){ // mahorRebalancing from empty graph
    size_t new_ct;
    tie(new_ct, DGnew) = DBTGraph::majorRebalancing(updates, s, e, n, DG.get_block_size(), P, true);
    return make_tuple(new_ct, true, DGnew);
  }

  t.start(); //step 1
  UpdatesT updates_final = DBTInternal::Preprocessing<Graph, EdgeT, UT>(DG, updates, s, e);
  m = updates_final.size();
  t.next("1. preprocess");
  DBTInternal::PrintFunctionItem("1.", "valid b", m);

  t.start(); //toCSR
  UpdatesT edges = UpdatesT::no_init(2*m);
  pbbs::sequence<DBTGraph::VtxUpdate> vtxNew = DBTInternal::toCSR(DG, vtxMap, updates_final, edges, DG.num_vertices());
  // pbbs::sequence<DBTGraph::VtxUpdate> vtxNew = result.first;
  // pbbs::sequence<size_t> vtxMap = result.second;
  t.next("count degrees");

  auto insertDegrees = pbbs::delayed_sequence<size_t, DBTGraph::VtxUpdateInsDeg>(vtxNew.size(), DBTGraph::VtxUpdateInsDeg(vtxNew));
  auto monoid = pbbslib::addm<size_t>();
  m_ins = pbbs::reduce(insertDegrees, monoid) / 2;
  if(DG.majorRebalance(m_ins, m-m_ins)){
    size_t new_ct;
    tie(new_ct, DGnew) = DBTGraph::majorRebalancing(DG, edges,vtxNew, vtxMap, n, P);
    return make_tuple(new_ct, true, DGnew);
  }else{
  // insertion must be before deletion, because when resizing write OLD_EDGE into tables
  t.start(); //step 2 mark insert,  some array moves to tables
  par_for(0, vtxNew.size(), [&] (size_t i) {
    DG.markEdgeInsertion(vtxNew[i], edges.slice(vtxNew[i].offset,      vtxNew[i].insOffset()));
    DG.markEdgeDeletion(vtxNew[i],  edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
  });
  t.next("2. 3. mark insertions + deletions");

  t.start(); //step 4 and 5 update insertions  and deletions
  // loop over the low degree vertices w, process if the other is high
  // only process each edge once
  par_for(0, vtxNew.size(), [&] (size_t i) {
    DG.updateTable(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].end()));
  });
  t.next("4. 5. update insertions and deletions");

  t.start(); //step 6. count triangles 
  DBTGraph::TriangleCounts tc = DBTGraph::TriangleCounts();
  // updates final has one copy of each edge
  // for each edge, count the delta triangle caused by counting wedges
  par_for(0, updates_final.size(), [&] (size_t i) {
    EdgeT e = updates_final[i].first;
    bool flag = updates_final[i].second;
    DBTGraph::VtxUpdate u = vtxNew[vtxMap[e.first]]; // must make a copy here, because countTriangles might swap variables
    DBTGraph::VtxUpdate v = vtxNew[vtxMap[e.second]]; // must make a copy here, because countTriangles might swap variables
    DG.countTriangles(u,v,flag, tc);
  });
  pbbs::sequence<size_t> triCounts = tc.report();  //TODO: reuse
  delta_triangles_pos = triCounts[0] + triCounts[1]/2 + triCounts[2]/3;
  delta_triangles_neg = triCounts[3] + triCounts[4]/2 + triCounts[5]/3;
  // DBTInternal::PrintFunctionItem("6.", "# tri +", delta_triangles_pos);
  // DBTInternal::PrintFunctionItem("6.", "# tri -", delta_triangles_neg);
  tc.clear();
  triCounts.clear();
  t.next("6. count triangles");

  t.start(); //  first cleanup wedge tables, then re-mark inserts to OLD_EDGE, then remove
  par_for(0, vtxNew.size(), [&] (size_t i) { //cleanup T, called before tables are cleaned up
    DG.cleanUpTable(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].end()), false);
  });
  par_for(0, vtxNew.size(), [&] (size_t i) { //cleanup T, delete 0 wedges
    DG.cleanUpTable(vtxNew[i], edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()), true); //TODO: change to tombstone
  });
  par_for(0, vtxNew.size(), [&] (size_t i) { // remark inserts, must be before remove deletes
    DG.cleanUpEdgeInsertion(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].insOffset()));
  });
  par_for(0, vtxNew.size(), [&] (size_t i) { // remove deletes
    DG.cleanUpEdgeDeletion(vtxNew[i], edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
  });
  t.next("7. clean up tables");
  
  t.start(); //  minor rebalancing 
  DBTGraph::minorRebalancing(DG, vtxNew, vtxMap);
  t.next("8. 9. update degree + minor rebalancing");
  }

  par_for(0, vtxNew.size(), [&] (size_t i) { //TODO: reuse
    vtxMap[vtxNew[i].id] = EMPTYVMAP;
  });

  updates_final.clear(); // preprocessed edges
  edges.clear(); // sorted preprocessed edges in both direction
  vtxNew.clear(); 
  t.stop(); t.reportTotal("DCountTime");


  // DG.checkStatus();

  return make_tuple(C0 + delta_triangles_pos - delta_triangles_neg, false, DGnew);

}

//edges already randomly shuffled
template <class Graph, class F, class UT>
inline size_t dynamicBatches(DBTGraph::DyGraph<Graph> DG, std::vector<UT>& edges, int num_batch, commandLine& P, size_t n, size_t batch_offset, size_t C0) {
  size_t batch_size = edges.size()/num_batch;
  std::cout << "Batch Size " << batch_size << std::endl;
  bool use_new = false; // use new graph (major rebalanced) or old graph
  bool switched = false; 
  size_t count = C0;
  // tuple<size_t, bool, DSymGraph *> result;
  DSymGraph DGold;
  DSymGraph *DGnew;
  pbbs::sequence<size_t> vtxMap = pbbs::sequence<size_t>(n, EMPTYVMAP);
  for(int i = batch_offset; i< num_batch; ++i){
    size_t batch_end = min(batch_size  * (i+1), edges.size());
    timer t; t.start();
    if(switched){
      tie(count, use_new, DGnew) = Dynamic_Triangle_Helper<DBTGraph::SymGraph, F, UT>(DGold, vtxMap, edges, count, P, n, i*batch_size, batch_end);
    }else{
      tie(count, use_new, DGnew) = Dynamic_Triangle_Helper<Graph, F, UT>(DG, vtxMap, edges, count, P, n, i*batch_size, batch_end);
    }
    if(use_new){
      if(switched){DGold.del();}else{DG.del();}
      DGold = *DGnew;
      switched = true;
    }
    std::cout << "### Batch " << i << std::endl;
    std::cout << "### Num triangles = " << count << "\n";
    t.stop();t.reportTotal("");
    std::cout << std::endl;
  }
  vtxMap.clear();
  return count;
}

}  //namespace DBTInternal

// if es flag is there assume edges is sorted and there is no duplicates
template <class Graph, class F, class UT>
inline size_t Dynamic_Triangle(Graph& G, std::vector<UT>& updates, const F& f, int batch_num, commandLine& P) {
  bool empty_graph = P.getOptionValue("-eg"); 
  size_t block_size = P.getOptionLongValue("-blocksize", 5);        
  timer t;t.start();

  size_t batch_offset = P.getOptionLongValue("-bo", 0);  
  // size_t batch_end = min(P.getOptionLongValue("-be", updates.size()),  updates.size());  

  if(!empty_graph){
    auto C0 = P.getOptionIntValue("-trict", 0);
    DBTGraph::DyGraph<Graph> DG = DBTGraph::DyGraph<Graph>(block_size, G, G.num_vertices());
    t.next("Init");
    return DBTInternal::dynamicBatches<Graph, F, UT>(DG, updates, batch_num, P, G.num_vertices(), batch_offset, C0);
  }

  size_t n = P.getOptionLongValue("-n", 0); 
  // bool sorted_inserts = P.getOptionValue("-alli"); 

  // if(sorted_inserts){ // all edges are inserts updates
  //   DBTInternal::DSymGraph DG = DBTInternal::DSymGraph(block_size, n); 
  //   DBTInternal::DSymGraph DGnew;
  //   t.next("Init");
  //   size_t new_ct = DBTInternal::dynamicBatches(DG, updates, batch_num, P, n, batch_offset);
  //   return new_ct;
  // }

  size_t batch_size = updates.size()/batch_num;

  DBTGraph::SymGraph G2 = DBTInternal::edge_list_to_symmetric_graph(updates, n, 0, batch_offset * batch_size);
  size_t C0 = Triangle(G2, f, "degree", P);  //TODO: which ordering?, how to ini commandline object?
  t.next("Init");
  DBTInternal::DSymGraph DG = DBTInternal::DSymGraph(block_size, G2, n);
  G2.del();
  t.next("Build DG");
  return DBTInternal::dynamicBatches<DBTGraph::SymGraph, F, UT>(DG, updates, batch_num, P, n, batch_offset, C0);
}

}  // namespace gbbs
