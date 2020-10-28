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

// #define DBT_TOMB_MERGE
#define DBT_USING_TOMB 

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

template <class Graph, class F>
inline tuple<size_t, bool, DSymGraph *> Dynamic_Triangle_Helper(DBTGraph::DyGraph<Graph>* DG, pbbs::sequence<size_t>& vtxMap, const vector<gbbs_io::Edge<int>>& updates, size_t C0, commandLine& P, size_t n, size_t s, size_t e) {
  using EdgeT = DBTGraph::EdgeT;
  using UpdatesT = pbbs::sequence<pair<EdgeT, bool>>;
  using UT = gbbs_io::Edge<int>;
  timer t; t.start();
  size_t m, m_ins;
  size_t delta_triangles_pos, delta_triangles_neg;
  DSymGraph *DGnew;
  size_t new_ct;
    
  if(DG->num_edges() == 0){ // majorRebalancing from empty graph
    //size_t new_ct;
    vector<gbbs_io::Edge<int>> updates_ins;
    for (size_t i = s; i< e; ++i) {
      if(updates[i].weight == 1 && updates[i].from != updates[i].to) updates_ins.emplace_back(updates[i].from, updates[i].to, 1);
    }
    t.next("1. preprocess");
    DBTInternal::PrintFunctionItem("1.", "valid b",  updates_ins.size());
    tie(new_ct, DGnew) = DBTGraph::majorRebalancing(updates_ins, 0, updates_ins.size(), n, DG->get_block_size(), P, true);
    return make_tuple(new_ct, true, DGnew);
  }

  // t.start(); //step 1
  UpdatesT updates_final = DBTInternal::Preprocessing<Graph, EdgeT, UT>(DG, updates, s, e);
  m = updates_final.size();
  t.next("1. preprocess");
  DBTInternal::PrintFunctionItem("1.", "valid b", m);

  // t.start(); //toCSR
  UpdatesT edges = UpdatesT::no_init(2*m);
  pbbs::sequence<DBTGraph::VtxUpdate> vtxNew = DBTInternal::toCSR(DG, vtxMap, updates_final, edges, DG->num_vertices()); // fill vtxMap and edges
  t.next("count degrees");

  auto insertDegrees = pbbs::delayed_sequence<size_t, DBTGraph::VtxUpdateInsDeg>(vtxNew.size(), DBTGraph::VtxUpdateInsDeg(vtxNew));
  auto monoid = pbbslib::addm<size_t>();
  m_ins = pbbs::reduce(insertDegrees, monoid) / 2;
  bool major_rebalanced = false;
  if(DG->majorRebalance(m_ins, m-m_ins)){
    major_rebalanced = true;;
    tie(new_ct, DGnew) = DBTGraph::majorRebalancing(DG, edges,vtxNew, vtxMap, n, P);
    // return make_tuple(new_ct, true, DGnew);
  }else{
  // insertion must be before deletion, because when resizing write OLD_EDGE into tables
  // t.start(); //step 2 mark insert,  some array moves to tables
  // par_for(0, vtxNew.size(), [&] (size_t i) {
  // int nworkers = pbbs::num_workers();
   //pbbs::set_num_workers(1);
   parallel_for(0, vtxNew.size(), [&] (size_t i) {
     DG->markEdgeInsertion(vtxNew[i], edges.slice(vtxNew[i].offset,      vtxNew[i].insOffset()));
     DG->markEdgeDeletion(vtxNew[i],  edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
   }, 1);
  //for(size_t i = 0; i< vtxNew.size(); ++i){
  //  DG->markEdgeInsertion(vtxNew[i], edges.slice(vtxNew[i].offset,      vtxNew[i].insOffset()));
  //  DG->markEdgeDeletion(vtxNew[i],  edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
  //}
  t.next("2. 3. mark insertions + deletions");
  //  pbbs::set_num_workers(nworkers);

  // t.start(); //step 4 and 5 update insertions  and deletions
  // loop over the low degree vertices w, process if the other is high
  // only process each edge once
  par_for(0, vtxNew.size(), [&] (size_t i) {
    DG->updateTable(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].end()));
  });
  t.next("4. 5. update insertions and deletions");

  // t.start(); //step 6. count triangles 
  DBTGraph::TriangleCounts tc = DBTGraph::TriangleCounts();
  // updates final has one copy of each edge
  // for each edge, count the delta triangle caused by counting wedges
  par_for(0, updates_final.size(), [&] (size_t i) {
    EdgeT elocal = updates_final[i].first;
    bool flag = updates_final[i].second;
    DBTGraph::VtxUpdate u = vtxNew[vtxMap[elocal.first]]; // must make a copy here, because countTriangles might swap variables
    DBTGraph::VtxUpdate v = vtxNew[vtxMap[elocal.second]]; // must make a copy here, because countTriangles might swap variables
    DG->countTriangles(u,v,flag, tc);
  });
  pbbs::sequence<size_t> triCounts = tc.report();  //TODO: reuse
  delta_triangles_pos = triCounts[0] + triCounts[1]/2 + triCounts[2]/3;
  delta_triangles_neg = triCounts[3] + triCounts[4]/2 + triCounts[5]/3;
  DBTInternal::PrintFunctionItem("6.", "# tri +", delta_triangles_pos);
  DBTInternal::PrintFunctionItem("6.", "# tri -", delta_triangles_neg);
  tc.clear();
  triCounts.clear();
  t.next("6. count triangles");

  // t.start(); //  first cleanup wedge tables, then re-mark inserts to OLD_EDGE, then remove
#ifdef DBT_TOMB_MERGE
  par_for(0, vtxNew.size(), [&] (size_t i) { //cleanup T and delete if count is 0
    DG->cleanUpTable(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].end()));
  });
#else
  par_for(0, vtxNew.size(), [&] (size_t i) { //cleanup T, called before tables are cleaned up
    DG->cleanUpTable(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].end()), false);
  });
  par_for(0, vtxNew.size(), [&] (size_t i) { //cleanup T, delete 0 wedges
    DG->cleanUpTable(vtxNew[i], edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()), true);
  });
#endif
  par_for(0, vtxNew.size(), [&] (size_t i) { // remark inserts, must be before remove deletes
    DG->cleanUpEdgeInsertion(vtxNew[i], edges.slice(vtxNew[i].offset, vtxNew[i].insOffset()));
  });
  par_for(0, vtxNew.size(), [&] (size_t i) { // remove deletes
    DG->cleanUpEdgeDeletion(vtxNew[i], edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
  });
  t.next("7. clean up tables");
  
  // t.start(); //  minor rebalancing 
  DBTGraph::minorRebalancing(DG, vtxNew, vtxMap);
  t.next("8. 9. update degree + minor rebalancing");
  DG->updateNumEdges(m_ins, m-m_ins);
  new_ct = C0 + delta_triangles_pos - delta_triangles_neg; 
  } //end else (nort major rebalancing)

  
  par_for(0, vtxNew.size(), pbbslib::kSequentialForThreshold, [&] (size_t i) {
    vtxMap[vtxNew[i].id] = EMPTYVMAP;
  });

  updates_final.clear(); // preprocessed edges
  edges.clear(); // sorted preprocessed edges in both direction
  vtxNew.clear(); 
  t.stop(); t.reportTotal("DCountTime");


  // DG->checkStatus();

  return make_tuple(new_ct, major_rebalanced, DGnew);

}

//edges already randomly shuffled
template <class Graph, class F>
inline size_t dynamicBatches(DBTGraph::DyGraph<Graph>* DG, const vector<gbbs_io::Edge<int>>& edges, size_t batch_size, commandLine& P, 
                              size_t n, int batch_offset, size_t C0, bool all_del) {

  size_t num_batch = edges.size()/batch_size;
  std::cout << "Batch Size " << batch_size << std::endl;
  std::cout <<""<< std::endl;
  bool use_new = false; // use new graph (major rebalanced) or old graph
  bool switched = false; 
  size_t count = C0;
  int batch_proc = P.getOptionLongValue("-bp", num_batch+1); 
  DSymGraph *DGold = new DSymGraph();
  DSymGraph *DGnew;
  pbbs::sequence<size_t> vtxMap = pbbs::sequence<size_t>(n, EMPTYVMAP);
  for(int i = batch_offset; i<= batch_offset+ batch_proc; ++i){
    size_t batch_start = i*batch_size;
    if(all_del){
      batch_start = (num_batch-i) * batch_size;
      // j--;
    }
    size_t batch_end = min(batch_size + batch_start, edges.size());
    if(batch_end <= batch_start || batch_start < 0) continue;
    timer t; t.start();
    if(switched){
      tie(count, use_new, DGnew) = Dynamic_Triangle_Helper<DBTGraph::SymGraph, F>(DGold, vtxMap, edges, count, P, n, batch_start, batch_end);
    }else{
      tie(count, use_new, DGnew) = Dynamic_Triangle_Helper<Graph, F>(DG, vtxMap, edges, count, P, n, batch_start, batch_end);
    }
    if(use_new){
      if(switched){DGold->del();}else{DG->del();}
      DGold = DGnew;
      switched = true;
    }

    std::cout << "### Batch " << i << " [" << batch_start << " " << batch_end << "]" << std::endl;
    std::cout << "### Num triangles = " << count << "\n";
    t.stop();t.reportTotal("");
    std::cout << std::endl;
  }
  vtxMap.clear();
  return count;
}

}  //namespace DBTInternal

// if es flag is there assume edges is sorted and there is no duplicates
template <class Graph, class F>
inline size_t Dynamic_Triangle(Graph& G, const vector<gbbs::gbbs_io::Edge<int>>& updates, const F& f, size_t batch_size, commandLine& P) {
  auto C0 = P.getOptionIntValue("-trict", 0);
  bool start_graph = P.getOptionValue("-sg"); 
  bool run_static = P.getOptionValue("-static");
  bool run_static_mix = P.getOptionValue("-staticmix"); 
  int batch_offset = P.getOptionLongValue("-bo", 0); 
  size_t block_size = P.getOptionLongValue("-blocksize", 10000);  
  size_t n = P.getOptionLongValue("-n", 0); 
  int weight = P.getOptionIntValue("-w", 1);    
  timer t;t.start();


  if(start_graph){
    n = G.num_vertices();
    DBTGraph::DyGraph<Graph> * DG = new DBTGraph::DyGraph<Graph>(block_size, G, n);
    t.next("Init");
    return DBTInternal::dynamicBatches<Graph, F>(DG, updates, batch_size, P, n, batch_offset, C0, weight == 2);
  }

  // size_t batch_size = updates.size()/batch_num;
  if(run_static_mix){
    return DBTInternal::staticCountMixed(updates, batch_size, P, n, batch_offset);
  }

  if(run_static){
    return DBTInternal::staticCount(updates, batch_size, P, n, batch_offset);
  }

  if(weight == 0 && P.getOptionValue("-shuffle")){
    C0 = Triangle(G, f, "degree", P); 
    t.next("Static Count");
    DBTGraph::DyGraph<Graph> * DG = new DBTGraph::DyGraph<Graph>(block_size, G, n);
    t.next("Build Dynamic Graph");
    return DBTInternal::dynamicBatches<Graph, F>(DG, updates, batch_size, P, n, batch_offset, C0, false);
  }

  size_t init_graph_end = min(batch_offset * batch_size, updates.size());
  if(weight == 2) init_graph_end =  updates.size();
  DBTGraph::SymGraph G2 = DBTInternal::edge_list_to_symmetric_graph(updates, n, 0, init_graph_end);
  t.next("Build Graph");
  C0 = Triangle(G2, f, "degree", P); 
  t.next("Static Count");
  DBTInternal::DSymGraph * DG = new DBTInternal::DSymGraph(block_size, G2, n);
  G2.del();
  t.next("Build Dynamic Graph");
  return DBTInternal::dynamicBatches<DBTGraph::SymGraph, F>(DG, updates, batch_size, P, n, batch_offset, C0, weight==2);
}

}  // namespace gbbs
