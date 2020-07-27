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
#include "benchmark.h"


namespace gbbs {
using namespace std;

template <class Graph, class EdgeT>
inline bool dupEdge(DBTGraph::DyGraph<Graph> G, pair<EdgeT, bool> e){
  if (e.first.first >= G.n || e.first.second >= G.n){
      cout << "edge out of bound " << endl;
      exit(1);
  }
  return G.haveEdge(e.first) == e.second;
}

template <class Graph, class EdgeT>
inline pbbs::sequence<pair<EdgeT, bool>> Preprocessing(DBTGraph::DyGraph<Graph> G, pbbs::sequence<pair<EdgeT, bool>> updates){
  size_t n = updates.size();
  // u < v
  parallel_for(0, n, [&](size_t i) {
    uintE u = updates[i].first.first;
    uintE v = updates[i].first.second;
    if(u > v){
      updates[i].first.first = v;
      updates[i].first.second = u;
    }

  }, 1);

  // nullify, only the chronologically last update
  // first find all updates on the same edge by hashing the updates 
  //by the edge that it is being performed on. 
  //parallel maximum-finding algorithm
  pbbs::sequence<size_t> inds = pbbs::sequence<size_t>::no_init(n);
  par_for(0, n, [&] (size_t i) {
      inds[i] = i;
  });

  pbbs::sample_sort_inplace(inds.slice(),  //check
    [&](const size_t i, const size_t j) {
      return updates[i].first < updates[j].first;
      }, true);

  pbbs::sequence<size_t> flag = pbbs::sequence<size_t>::no_init(n+1);
  par_for(0, n-1, [&] (size_t i) {
      if(updates[inds[i]].first != updates[inds[i+1]].first){
        flag[i] = 1;
      }else{
        flag[i] = 0;
      }
  });
  flag[n-1] = 1;
  flag[n] = 1;
  
  auto monoid = pbbslib::addm<size_t>();
  size_t new_n = pbbs::scan_inplace(flag.slice(), monoid);
  new_n --;

  pbbs::sequence<pair<EdgeT, bool>> updates_valid = pbbs::sequence<pair<EdgeT, bool>>::no_init(new_n);
  par_for(0, n, [&] (size_t i) {
      if(flag[i] != flag[i+1]){
        updates_valid[flag[i]] = updates[inds[i]];
      }
  });
  
  updates.clear();
  inds.clear();

  // remove inserts/deletes in/notin graph, approximate compaction
  pbbs::sequence<pair<EdgeT, bool>> updates_final = pbbs::filter(updates_valid, [&] (pair<EdgeT, bool> e) {return !dupEdge(G, e);});

  flag.clear();
  updates_valid.clear();
  return updates_final;
}



template <class Graph, class F>
inline size_t Triangle(Graph& G, const F& f, commandLine& P) {
  // auto C0 = P.getOptionIntValue("-c", 0);
  timer t;
  t.start();
  using EdgeT = DBTGraph::EdgeT;
  DBTGraph::DyGraph DG = DBTGraph::DyGraph(3, G);
  t.stop();
  t.reportTotal("init");

  // pbbs::sequence<pair<EdgeT, bool>> updates = UTIL::generateEdgeUpdates<EdgeT>(DG.n, 5);
  // pbbs::sequence<pair<EdgeT, bool>> updates_final = Preprocessing(DG, updates);
  // updates->clear();
  // delete updates;
  // UTIL::PrintPairSeq(updates);
  cout << "done" << endl;

  return 0;
}

}  // namespace gbbs
