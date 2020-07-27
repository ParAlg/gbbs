#pragma once

#include "gbbs/gbbs.h"
#include "pbbslib/sample_sort.h"
#include "pbbslib/monoid.h"
#include "dynamic_graph.h"

namespace gbbs {
using namespace std;

template <class Graph, class EdgeT>
inline bool dupEdge(const DBTGraph::DyGraph<Graph> &G, pair<EdgeT, bool> &e){
  if (e.first.first >= G.n || e.first.second >= G.n){
      cout << "edge out of bound " << endl;
      exit(1);
  }
  return G.haveEdge(e.first) == e.second;
}

template <class Graph, class EdgeT>
inline pbbs::sequence<pair<EdgeT, bool>> Preprocessing(const DBTGraph::DyGraph<Graph> &G, pbbs::sequence<pair<EdgeT, bool>> &updates){
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


}