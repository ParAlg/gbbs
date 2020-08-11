#pragma once

#include "gbbs/gbbs.h"
#include "pbbslib/sample_sort.h"
#include "pbbslib/monoid.h"
#include "dynamic_graph.h"

namespace gbbs {
using namespace std;

template <class Graph, class EdgeT>
inline bool dupEdge(const DBTGraph::DyGraph<Graph> &G, const pair<EdgeT, bool> &e){
  if (e.first.first >= G.num_vertices() || e.first.second >= G.num_vertices()){
      cout << "edge out of bound " << endl;
      exit(1);
  }
  return G.haveEdge(e.first) == e.second;
}

// //check duplicated edge, mark delete if exiting and need to delete
// template <class Graph, class EdgeT>
// inline bool dupEdgeDel(DBTGraph::DyGraph<Graph> &G, const pair<EdgeT, bool> &e){
//   if (e.first.first >= G.num_vertices() || e.first.second >= G.num_vertices()){
//       cout << "edge out of bound " << endl;
//       exit(1);
//   }
//   return G.haveEdgeDel(e.first, e.second) == e.second;
// }

template <class EdgeT, class UT>
inline pair<EdgeT, bool> toMyUpdateEdgeT(UT e){
  return make_pair(EdgeT(e.from, e.to), e.weight == 1);
}

template <class Graph, class EdgeT, class UT>
inline pbbs::sequence<pair<EdgeT, bool>> Preprocessing(DBTGraph::DyGraph<Graph> &G, std::vector<UT> &updates){
  size_t n = updates.size();
  // u < v
  parallel_for(0, n, [&](size_t i) {
    uintE u = updates[i].from;
    uintE v = updates[i].to;
    if(u > v){
      updates[i].to = v;
      updates[i].from = u;
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
      return make_pair(updates[i].from, updates[i].to) < make_pair(updates[j].from, updates[j].to);
      }, true);

  pbbs::sequence<size_t> flag = pbbs::sequence<size_t>::no_init(n+1);
  par_for(0, n-1, [&] (size_t i) {
      if(updates[inds[i]].from != updates[inds[i+1]].from || updates[inds[i]].to != updates[inds[i+1]].to){
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
        updates_valid[flag[i]] = toMyUpdateEdgeT<EdgeT, UT>(updates[inds[i]]);
      }
  });
  
  // updates.clear();
  inds.clear();

  // remove inserts/deletes in/notin graph
  pbbs::sequence<pair<EdgeT, bool>> updates_final = pbbs::filter(updates_valid, [&] (const pair<EdgeT, bool>& e) {return !dupEdge(G, e);});
  // pbbs::sequence<pair<EdgeT, bool>> updates_final = pbbs::filter(updates_valid, [&] (const pair<EdgeT, bool>& e) {return !dupEdgeDel(G, e);});

  flag.clear();
  updates_valid.clear();
  return updates_final;
}


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
pair<pbbs::sequence<DBTGraph::VtxUpdate>, pbbs::sequence<size_t>> toCSR(DBTGraph::DyGraph<Graph>& G, pbbs::sequence<pair<EdgeT,bool>> &edgesIn, pbbs::sequence<pair<EdgeT,bool>> &edges, size_t n){
  size_t m = edgesIn.size();
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

  par_for(1, 2*m, [&] (size_t i) {
  if(flag[i-1]!=flag[i]){
    uintE u = getFirst(edges,i);
    vtxNew[flag[i]] = DBTGraph::VtxUpdate(u,i);
    vtxMap[u] = flag[i];
  }});
  uintE u = getFirst(edges,0);
  vtxNew[0] = DBTGraph::VtxUpdate(u,0);
  vtxMap[u] = 0;

  //count D and insert D
  par_for(0, 2*m-1, [&] (size_t i) {
  if(getFirst(edges,i) == getFirst(edges,i+1) && edges[i].second && !edges[i+1].second){
    uintE u = getFirst(edges,i);
    vtxNew[vtxMap[u]].insert_degree = i + 1 - vtxNew[vtxMap[u]].offset;
  }else if(getFirst(edges,i) != getFirst(edges,i+1)){
    uintE u = getFirst(edges,i);
    uintE next_v = getFirst(edges,i+1);
    vtxNew[vtxMap[u]].setDeg(vtxNew[vtxMap[next_v]].offset - vtxNew[vtxMap[u]].offset);
    if(edges[i].second)vtxNew[vtxMap[u]].insert_degree = vtxNew[vtxMap[u]].degree;
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
      size_t s2 = vtxNew[i].insOffset();
      size_t e = vtxNew[i].end();
      vtxNew[i].insert_low_degree = pbbslib::reduce(flag.slice(s,s2 ), monoid);
      vtxNew[i].delete_low_degree = pbbslib::reduce(flag.slice(s2,e ), monoid);
    });

    flag.clear();
    return make_pair(vtxNew,vtxMap);
}

}