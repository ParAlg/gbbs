#pragma once

#include "gbbs/gbbs.h"
#include "pbbslib/sample_sort.h"
#include "pbbslib/monoid.h"
#include "dynamic_graph.h"

namespace gbbs {
using namespace std;


namespace DBTInternal{

template <class weight_type>
symmetric_graph<symmetric_vertex, weight_type>
edge_list_to_symmetric_graph(const std::vector<gbbs_io::Edge<weight_type>>& edge_list, size_t num_vertices) {
  using edge_type = typename symmetric_vertex<weight_type>::edge_type;

  if (edge_list.empty()) {
    return symmetric_graph<symmetric_vertex, weight_type>{};
  }

  pbbs::sequence<gbbs_io::Edge<weight_type>> edges_both_directions(2 * edge_list.size());
  par_for(0, edge_list.size(), [&](const size_t i) {
      const gbbs_io::Edge<weight_type>& edge = edge_list[i];
      edges_both_directions[2 * i] = edge;
      edges_both_directions[2 * i + 1] =
        gbbs_io::Edge<weight_type>{edge.to, edge.from, edge.weight};
  });
  constexpr auto compare_endpoints = [](
      const gbbs_io::Edge<weight_type>& left,
      const gbbs_io::Edge<weight_type>& right) {
    return std::tie(left.from, left.to) < std::tie(right.from, right.to);
  };
  pbbs::sequence<gbbs_io::Edge<weight_type>> t_edges =
    pbbs::remove_duplicates_ordered(edges_both_directions, compare_endpoints);
  pbbs::sequence<gbbs_io::Edge<weight_type>> edges = pbbs::filter(t_edges, [&] (const gbbs_io::Edge<weight_type>& e) {return e.from  != e.to;});
  t_edges.clear();
  const size_t num_edges = edges.size();
  // const size_t num_vertices = internal::get_num_vertices_from_edges(edges);
  pbbs::sequence<vertex_data> vertex_data =
    gbbs_io::internal::sorted_edges_to_vertex_data_array(num_vertices, edges);

  edge_type* edges_array = pbbs::new_array_no_init<edge_type>(num_edges);
  par_for(0, num_edges, [&](const size_t i) {
    const gbbs_io::Edge<weight_type>& edge = edges[i];
    edges_array[i] = std::make_tuple(edge.to, edge.weight);
  });

  auto vertex_data_array = vertex_data.to_array();
  return symmetric_graph<symmetric_vertex, weight_type>{
    vertex_data_array,
    num_vertices,
    num_edges,
    [=] () { pbbslib::free_arrays(vertex_data_array, edges_array); },
    edges_array};
}

// return true if an insert edge is in graph or a delete edge is not in graph
template <class Graph, class EdgeT>
inline bool dupEdge(const DBTGraph::DyGraph<Graph> &G, const pair<EdgeT, bool> &e){
  // if (e.first.first >= G.num_vertices() || e.first.second >= G.num_vertices()){
  //     cout << "edge out of bound " << endl;
  //     abort();
  // }
  if(e.first.first == e.first.second){return true;}
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

// change graphio edge type into our edge type 
// if (from > to), swap 
template <class EdgeT, class UT>
inline pair<EdgeT, bool> toMyUpdateEdgeT(UT e){
  bool rev = e.to < e.from;
  if (rev) return make_pair(EdgeT(e.to, e.from), e.weight == 1);
  return make_pair(EdgeT(e.from, e.to), e.weight == 1);
}

// given raw updates, give valid updates
// remove duplicates, leave only last update that's in/not in graph
template <class Graph, class EdgeT, class UT>
inline pbbs::sequence<pair<EdgeT, bool>> Preprocessing(DBTGraph::DyGraph<Graph> &G, std::vector<UT> &t_updates){
  size_t n = t_updates.size();
  
  // change to our type
  pbbs::sequence<pair<EdgeT, bool>> updates(n);
  par_for(0, n, [&](const size_t i) {
      updates[i] = toMyUpdateEdgeT<EdgeT, UT>(t_updates[i]);
  });
  t_updates.clear();

  // nullify, leave only the chronologically last update
  // sort indices instead of edges directly
  pbbs::sequence<size_t> inds = pbbs::sequence<size_t>::no_init(n);
  par_for(0, n, [&] (size_t i) {inds[i] = i;});

  pbbs::sample_sort_inplace(inds.slice(),  //check
    [&](const size_t i, const size_t j) {
      return updates[i].first < updates[j].first;
      }, true);

  pbbs::sequence<size_t> flag = pbbs::sequence<size_t>::no_init(n+1); // flag[i] == 1 if i is the last update of edge inds[i]
  par_for(0, n-1, [&] (const size_t i) {
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

  // remove inserts/deletes in/notin graph
  pbbs::sequence<pair<EdgeT, bool>> updates_final = pbbs::filter(updates_valid, [&] (const pair<EdgeT, bool>& e) {return !dupEdge(G, e);});
  // pbbs::sequence<pair<EdgeT, bool>> updates_final = pbbs::filter(updates_valid, [&] (const pair<EdgeT, bool>& e) {return !dupEdgeDel(G, e);});

  flag.clear();
  updates_valid.clear();
  return updates_final;
}

// edges is edges in both directions
// sort edges by (edges, bool) and put info into vtxNew and vtxMap
// true is before false
// vtxNew is filled with offset, degree, and insert degree
template <class EdgeT, class VTX>
inline void computeOffsets(pbbs::sequence<pair<EdgeT,bool>> &edges, pbbs::sequence<VTX> &vtxNew, pbbs::sequence<size_t> &vtxMap, pbbs::sequence<size_t> flag = NULL){
  pbbs::sample_sort_inplace(edges.slice(), [&](const pair<EdgeT,bool>& i, const pair<EdgeT,bool>& j) {
    if(i.first.first == j.first.first) return i.second && !j.second;
      return i.first.first < j.first.first; 
    });

  size_t edgeL = edges.size();
  bool clearflag = false;
  if(flag == NULL){
    flag = pbbs::sequence<size_t>::no_init(edgeL+1);
    clearflag = true;
  }
  //find offsets of vertices
  par_for(0, edgeL-1, [&] (size_t i) {
    if(DBTGraph::getFirst(edges,i) != DBTGraph::getFirst(edges,i+1)){flag[i] = 1;
    }else{flag[i] = 0;}});
  flag[edgeL-1] = 1;
  flag[edgeL] = 1;
  auto monoid = pbbslib::addm<size_t>();
  size_t numVtx = pbbs::scan_inplace(flag.slice(), monoid) - 1 ;
  vtxNew =  pbbs::sequence<VTX>::no_init(numVtx);

  // compute offsets
  par_for(1, edgeL, [&] (size_t i) {
  if(flag[i-1]!=flag[i]){
    uintE u = DBTGraph::getFirst(edges,i);
    vtxNew[flag[i]] = VTX(u,i);
    vtxMap[u] = flag[i];
  }});
  uintE u = DBTGraph::getFirst(edges,0);
  vtxNew[0] = VTX(u,0);
  vtxMap[u] = 0;

  //count D and insert D
  par_for(0, edgeL-1, [&] (size_t i) {
  if(DBTGraph::getFirst(edges,i) == DBTGraph::getFirst(edges,i+1) && edges[i].second && !edges[i+1].second){
    uintE u = DBTGraph::getFirst(edges,i);
    vtxNew[vtxMap[u]].setInsDeg(i + 1 - vtxNew[vtxMap[u]].offset);
  }else if(DBTGraph::getFirst(edges,i) != DBTGraph::getFirst(edges,i+1)){
    uintE u = DBTGraph::getFirst(edges,i);
    uintE next_v = DBTGraph::getFirst(edges,i+1);
    vtxNew[vtxMap[u]].setDeg(vtxNew[vtxMap[next_v]].offset - vtxNew[vtxMap[u]].offset);
    if(edges[i].second)vtxNew[vtxMap[u]].setInsDeg(vtxNew[vtxMap[u]].degree);
  }
  });
  vtxNew[numVtx-1].setDeg(edgeL - vtxNew[numVtx-1].offset);
  if(edges[edgeL-1].second) vtxNew[numVtx-1].setInsDeg(vtxNew[numVtx-1].degree);

  if(clearflag) flag.clear();
}


//TODO: not keeping vtxMap if not used later
template <class Graph, class EdgeT>
pair<pbbs::sequence<DBTGraph::VtxUpdate>, pbbs::sequence<size_t>> toCSR(DBTGraph::DyGraph<Graph>& G, pbbs::sequence<pair<EdgeT,bool>> &edgesIn, pbbs::sequence<pair<EdgeT,bool>> &edges, size_t n){
  size_t m = edgesIn.size();
  pbbs::sequence<DBTGraph::VtxUpdate> vtxNew;
  pbbs::sequence<size_t> vtxMap = pbbs::sequence<size_t>(n, EMPTYVMAP);
  // pbbs::sequence<pair<EdgeT,bool>> edges = pbbs::sequence<pair<EdgeT,bool>>::no_init(2*m);
  pbbs::sequence<size_t> flag = pbbs::sequence<size_t>::no_init(2*m+1);
  auto monoid = pbbslib::addm<size_t>();


  //double edges
  par_for(0, m, [&] (size_t i) {
    edges[2*i] = edgesIn[i];
    edges[2*i+1] = make_pair(EdgeT(DBTGraph::getSecond(edgesIn,i), DBTGraph::getFirst(edgesIn,i)), edgesIn[i].second);
  });

  computeOffsets<EdgeT, DBTGraph::VtxUpdate>(edges, vtxNew, vtxMap, flag);

  //count lowD
    par_for(0, 2*m, [&] (size_t i) {
      flag[i] = G.is_low_v(DBTGraph::getSecond(edges,i));
    });
    par_for(0, vtxNew.size(), [&] (size_t i) {
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
}