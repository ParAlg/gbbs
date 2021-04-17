#pragma once

#include "gbbs/gbbs.h"
#include "pbbslib/sample_sort.h"
#include "pbbslib/monoid.h"
#include "dynamic_graph.h"

namespace gbbs {
using namespace std;


namespace DBTInternal{

template <class weight_type>
symmetric_graph<symmetric_vertex, gbbs::empty>
edge_list_to_symmetric_graph(const std::vector<gbbs_io::Edge<weight_type>>& edge_list, size_t num_vertices, size_t s, size_t e) {
  using edge_type = typename symmetric_vertex<gbbs::empty>::edge_type;

  size_t edgelistsize = e-s;
  if (edge_list.empty() || edgelistsize == 0) {
    return symmetric_graph<symmetric_vertex, gbbs::empty>{};
  }

  sequence<gbbs_io::Edge<gbbs::empty>> edges_both_directions(2 * edgelistsize);
  par_for(0, edgelistsize, kDefaultGranularity, [&](const size_t i) {
      const gbbs_io::Edge<weight_type>& edge = edge_list[s+i];
      edges_both_directions[2 * i] = gbbs_io::Edge<gbbs::empty>{edge.from, edge.to, gbbs::empty()};
      edges_both_directions[2 * i + 1] =
        gbbs_io::Edge<gbbs::empty>{edge.to, edge.from, gbbs::empty()};
  });
  constexpr auto compare_endpoints = [](
      const gbbs_io::Edge<gbbs::empty>& left,
      const gbbs_io::Edge<gbbs::empty>& right) {
    return std::tie(left.from, left.to) < std::tie(right.from, right.to);
  };
  sequence<gbbs_io::Edge<gbbs::empty>> t_edges = pbbslib::remove_duplicates_ordered(edges_both_directions, compare_endpoints);
  sequence<gbbs_io::Edge<gbbs::empty>> edges = pbbslib::filter(t_edges, [&] (const gbbs_io::Edge<gbbs::empty>& ee) {return ee.from  != ee.to;});
  t_edges.clear();
  const size_t num_edges = edges.size();
  // const size_t num_vertices = internal::get_num_vertices_from_edges(edges);
  vertex_data* data =
    gbbs_io::internal::sorted_edges_to_vertex_data_array(num_vertices, edges);

  edge_type* edges_array = pbbslib::new_array_no_init<edge_type>(num_edges);
  par_for(0, num_edges, kDefaultGranularity, [&](const size_t i) {
    const gbbs_io::Edge<gbbs::empty>& edge = edges[i];
    edges_array[i] = std::make_tuple(edge.to, edge.weight);
  });
  edges.clear();
  return symmetric_graph<symmetric_vertex, gbbs::empty>{
    data,
    num_vertices,
    num_edges,
    [=] () { pbbslib::free_arrays(data, edges_array); },
    edges_array};
}

// return true if an insert edge is in graph or a delete edge is not in graph
template <class Graph, class EdgeT>
inline bool dupEdge(const DBTGraph::DyGraph<Graph> *G, const pair<EdgeT, bool> &e){
  // if (e.first.first >= G->num_vertices() || e.first.second >= G->num_vertices()){
  //     cout << "edge out of bound " << endl;
  //     abort();
  // }
  if(e.first.first == e.first.second){return true;}
  return G->haveEdge(e.first) == e.second;
}

// //check duplicated edge, mark delete if exiting and need to delete
// template <class Graph, class EdgeT>
// inline bool dupEdgeDel(DBTGraph::DyGraph<Graph> &G, const pair<EdgeT, bool> &e){
//   if (e.first.first >= G->num_vertices() || e.first.second >= G->num_vertices()){
//       cout << "edge out of bound " << endl;
//       exit(1);
//   }
//   return G->haveEdgeDel(e.first, e.second) == e.second;
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
inline sequence<pair<EdgeT, bool>> Preprocessing(DBTGraph::DyGraph<Graph> *G, const std::vector<UT> &t_updates, size_t s, size_t e){
  size_t n = e-s;//t_updates.size();

  // change to our type
  sequence<pair<EdgeT, bool>> updates(n);
  par_for(0, n, kDefaultGranularity, [&](const size_t i) {
      updates[i] = toMyUpdateEdgeT<EdgeT, UT>(t_updates[s+i]);
  });

  // nullify, leave only the chronologically last update
  // sort indices instead of edges directly
  sequence<size_t> inds = sequence<size_t>::uninitialized(n);
  par_for(0, n, kDefaultGranularity, [&] (size_t i) {inds[i] = i;});

  pbbslib::sample_sort_inplace(make_slice(inds),  //check
    [&](const size_t i, const size_t j) {
      return updates[i].first < updates[j].first;
      }, true);

  sequence<size_t> flag = sequence<size_t>::uninitialized(n+1); // flag[i] == 1 if i is the last update of edge inds[i]
  par_for(0, n-1, kDefaultGranularity, [&] (const size_t i) {
      if(updates[inds[i]].first != updates[inds[i+1]].first){
        flag[i] = 1;
      }else{
        flag[i] = 0;
      }
  });
  flag[n-1] = 1;
  flag[n] = 1;

  auto monoid = pbbslib::addm<size_t>();
  size_t new_n = pbbslib::scan_inplace(make_slice(flag), monoid);
  new_n --;

  sequence<pair<EdgeT, bool>> updates_valid = sequence<pair<EdgeT, bool>>::uninitialized(new_n);
  par_for(0, n, kDefaultGranularity, [&] (size_t i) {
      if(flag[i] != flag[i+1]){
        updates_valid[flag[i]] = updates[inds[i]];
      }
  });

  updates.clear();
  inds.clear();

  // remove inserts/deletes in/notin graph
  sequence<pair<EdgeT, bool>> updates_final = pbbslib::filter(updates_valid, [&] (const pair<EdgeT, bool>& eee) {return !dupEdge(G, eee);});
  // sequence<pair<EdgeT, bool>> updates_final = pbbslib::filter(updates_valid, [&] (const pair<EdgeT, bool>& e) {return !dupEdgeDel(G, e);});

  flag.clear();
  updates_valid.clear();
  return updates_final;
}

// edges is edges in both directions
// sort edges by (edges, bool) and put info into vtxNew and vtxMap
// true is before false
// vtxNew is filled with offset, degree, and insert degree
template <class EdgeT, class VTX>
inline sequence<VTX> computeOffsets(pbbslib::range<pair<EdgeT,bool> *> edges, pbbslib::range<size_t *> vtxMap, sequence<size_t> &flag ){
  pbbslib::sample_sort_inplace(edges, [&](const pair<EdgeT,bool>& i, const pair<EdgeT,bool>& j) {
    if(i.first.first == j.first.first) return i.second && !j.second;
      return i.first.first < j.first.first;
    });

  size_t edgeL = edges.size();
  bool clearflag = false;
  if(flag.empty()){
    flag = sequence<size_t>::uninitialized(edgeL+1);
    clearflag = true;
  }
  //find offsets of vertices
  par_for(0, edgeL-1, kDefaultGranularity, [&] (size_t i) {
    if(DBTGraph::getFirst(edges,i) != DBTGraph::getFirst(edges,i+1)){flag[i] = 1;
    }else{flag[i] = 0;}});
  flag[edgeL-1] = 1;
  flag[edgeL] = 1;
  auto monoid = pbbslib::addm<size_t>();
  size_t numVtx = pbbslib::scan_inplace(make_slice(flag), monoid) - 1 ;
  sequence<VTX> vtxNew =  sequence<VTX>::uninitialized(numVtx);

  // compute offsets
  par_for(1, edgeL, kDefaultGranularity, [&] (size_t i) {
  if(flag[i-1]!=flag[i]){
    uintE u = DBTGraph::getFirst(edges,i);
    vtxNew[flag[i]] = VTX(u,i);
    vtxMap[u] = flag[i];
  }});
  uintE uu = DBTGraph::getFirst(edges,0);
  vtxNew[0] = VTX(uu,0);
  vtxMap[uu] = 0;

  //count D and insert D
  par_for(0, edgeL-1, kDefaultGranularity, [&] (size_t i) {
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
  return vtxNew;
}


//TODO: not keeping vtxMap if not used later
template <class Graph, class EdgeT>
sequence<DBTGraph::VtxUpdate> toCSR(DBTGraph::DyGraph<Graph>* G, sequence<size_t>& vtxMap, sequence<pair<EdgeT,bool>> &edgesIn, sequence<pair<EdgeT,bool>> &edges, size_t n){
  size_t m = edgesIn.size();
  sequence<DBTGraph::VtxUpdate> vtxNew;
  // sequence<pair<EdgeT,bool>> edges = sequence<pair<EdgeT,bool>>::uninitialized(2*m);
  sequence<size_t> flag = sequence<size_t>::uninitialized(2*m+1);
  auto monoid = pbbslib::addm<size_t>();


  //double edges
  par_for(0, m, kDefaultGranularity, [&] (size_t i) {
    edges[2*i] = edgesIn[i];
    edges[2*i+1] = make_pair(EdgeT(DBTGraph::getSecond(edgesIn,i), DBTGraph::getFirst(edgesIn,i)), edgesIn[i].second);
  });

  vtxNew = computeOffsets<EdgeT, DBTGraph::VtxUpdate>(make_slice(edges), make_slice(vtxMap), flag);

  //count lowD
    par_for(0, 2*m, kDefaultGranularity, [&] (size_t i) {
      flag[i] = G->is_low_v(DBTGraph::getSecond(edges,i));
    });
    par_for(0, vtxNew.size(), kDefaultGranularity, [&] (size_t i) {
      size_t s = vtxNew[i].offset;
      size_t s2 = vtxNew[i].insOffset();
      size_t e = vtxNew[i].end();
      vtxNew[i].insert_low_degree = pbbslib::reduce(flag.cut(s,s2 ), monoid);
      vtxNew[i].delete_low_degree = pbbslib::reduce(flag.cut(s2,e ), monoid);
    });

    flag.clear();
    return vtxNew;
}

template <class Graph, class UT>
void compare(DBTGraph::DyGraph<Graph>* DG, const std::vector<UT>& edges, size_t s, size_t e, size_t n, commandLine& P){
  using W = gbbs::empty;
  using vertex_type = symmetric_vertex<W>;
  using edge_type = vertex_type::edge_type; //std::tuple<uintE, W>
  size_t num_vertices = n;
  auto monoid = pbbslib::addm<size_t>();

  // DBTGraph::SymGraph G = edge_list_to_symmetric_graph(edges, n, 0, e);

  // count new degrees
  vertex_data* vertex_data_array = pbbslib::new_array_no_init<vertex_data>(num_vertices);
  sequence<size_t> newDegrees = sequence<size_t>::from_function(num_vertices, [&](const size_t i) {
    return DG->get_degree(i);
  });
  size_t num_edges = pbbslib::scan_inplace(make_slice(newDegrees), monoid);

  par_for(0, num_vertices-1, kDefaultGranularity, [&](const size_t i) {
    vertex_data_array[i].degree = newDegrees[i+1]-newDegrees[i];
    vertex_data_array[i].offset = newDegrees[i];
  });
  vertex_data_array[num_vertices-1].degree = num_edges - newDegrees[num_vertices-1];
  vertex_data_array[num_vertices-1].offset = newDegrees[num_vertices-1];
  newDegrees.clear();
  // for(size_t i=0; i< num_vertices; ++i) {
  //   if(vertex_data_array[i].degree != G.v_data[i].degree){
  //     cout << "degree wrong! " << i << endl;
  //   }
  // }

  // put edges to array, first old edges, then new edges
  auto edges_seq = pbbslib::new_array_no_init<edge_type>(num_edges);
  for(size_t i=0; i< num_edges; ++i) {
    get<0>(edges_seq[i]) = 345;
  }

  // insert from tables
  auto edges_slice = pbbslib::make_range(edges_seq, num_edges);
  for(uintE u=0; u < num_vertices; ++u) {
    size_t offset = vertex_data_array[u].offset;
    DG->get_neighbors_major(u, edges_slice, offset);
    pbbslib::sample_sort_inplace(edges_slice.slice(offset, offset + vertex_data_array[u].degree),
    [&](const edge_type& a, const edge_type& b) {
      return get<0>(a) < get<0>(b);
    });
  }


  // for(size_t i=0; i< num_edges; ++i) {
  //   if(get<0>(edges_seq[i]) != get<0>(G.e0[i])){
  //     cout << "neighbor wrong! " << i << endl;
  //   }
  // }

  DBTGraph::SymGraph G2 = DBTGraph::SymGraph(
    vertex_data_array,
    num_vertices,
    num_edges,
    [=] () { pbbslib::free_arrays(vertex_data_array, edges_seq); },
    edges_seq);

    auto f = [&] (uintE u, uintE v, uintE w) { };
  size_t c = Triangle_degree_ordering(G2, f);
  cout << c << endl;
}


}
}