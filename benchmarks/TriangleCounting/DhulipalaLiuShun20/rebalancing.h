#pragma once

#include <tuple>
#include "gbbs/gbbs.h"
#include "pbbslib/monoid.h"
#include "sparse_table.h"
// #include "gbbs/macros.h"
#include "shared.h"
#include "dynamic_graph.h"
#include "preprocess.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"

using namespace std;

namespace gbbs{
namespace DBTGraph{

// if DGnew is NULL, do not generate new DyGraph
template <class UT>
size_t majorRebalancing(const std::vector<UT>& updates, size_t s, size_t e, size_t n, size_t block_size, DyGraph<SymGraph> &DGnew, commandLine& P, bool build_new = true){

  // vector<gbbs::gbbs_io::Edge<pbbs::empty>> edges = getEdgeVec(updates, s, e);
  SymGraph G = DBTInternal::edge_list_to_symmetric_graph(updates, n, s, e);

  // count triangles
  auto f = [&] (uintE u, uintE v, uintE w) { };
  size_t c = Triangle(G, f, "degree", P);  //TODO: which ordering?, how to ini commandline object?

  // convert to new grpah
  if(build_new) DGnew = DyGraph(block_size, G, n);
  // delete G;
  
  return c;
}

// edges needs to have size 2m
// edges are preprocessed and sorted
template <class Graph>
size_t majorRebalancing(DyGraph<Graph>& DG, DyGraph<SymGraph> &DGnew, pbbs::sequence<pair<EdgeT,bool>> &edges, 
                              pbbs::sequence<VtxUpdate> &vtxNew,  pbbs::sequence<size_t> &vtxMap, size_t num_vertices, commandLine& P){

  using W = pbbslib::empty;
  using vertex_type = symmetric_vertex<W>;
  using edge_type = vertex_type::edge_type; //std::tuple<uintE, W>

  auto monoid = pbbslib::addm<size_t>();

  if( DG.num_edges()  !=  0){
    //mark deletions in tables
    par_for(0, vtxNew.size(), [&] (size_t i) {
      DG.markEdgeDeletion(vtxNew[i],  edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
    });
  }


  // count new degrees
  vertex_data *vertex_data_array = pbbs::new_array_no_init<vertex_data>(num_vertices);
  pbbs::sequence<size_t> newDegrees = pbbs::sequence<size_t>(num_vertices, [&](const size_t i) {
    if(vtxMap[i] == EMPTYVMAP){
      return DG.get_degree(i);
    }else{
      return DG.get_new_degree(vtxNew[vtxMap[i]]);
    }
  }); 
  size_t num_edges = pbbs::scan_inplace(newDegrees.slice(), monoid);  

  par_for(0, num_vertices-1, [&](const size_t i) {
    vertex_data_array[i].degree = newDegrees[i+1]-newDegrees[i];
    vertex_data_array[i].offset = newDegrees[i];
  });    
  vertex_data_array[num_vertices-1].degree = num_edges - newDegrees[num_vertices-1];
  vertex_data_array[num_vertices-1].offset = newDegrees[num_vertices-1];
  newDegrees.clear();

  // put edges to array, first old edges, then new edges
  edge_type *edges_array = pbbs::new_array_no_init<edge_type>(num_edges);
  // insert from arrays
  par_for(0, vtxNew.size(), [&] (const size_t i) {
    VtxUpdate u =  vtxNew[i];
    size_t offset = vertex_data_array[u.id].offset;
    offset += DG.get_degree(u.id) - u.delDeg();
    par_for(0, u.insert_degree, [&](const size_t j){
      edges_array[offset + j] = std::make_tuple(getSecond(edges, u.offset + j), pbbs::empty());
    });
  });

  pbbs::sequence<edge_type> edges_seq = pbbs::sequence<edge_type>(edges_array, num_edges);

  // insert from tables 
  if( DG.num_edges()  !=  0){
  par_for(0, num_vertices, [&](const size_t u) {
    size_t offset = vertex_data_array[u].offset;
    DG.get_neighbors_major(u, edges_seq, offset);
    pbbs::sample_sort_inplace(edges_seq.slice(offset, offset + vertex_data_array[u].degree), 
    [&](const edge_type& a, const edge_type& b) {
      return get<0>(a) < get<0>(b);
    });
  });}



  // make graph, count triangles
  SymGraph G = SymGraph(
    vertex_data_array,
    num_vertices,
    num_edges,
    [=] () { pbbslib::free_arrays(vertex_data_array, edges_array); },
    edges_array);

  // count triangles
  auto f = [&] (uintE u, uintE v, uintE w) { };
  size_t c = Triangle(G, f, "degree", P);  //TODO: which ordering?, how to ini commandline object?

  // convert to new grpah
  size_t block_size = DG.get_block_size();
  DGnew = DyGraph(block_size, G, num_vertices);
  // delete G;
  
  return c;
}

/* Perform minor rebalancing and update the degrees arrays and status arrays in DG */
template <class Graph>
size_t minorRebalancing(DyGraph<Graph>& DG, pbbs::sequence<VtxUpdate>& vtxNew, pbbs::sequence<size_t>& vtxMap){
    size_t n = DG.num_vertices();
    // size_t numVtx = vtxNew.size();
    auto monoid = pbbslib::addm<size_t>();

    //  ============================= find vertices that change low/high status   =============================
    //TODO: can optimize
    pbbs::sequence<VtxUpdate> vtxChange = pbbslib::filter(vtxNew, [&](const VtxUpdate &u){
      return DG.change_status(u);
    });
    if(vtxChange.size() == 0) return DG.num_vertices_low();
    pbbs::sequence<bool> flag = pbbs::sequence<bool>(vtxChange.size(), [&](const size_t i){
      return DG.is_high_v(vtxChange[i].id);
    });

    // [0,numLtoH) are vertices that change from L to H, the rest from H to L
    pair<pbbs::sequence<VtxUpdate>, size_t> vtxChangeLHsize = pbbs::split_two(vtxChange, flag);
    pbbs::sequence<VtxUpdate> vtxChangeLH =  vtxChangeLHsize.first; // LtoH, then HtoL 
    size_t numLtoH =  vtxChangeLHsize.second;
    size_t numHtoL =  vtxChange.size() - numLtoH;

    vtxChange.clear();
    flag.clear();

    //  ============================= update Wedge Table, H to L ============================= 
    par_for(numLtoH, vtxChangeLH.size(), [&] (const size_t i) { 
      VtxUpdate u = vtxChangeLH[i];
      DG.minorRblDeleteWedge(u, vtxNew, vtxMap);
    });

    //  ============================= Count Rbled Degrees =============================
    pbbs::sequence<size_t> newDegrees = pbbs::sequence<size_t>(vtxChangeLH.size(), [&](size_t i){return DG.get_new_degree(vtxChangeLH[i]);}); //TOCO: can optimize to delayed seq
    size_t rblN = pbbs::scan_inplace(newDegrees.slice(), monoid);   
    pbbs::sequence<pair<EdgeT,bool>> rblEdges = pbbs::sequence<pair<EdgeT,bool>>::no_init(rblN); 
    pbbs::sequence<VtxRbl> vtxRbl = pbbs::sequence<VtxRbl>();
    pbbs::sequence<size_t> vtxRblMap = pbbs::sequence<size_t>(n, EMPTYVMAP);
    par_for(0,vtxChangeLH.size(),[&](size_t i){
      size_t ngh_s = newDegrees[i];
      size_t ngh_e = rblN;
      if(i < vtxChangeLH.size()-1) ngh_e = newDegrees[i+1];
      if(i < numLtoH){
      DG.template get_neighbors_minor<pair<EdgeT,bool>, MakeEdgeLtoH<typename DyGraph<Graph>::SetT>>(vtxChangeLH[i], rblEdges, ngh_s, ngh_e, true); 
      }else{
      DG.template get_neighbors_minor<pair<EdgeT,bool>, MakeEdgeHtoL<typename DyGraph<Graph>::SetT>>(vtxChangeLH[i], rblEdges, ngh_s, ngh_e, false); 
      }
    });
    DBTInternal::computeOffsets<EdgeT, VtxRbl>(rblEdges, vtxRbl, vtxRblMap);

    //  ============================= Reisze lower table ============================= 
    // rezize bottom table for all v (L to H and H to L)
    // vertices using blocks do not need to update lower size table bc low and high ngh are stored together
    // TODO: resize to account for delete edges

    par_for(0,vtxRbl.size(),[&](size_t i){
      uintE v = vtxRbl[i].id;
      if(DG.use_block_v(v)) return;
      size_t deltaLtoH = vtxRbl[i].LtoH;
      size_t deltaHtoL = vtxRbl[i].getHtoL();
      VtxUpdate vobj = VtxUpdate(v);
      if(vtxMap[v] != EMPTYVMAP) vobj = vtxNew[vtxMap[v]];
      if(deltaLtoH > deltaHtoL) DG.minorRblResizeBottomTable(vobj, deltaLtoH - deltaHtoL, true);
      if(deltaHtoL > deltaLtoH) DG.minorRblResizeBottomTable(vobj, deltaHtoL - deltaLtoH, false);
    });

    //  ============================= Move between lower tables ============================= 

    //delete from bottom table
    par_for(0,vtxRbl.size(),[&](const size_t i){
      VtxRbl v = vtxRbl[i];
      DG.minorRblMoveBottomTable(v.id, rblEdges.slice(v.offset, v.end()), DG.is_low_v(v.id), true);
    });

    //insert to bottom table
    par_for(0,vtxRbl.size(), [&](const size_t i){
      VtxRbl v = vtxRbl[i];
      DG.minorRblMoveBottomTable(v.id, rblEdges.slice(v.offset, v.end()), DG.is_low_v(v.id), false);
    });


    //  ============================= Reisze top table ============================= 
    size_t newLowNum = DG.minorRblResizeTop(numHtoL, numLtoH);

    //  ============================= Move between top tables ============================= 
    // move top table LtoH
    par_for(0,  numLtoH, [&] (const size_t i) { 
      DG.minorRblMoveTopTable(vtxChangeLH[i], true, false);});

    par_for(0,  numLtoH, [&] (const size_t i) { //update  status
      DG.minorRblMoveTopTable(vtxChangeLH[i], true, true);});

    // move top table HtoL
    par_for(numLtoH,  vtxChangeLH.size(), [&] (const size_t i) { 
      DG.minorRblMoveTopTable(vtxChangeLH[i], false, false);});

    par_for(numLtoH,  vtxChangeLH.size(), [&] (const size_t i) { //update  status
      DG.minorRblMoveTopTable(vtxChangeLH[i], false, true);});

    //  ============================= Update Degrees ============================= 
    timer t;
    t.start(); //  update degrees


    //update degree and low degree and lowNum,
    
    par_for(0, vtxNew.size(), [&] (const size_t i) { // update degrees and low degrees from inserts/deletes
      DG.updateDegrees(vtxNew[i]);
    });

    par_for(0, vtxRbl.size(), [&] (const size_t i) { // update low degrees from rebalancing
      DG.updateDegrees(vtxRbl[i]);
    });

    par_for(0, vtxNew.size(), [&] (const size_t i) { // pack table to arrays if new degree is low enough, called before downSizeTablesDeletes
      DG.downSizeTables(vtxNew[i]);
    });
    par_for(0, vtxNew.size(), [&] (const size_t i) { // delete tables packed and change status
      DG.downSizeTablesDeletes(vtxNew[i]);
    });

    DG.set_vertices_low(newLowNum);

    t.stop();t.reportTotal("8. update degrees");


    //  ============================= Update Wedge Table, L to H ============================= 
    par_for(0, numLtoH, [&] (const size_t i) { 
      VtxUpdate u = vtxChangeLH[i];
      if(DG.use_block_v(u.id)) return;
      DG.minorRblInsertWedge(u);
    });

  // flag.clear();
  // vtxChange.clear();
  vtxChangeLH.clear();
  rblEdges.clear(); 
  vtxRbl.clear();
  vtxRblMap.clear();

  return newLowNum;

}

}}