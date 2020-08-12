#pragma once

#include <tuple>
#include "gbbs/gbbs.h"
#include "pbbslib/monoid.h"
// #include "gbbs/pbbslib/sparse_table.h"
#include "sparse_table.h"
#include "set.h"
#include "gbbs/macros.h"
#include "dynamic_graph.h"
#include "shared.h"
#include "preprocess.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"

using namespace std;

namespace gbbs{
namespace DBTGraph{

// edges needs to have size 2m
// edges are preprocessed and sorted
template <class Graph, class EdgeT>
DyGraph<Graph> majorRebalancing(DyGraph<Graph>& DG, pbbs::sequence<pair<EdgeT,bool>> &edges, 
                              pbbs::sequence<VtxUpdate> vtxNew,  pbbs::sequence<size_t> vtxMap, commandLine& P){
  // using W = G::weight_type;
  // using vertex_type = Graph::vertex_type;
  using W = pbbslib::empty;
  using vertex_type = symmetric_vertex<W>;
  using edge_type = vertex_type::edge_type; //std::tuple<uintE, W>

  size_t num_vertices = DG.num_vertices();
  auto monoid = pbbslib::addm<size_t>();

  //mark deletions in tables
  par_for(0, vtxNew.size(), [&] (size_t i) {
    DG.markEdgeDeletion(vtxNew[i],  edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
  });

  // count new degrees
  pbbs::sequence<vertex_data> *vertex_data_array = pbbs::new_array_no_init<vertex_type>(num_vertices);
  pbbs::sequence<size_t> newDegrees = pbbs::sequence<size_t>(num_vertices, [&](const size_t i) {
    if(vertexMap[i] == EMPTYVMAP){
      return DG.get_degree(i);
    }else{
      return DG.get_new_degree(vtxNew[vertexMap[i]]);
    }
  }); 
  size_t num_edges = pbbs::scan_inplace(newDegrees, monoid);  

  par_for(0, num_vertices-1, [&](const size_t i) {
    vertex_data_array[i].degree = newDegrees[i+1]-newDegrees[i];
    vertex_data_array[i].offset = newDegrees[i];
  });    
  vertex_data_array[num_vertices-1].degree = num_edges - newDegrees[num_vertices-1];
  newDegrees.clear();

  // put edges to array, first old edges, then new edges
  edge_type *edges_array = pbbs::new_array_no_init<edge_type>(num_edges);
  // insert from arrays
  par_for(0, vtxNew.size(), [&] (const size_t i) {
    VtxUpdate u =  vtxNew[i];
    size_t offset = vertex_data_array[u.id].offset;
    offset += DG.get_degree(u.id) - u.delDeg();
    par_for(0, u.insert_degree, [&](const size_t j){
      edges_array[offset + j] = edges[u.offset + j];
    });
  });

  // insert from tables 
  par_for(0, num_vertices, [&](const size_t u) {
    DG.get_neighbors_major(u, edges_array, vertex_data_array[u].offset);
  }

  // make graph, count triangles
  Graph G = Graph(
    vertex_data_array,
    num_vertices,
    num_edges,
    [=] () { pbbslib::free_arrays(vertex_data_array, edges_array); },
    edges_array);

  // count triangles
  auto f = [&] (uintE u, uintE v, uintE w) { };
  size_t c = Triangle(G, f, "degree", P);  //TODO: which ordering?, how to ini commandline object?
  G.clear();

  // convert to new grpah
  DyGraph<Graph> newDG = DyGraph(DG.get_block_size(), G);
  DG.clear();
  
  return newDG;
}

/* Perform minor rebalancing and update the degrees arrays and status arrays in DG */
template <class Graph, class EdgeT>
size_t minorRebalancing(DyGraph<Graph>& DG, pbbs::sequence<VtxUpdate>& vtxNew, pbbs::sequence<size_t>& vtxMap){
    size_t n = DG.num_vertices();
    size_t numVtx = vtxNew.size();
    auto monoid = pbbslib::addm<size_t>();

    //  ============================= find vertices that change low/high status   =============================
    //TODO: can optimize
    pbbs::sequence<VtxUpdate> vtxChange = pbbslib::filter(vtxNew, [&](VtxUpdate u){
      return DG.change_status(u);
    });
    pbbs::sequence<bool> flag = pbbs::sequence<bool>(vtxChange.size(), [&](VtxUpdate u){
      return DG.is_high_v(u.id);
    });

    // [0,numLtoH) are vertices that change from L to H, the rest from H to L
    pair<pbbs::sequence<VtxUpdate>, size_t> vtxChangeLHsize = pbbslib::split_two(vtxChange, flag);
    pbbs::sequence<VtxUpdate> vtxChangeLH =  vtxChangeLHsize.first; // LtoH, then HtoL 
    size_t numLtoH =  vtxChangeLHsize.second;
    size_t numHtoL =  vtxChange.size() - numLtoH;

    vtxChange.clear();
    flag.clear();

    //  ============================= update Wedge Table, H to L ============================= 
    par_for(numLtoH, vtxChangeLH.size(), [&] (size_t i) { 
      VtxUpdate u = vtxChangeLH[i];
      DG.minorRblDeleteWedge(u));
    });

    //  ============================= Count Rbled Degrees =============================
    pbbs::sequence<size_t> newDegrees = pbbs::sequence<size_t>(vtxChangeLH.size(), [&](size_t i){return DG.get_new_degree(vtxChangeLH[i]);}); //TOCO: can optimize to delayed seq
    size_t rblN = pbbs::scan_inplace(newDegrees, monoid);   
    pbbs::sequence<pair<EdgeT,bool>> rblEdges = pbbs::sequence<pair<EdgeT,bool>>::no_init(rblN); 
    pbbs::sequence<VtxRbl> vtxRbl;
    pbbs::sequence<size_t> vtxRblMap = pbbs::sequence<size_t>(n, EMPTYVMAP);
    par_for(0,vtxChangeLH.size(),[&](size_t i){
      size_t ngh_s = newDegrees[i];
      size_t ngh_e = rblN;
      if(i < vtxChangeLH.size()-1) ngh_e = newDegrees[i+1];
      if(i < numLtoH){
      DG.get_neighbors_minor<pair<EdgeT,bool>, MakeEdgeLtoH<DyGraph<Graph>::SetT>>(vtxChangeLH[i], rblEdges, ngh_s, ngh_e, i<numLtoH, false); 
      }else{
      DG.get_neighbors_minor<pair<EdgeT,bool>, MakeEdgeHtoL<DyGraph<Graph>::SetT>>(vtxChangeLH[i], rblEdges, ngh_s, ngh_e, i<numLtoH, false); 
      }
    });
    DBTInternal::computeOffsets<EdgeT, VtxRbl>(rblEdges, vtxRbl, vtxRblMap);
    size_t tblVtxN = vtxRbl.size();

    //  ============================= Reisze lower table ============================= 
    // rezize bottom table for all v (L to H and H to L)
    // vertices using blocks do not need to update lower size table bc low and high ngh are stored together
    // TODO: resize to account for delete edges

    par_for(0,tblVtxN,[&](size_t i){
      uintE v = vtxRbl[i].id;
      if(DG.use_block_v(v)) return;
      size_t deltaLtoH = vtxRbl[i].LtoH;
      size_t deltaHtoL = vtxRbl[i].getHtoL();
      VtxUpdate vobj;
      if(vertexMap[v] == EMPTYVMAP){vobj = VtxUpdate(v);
      }else{vobj = vtxNew[vertexMap[v]];}
      if(deltaLtoH > 0) DG.minorRblResizeBottomTable(vobj, deltaLtoH, true);
      if(deltaHtoL > 0) DG.minorRblResizeBottomTable(vobj, deltaHtoL, false);
    });

    //  ============================= Move between lower tables ============================= 

    //delete from bottom table
    par_for(0,tblVtxN,[&](size_t i){
      VtxRbl v = vtxRbl[i];
      DG.minorRblMoveBottomTable(v.id, rblEdges.slice(v.offset, v.end()), DG.is_low_v(v.id), true)
    });

    //insert to bottom table
    par_for(0,tblVtxN, [&](size_t i){
      VtxRbl v = vtxRbl[i];
      DG.minorRblMoveBottomTable(v.id, rblEdges.slice(v.offset, v.end()), DG.is_low_v(v.id), false)
    });


    //  ============================= Reisze top table ============================= 
    size_t newLowNum = DG.minorRblResizeTop(numHtoL, numLtoH);

    //  ============================= Move between top tables ============================= 
    // move top table LtoH
    par_for(0,  numLtoH, [&] (size_t i) { //update  status
      VtxUpdate u = vtxChangeLH[i];      
      DG.minorRblMoveTopTable(u, true, false);
    });

    par_for(0,  numLtoH, [&] (size_t i) { //update  status
      VtxUpdate u = vtxChangeLH[i];      
      DG.minorRblMoveTopTable(u, true, true);
    });

    // move top table HtoL
    par_for(numLtoH,  vtxChangeLH.size(), [&] (size_t i) { //update  status
      VtxUpdate u = vtxChangeLH[i];      
      DG.minorRblMoveTopTable(u, false, false);
    });

    par_for(numLtoH,  vtxChangeLH.size(), [&] (size_t i) { //update  status
      VtxUpdate u = vtxChangeLH[i];      
      DG.minorRblMoveTopTable(u, false, true);
    });

    //  ============================= Update Degrees ============================= 
    t.start(); //  update degrees


    //update degree and low degree and lowNum,
    
    par_for(0, vtxNew.size(), [&] (size_t i) { // update degrees and low degrees from inserts/deletes
      DG.updateDegrees(vtxNew[i]);
    });

    par_for(0, vtxRbl.size(), [&] (size_t i) { // update low degrees from rebalancing
      DG.updateDegrees(vtxRbl[i]);
    });

    par_for(0, vtxNew.size(), [&] (size_t i) { // pack table to arrays if new degree is low enough
      DG.downSizeTables(vtxNew[i]);
    });
    par_for(0, vtxNew.size(), [&] (size_t i) { // delete tables packed
      DG.downSizeTablesDeletes(vtxNew[i]);
    });

    DG.set_vertices_low(newLowNum);

    t.stop();t.reportTotal("8. update degrees");


    //  ============================= Update Wedge Table, L to H ============================= 
    par_for(0, numLtoH, [&] (size_t i) { 
      VtxUpdate u = vtxChangeLH[i];
      if(use_block_v(u.id)) return;
      DG.minorRblInsertWedge(u));
    });

  // flag.clear();
  // vtxChange.clear();
  vtxChangeLH.clear();
  Ngh.clear();
  offsetsNgh.clear();
  rblEdges.clear(); 
  vtxRbl.clear();
  vtxRblMap.clear();

  return newLowNum;

}

}}