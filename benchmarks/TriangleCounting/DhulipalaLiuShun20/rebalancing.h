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

using namespace std;

namespace gbbs{
namespace DBTGraph{

// edges needs to have size 2m
template <class Graph, class EdgeT>
inline size_t minorRebalancing(DBTGraph::DyGraph<Graph>& DG, pbbs::range<pair<EdgeT,bool>> &updates, 
                              pbbs::sequence<DBTGraph::VtxUpdate> vtxNew, pbbs::sequence<size_t> vtxMap){
    size_t n = DG.num_vertices();
    size_t numVtx = vtxNew.size();
    auto monoid = pbbslib::addm<size_t>();
    pbbs::sequence<bool> newStatus = pbbs::sequence<bool>::no_init(n);

    //todo: can optimize
    pbbs::sequence<DBTGraph::VtxUpdate> vtxChange = pbbslib::filter(vtxNew, [&](DBTGraph::VtxUpdate u){
      return DG.change_status(u.id, u.insert_degree, u.delDeg());
    });
    pbbs::sequence<bool> flag = pbbs::sequence<bool>(vtxChange.size(), [&](DBTGraph::VtxUpdate u){
      return DG.is_high_v(u.id);
    })

    pair<pbbs::sequence<DBTGraph::VtxUpdate>, size_t> vtxChangeLHsize = pbbslib::split_two(vtxChange, flag);
    pbbs::sequence<DBTGraph::VtxUpdate> vtxChangeLH =  vtxChangeLHsize.first; // LtoH, then HtoL 
    size_t numLtoH =  vtxChangeLHsize.second;
    size_t numHtoL =  vtxChange.size() - numL;

    vtxChange.clear();
    flag.clear();

    pbbs::sequence<size_t> newDegrees = pbbs::sequence<size_t>::no_init(vtxChangeLH.size());

    par_for(0,vtxChangeLH.size(),[&](size_t i){
      newDegrees[i] = DG.get_new_degree(vtxChangeLH[i]);
    });

    //reisze lower table
    size_t nghN = pbbs::scan_inplace(newDegrees,monoid);
    pbbs::sequence<EdgeT> Ngh = pbbs::sequence<EdgeT>::no_init(nghN);
    par_for(0,vtxChangeLH.size(),[&](size_t i){
      size_t ngh_s = newDegrees[i];
      size_t ngh_e = nghN;
      if(i < vtxChangeLH.size()-1) ngh_e = newDegrees[i+1];
      DG.get_neighbors(vtxChangeLH[i], Ngh, ngh_s, ngh_e, i<numLtoH);
    });
    pbbs::sample_sort_inplace(Ngh.slice(0,numLtoH),less<EdgeT>);
    pbbs::sample_sort_inplace(Ngh.slice(numLtoH,Ngh.size()),less<EdgeT>);


    pbbs::sequence<size_t> flag2 = pbbs::sequence<size_t>(nghN, 0);
    

    par_for(0,nghN-1,[&](size_t i){
      if(Ngh[i].first!=Ngh[i+1].first){flag2[i] = 1;}
    });
    flag2[numLtoH-1] = 1;
    flag2[nghN-1] = 1;
    size_t offsetsLength =  pbbs::scan_inplace(flag2.slice(), monoid);
    size_t HtoLOffset = flag2[numLtoH];
    pbbs::sequence<size_t> offsetsNgh = pbbs::sequence<size_t>::no_init(offsetsLength+1);


    offsetsNgh[0] = 0; 
    par_for(1,nghN,[&](size_t i){
      if(flag2[i]!=flag2[i-1]){
        offsetsNgh[flag2[i]] = i;
      }
    });
    offsetsNgh[nghN] = nghN;

    // rezize bottom table 
    par_for(0,offsetsLength,[&](size_t i){
      uintE v = Ngh[offsetsNgh[i]].first;
      if(use_block_v(v.id)) return;
      size_t delta = offsetsNgh[i+1] - offsetsNgh[i];
      DG.minorRblResizeBottomTable(vtxNew[vertexMap[v]], delta, i < HtoLOffset);
    });

    //delete from bottom table LtoH
    par_for(0,offsetsLength,[&](size_t i){
      uintE v = Ngh[offsetsNgh[i]].first;
      if(use_block_v(v.id)) return;
      uintE u = Ngh[offsetsNgh[i]].second;
      DG.minorRblMoveBottomTable(vtxNew[vertexMap[v]], Ngh.slice(offsetsNgh[i], offsetsNgh[i+1]), i < HtoLOffset, true)
    });

    //insert to bottom table HtoL
    par_for(0,offsetsLength, [&](size_t i){
      uintE v = Ngh[offsetsNgh[i]].first;
      if(use_block_v(v.id)) return;
      uintE u = Ngh[offsetsNgh[i]].second;
      DG.minorRblMoveBottomTable(vtxNew[vertexMap[v]], Ngh.slice(offsetsNgh[i], offsetsNgh[i+1]), i < HtoLOffset, false)
    });


    //resize top table
    size_t newLowNum = DG.minorRblResizeTop(numHtoL, numLtoH);

    // move top table LtoH
    par_for(0, numLtoH, [&] (size_t i) { //update  status
      DBTGraph::VtxUpdate u = vtxChangeLH[i];
      if(use_block_v(u.id)) return;
      DG.minorRblMoveTopTable(u);
    });

    // move top table HtoL
    par_for(numLtoH, vtxChangeLH.size(), [&] (size_t i) { //update  status
      DBTGraph::VtxUpdate u = vtxChangeLH[i];
      if(use_block_v(u.id)) return;
      DG.minorRblMoveTopTable(u);
    });

  // flag.clear();
  newStatus.clear();
  // vtxChange.clear();
  vtxChangeLH.clear();
  newDegrees.clear();
  flag2.clear();
  offsetsNgh.clear();

  return newLowNum;

}
    };


}