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


    template <class Graph>
    class UpdateDyGraph{
        size_t n;
        pbbs::sequence<size_t> vtxMap; //n
        pbbs::sequence<bool> newStatus; //nn
        pbbs::sequence<pair<uintE, size_t>> vtxNew; //nn



        UpdateDyGraph(DBTGraph::DyGraph<Graph>& DG){
            n = DG.num_vertices();
            newStatus = pbbs::sequence<bool>::no_init(n);


        }

// edges needs to have size 2m
template <class Graph, class EdgeT>
inline pbbs::sequence<DBTGraph::VtxUpdate> minorRebalancing(DBTGraph::DyGraph<Graph>& DG, pbbs::range<pair<EdgeT,bool>> &updates){
    // count top level resizing
    flag.shrink(2 * numVtx); // to H, then to L, 
    par_for(0, numVtx, [&] (size_t i) {flag[i] = 0;}
    par_for(0, numVtx, [&] (size_t i) {
      uintE  u = vtxMap[i].first;
      size_t d = newD[u];

      bool status = DG.status_v(u);//true if high
      if(DG.change_status(u, d)){
        newStatus[u] = !status;
        if(status){flag[numVtx + i] = 1;} else{flag[i] = 1;}
      }else{
        newStatus[u] =  status;
      }
    });

    cilk_spawn size_t toH = pbbslib::reduce(flag.slice(0, numVtx), monoid);
    size_t toL = pbbslib::reduce(flag.slice(numVtx, 2*numVtx), monoid);
    cilk_sync;

    //todo: tune top level load
    size_t delta = toH-toL;
    if(delta < 0 && DG.LL->size() < 1.1*(lowNum - delta)){
      // resize LL and LH top level
      size_t newLowNum = lowNum - delta;

      // insert changes

    }else if(delta > 0 && DG.HL->size() < 1.1*(n - lowNum + delta)){
     // resize HL and HH top level

      // insert changes
    }  

  // update graph D, lowD, lowNum
  flag.clear();

}
    };


}
}