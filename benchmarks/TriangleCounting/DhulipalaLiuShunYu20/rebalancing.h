#pragma once

#include <tuple>
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"
#include "dynamic_graph.h"
#include "gbbs/gbbs.h"
#include "pbbslib/monoid.h"
#include "preprocess.h"
#include "shared.h"
#include "sparse_table.h"

using namespace std;

namespace gbbs {
namespace DBTGraph {

// if build_new is false, do not update  DGnew to new graph, used in benchmark
// given edges updates[s,e] and vertices number  n
// build DGnew to be a graph with block_size
// return triangle counts in the graph
template <class UT>
tuple<size_t, DyGraph<SymGraph>*> majorRebalancing(
    const std::vector<UT>& updates, size_t s, size_t e, size_t n,
    size_t block_size, commandLine& P, bool build_new = true) {
  DyGraph<SymGraph>* DGnew;
  timer t;
  t.start();
  SymGraph G = DBTInternal::edge_list_to_symmetric_graph(updates, n, s, e);
  t.next("[MAJ] build G");

  t.start();
  auto f = [&](uintE u, uintE v, uintE w) {};  // count triangles
  size_t c = Triangle_degree_ordering(G, f);
  t.next("[MAJ] count tri");  

  // convert to new grpah
  if (build_new) {
    t.start();
    DGnew = new DyGraph<SymGraph>(block_size, G, n);
    t.next("[MAJ] update DG");
  }

  std::cout << "MAJOR REBALANCE deleting G" << std::endl;
  G.del();

  return make_tuple(c, DGnew);
}

// edges needs to have size 2m
// edges are preprocessed and sorted
// merge edges and DG  to a new graph DGnew
// return new triangle counts
template <class Graph>
tuple<size_t, DyGraph<SymGraph>*> majorRebalancing(
    DyGraph<Graph>* DG, pbbs::sequence<pair<EdgeT, bool>>& edges,
    pbbs::sequence<VtxUpdate>& vtxNew, pbbs::sequence<size_t>& vtxMap,
    size_t num_vertices, commandLine& P) {
  using W = gbbs::empty;
  using vertex_type = symmetric_vertex<W>;
  using edge_type = vertex_type::edge_type;  // std::tuple<uintE, W>
  DyGraph<SymGraph>* DGnew;
  auto monoid = pbbslib::addm<size_t>();

  timer md; md.start();
  if (DG->num_edges() != 0) {
    // mark deletions in tables
    parallel_for(0, vtxNew.size(),
                 [&](size_t i) {
                   DG->markEdgeDeletion(
                       vtxNew[i],
                       edges.slice(vtxNew[i].insOffset(), vtxNew[i].end()));
                 },
                 1);
  }
  md.stop(); md.reportTotal("Major Rebalance: mark deletion time");

  // count new degrees
  vertex_data* vertex_data_array =
      pbbs::new_array_no_init<vertex_data>(num_vertices);
  pbbs::sequence<size_t> newDegrees =
      pbbs::sequence<size_t>(num_vertices, [&](const size_t i) {
        if (vtxMap[i] == EMPTYVMAP) {
          return DG->get_degree(i);
        } else {
          return DG->get_new_degree(vtxNew[vtxMap[i]]);
        }
      });
  size_t num_edges = pbbs::scan_inplace(newDegrees.slice(), monoid);
  if (num_edges == 0) {
    DGnew = new DyGraph<SymGraph>(DG->get_block_size(), num_vertices);
    return make_tuple(0, DGnew);
  }

  par_for(0, num_vertices - 1, DBTGraph::smallTasksForThreshold,
          [&](const size_t i) {
            vertex_data_array[i].degree = newDegrees[i + 1] - newDegrees[i];
            vertex_data_array[i].offset = newDegrees[i];
          });
  vertex_data_array[num_vertices - 1].degree =
      num_edges - newDegrees[num_vertices - 1];
  vertex_data_array[num_vertices - 1].offset = newDegrees[num_vertices - 1];
  newDegrees.clear();

  // put edges to array, first old edges, then new edges
  edge_type* edges_array = pbbs::new_array_no_init<edge_type>(num_edges);
  // insert from arrays
  parallel_for(0, vtxNew.size(),
               [&](const size_t i) {
                 VtxUpdate u = vtxNew[i];
                 size_t offset = vertex_data_array[u.id].offset;
                 offset += DG->get_degree(u.id) - u.delDeg();
                 par_for(0, u.insert_degree, DBTGraph::smallTasksForThreshold,
                         [&](const size_t j) {
                           edges_array[offset + j] = std::make_tuple(
                               getSecond(edges, u.offset + j), gbbs::empty());
                         });
               },
               1);

  pbbs::sequence<edge_type> edges_seq =
      pbbs::sequence<edge_type>(edges_array, num_edges);

  // insert from tables
  if (DG->num_edges() != 0) {
    parallel_for(
        0, num_vertices,
        [&](const size_t u) {
          size_t offset = vertex_data_array[u].offset;
          DG->get_neighbors_major(u, edges_seq.slice(), offset);
          pbbs::sample_sort_inplace(
              edges_seq.slice(offset, offset + vertex_data_array[u].degree),
              [&](const edge_type& a, const edge_type& b) {
                return get<0>(a) < get<0>(b);
              });
        },
        1);
  }

  // make graph, count triangles
  timer sg; sg.start();
  SymGraph G =
      SymGraph(vertex_data_array, num_vertices, num_edges,
               [=]() { pbbslib::free_array(vertex_data_array); },
               edges_array);
  sg.stop(); sg.reportTotal("Major Rebalance: build sym graph time");

  // count triangles
  auto f = [&](uintE u, uintE v, uintE w) {};
  size_t c =
      Triangle_degree_ordering(G, f);

  timer ng; ng.start();
  // convert to new grpah
  size_t block_size = DG->get_block_size();
  DGnew = new DyGraph<SymGraph>(block_size, G, num_vertices);
  // delete G;
  ng.stop(); ng.reportTotal("build new graph time");
  G.del();

  return make_tuple(c, DGnew);
}

/* Perform minor rebalancing and update the degrees arrays and status arrays in
 * DG */
template <class Graph>
size_t minorRebalancing(DyGraph<Graph>* DG, pbbs::sequence<VtxUpdate>& vtxNew,
                        pbbs::sequence<size_t>& vtxMap) {
  size_t n = DG->num_vertices();
  auto monoid = pbbslib::addm<size_t>();
  pbbs::sequence<VtxUpdate> vtxChangeLH = pbbs::sequence<VtxUpdate>();
  pbbs::sequence<VtxRbl> vtxRbl = pbbs::sequence<VtxRbl>();
  pbbs::sequence<size_t> vtxRblMap = pbbs::sequence<size_t>();
  pbbs::sequence<pair<EdgeT, bool>> rblEdges =
      pbbs::sequence<pair<EdgeT, bool>>();
  size_t numLtoH = 0;
  size_t numHtoL = 0;
  size_t newLowNum = DG->num_vertices_low();

  //  ============================= find vertices that change low/high status
  //  =============================
  par_for(0, vtxNew.size(), DBTGraph::smallTasksForThreshold,
          [&](const size_t i) {
            if (DG->change_status(vtxNew[i])) vtxNew[i].change_status = true;
          });
  pbbs::sequence<VtxUpdate> vtxChange = pbbslib::filter(
      vtxNew, [&](const VtxUpdate& u) { return u.change_status; });

  if (vtxChange.size() != 0) {  //  ============================= continue if
                                //  there is changes. Otherwise go to degree
                                //  updates   =============================

    pbbs::sequence<bool> flag = pbbs::sequence<bool>(
        vtxChange.size(),
        [&](const size_t i) { return DG->is_high_v(vtxChange[i].id); });

    // [0,numLtoH) are vertices that change from L to H, the rest from H to L
    pair<pbbs::sequence<VtxUpdate>, size_t> vtxChangeLHsize =
        pbbs::split_two(vtxChange, flag);
    vtxChangeLH = move(vtxChangeLHsize.first);  // LtoH, then HtoL
    numLtoH = vtxChangeLHsize.second;
    numHtoL = vtxChange.size() - numLtoH;

    vtxChange.clear();
    flag.clear();

    //  ============================= update Wedge Table, remove
    //  =============================
    // remove HLH where u was high
    parallel_for(numLtoH, vtxChangeLH.size(),
                 [&](const size_t i) {
                   VtxUpdate u = vtxChangeLH[i];
                   DG->minorRblDeleteWedge(u, vtxNew, vtxMap);
                 },
                 1);

    // remove HLH where center w was low
    parallel_for(0, numLtoH,
                 [&](const size_t i) {  // called after deleting HLH
                   DG->minorRblDeleteWedgeCenter(vtxChangeLH[i]);
                 },
                 1);

    //  ============================= Count Rbled Degrees
    //  =============================
    pbbs::sequence<size_t> newDegrees =
        pbbs::sequence<size_t>(vtxChangeLH.size(), [&](size_t i) {
          return DG->get_new_degree(vtxChangeLH[i]);
        });  // TOCO: can optimize to delayed seq
    size_t rblN = pbbs::scan_inplace(newDegrees.slice(), monoid);
    rblEdges = pbbs::sequence<pair<EdgeT, bool>>::no_init(rblN);
    vtxRblMap = pbbs::sequence<size_t>(n, EMPTYVMAP);
    parallel_for(
        0, vtxChangeLH.size(),
        [&](size_t i) {
          size_t ngh_s = newDegrees[i];
          size_t ngh_e = rblN;
          if (i < vtxChangeLH.size() - 1) ngh_e = newDegrees[i + 1];
          if (i < numLtoH) {
            DG->template get_neighbors_minor<
                pair<EdgeT, bool>, MakeEdgeLtoH<typename DyGraph<Graph>::SetT>>(
                vtxChangeLH[i], rblEdges.slice(), ngh_s, ngh_e, true);
          } else {
            DG->template get_neighbors_minor<
                pair<EdgeT, bool>, MakeEdgeHtoL<typename DyGraph<Graph>::SetT>>(
                vtxChangeLH[i], rblEdges.slice(), ngh_s, ngh_e, false);
          }
        },
        1);
    auto flagg = pbbs::sequence<size_t>();  // dump array
    if (rblN != 0) {
      vtxRbl = DBTInternal::computeOffsets<EdgeT, VtxRbl>(
          rblEdges.slice(), vtxRblMap.slice(), flagg);
    }
    newDegrees.clear();

    par_for(0, vtxNew.size(), DBTGraph::smallTasksForThreshold,
            [&](const size_t i) {  // update degrees and low degrees from
                                   // inserts/deletes to tmp array
              DG->updateDegreesTmp(vtxNew[i]);
            });

    par_for(0, vtxRbl.size(), DBTGraph::smallTasksForThreshold,
            [&](const size_t
                    i) {  // update low degrees from rebalancing  to tmp array
              DG->updateDegreesTmp(vtxRbl[i]);
            });

    //  ============================= Reisze lower table
    //  =============================
    // rezize bottom table for all v (L to H and H to L)
    // vertices using blocks do not need to update lower size table bc low and
    // high ngh are stored together
    // TODO: resize to account for delete edges

    parallel_for(
        0, vtxRbl.size(),
        [&](size_t i) {
          uintE v = vtxRbl[i].id;
          if (DG->use_block_v(v)) return;
          size_t deltaLtoH = vtxRbl[i].LtoH;
          size_t deltaHtoL = vtxRbl[i].getHtoL();
          VtxUpdate vobj = VtxUpdate(v);
          if (vtxMap[v] != EMPTYVMAP) vobj = vtxNew[vtxMap[v]];
          if (deltaLtoH > deltaHtoL)
            DG->minorRblResizeBottomTable(vobj, deltaLtoH - deltaHtoL, true);
          if (deltaHtoL > deltaLtoH)
            DG->minorRblResizeBottomTable(vobj, deltaHtoL - deltaLtoH, false);
        },
        1);

    //  ============================= Move between lower tables
    //  =============================

    // delete from bottom table
    parallel_for(0, vtxRbl.size(),
                 [&](const size_t i) {
                   VtxRbl v = vtxRbl[i];
                   DG->minorRblMoveBottomTable(
                       v.id, rblEdges.slice(v.offset, v.end()), true);
                 },
                 1);

    // insert to bottom table
    parallel_for(0, vtxRbl.size(),
                 [&](const size_t i) {
                   VtxRbl v = vtxRbl[i];
                   DG->minorRblMoveBottomTable(
                       v.id, rblEdges.slice(v.offset, v.end()), false);
                 },
                 1);

    //  ============================= Reisze top table
    //  =============================
    newLowNum = DG->minorRblResizeTop(numHtoL, numLtoH);

    //  ============================= Move between top tables
    //  =============================
    
#ifdef DBT_USING_ARRAYTOP
    // move top table LtoH
    par_for(0, numLtoH, DBTGraph::smallTasksForThreshold, [&](const size_t i) {
      DG->minorRblMoveTopTable(vtxChangeLH[i], true);
    });
    // move top table HtoL
    par_for(numLtoH, vtxChangeLH.size(), DBTGraph::smallTasksForThreshold,
            [&](const size_t i) {
              DG->minorRblMoveTopTable(vtxChangeLH[i], false);
            });
#else
    // move top table LtoH
    par_for(0, numLtoH, DBTGraph::smallTasksForThreshold, [&](const size_t i) {
      DG->minorRblMoveTopTable(vtxChangeLH[i], true, false);
    });

    par_for(0, numLtoH, DBTGraph::smallTasksForThreshold,
            [&](const size_t i) {  // update  status
              DG->minorRblMoveTopTable(vtxChangeLH[i], true, true);
            });

    // move top table HtoL
    par_for(numLtoH, vtxChangeLH.size(), DBTGraph::smallTasksForThreshold,
            [&](const size_t i) {
              DG->minorRblMoveTopTable(vtxChangeLH[i], false, false);
            });

    par_for(numLtoH, vtxChangeLH.size(), DBTGraph::smallTasksForThreshold,
            [&](const size_t i) {  // update  status
              DG->minorRblMoveTopTable(vtxChangeLH[i], false, true);
            });
#endif

  }  // end if vtxChange.size() != 0;

  //  ============================= Update Degrees =============================
  timer t;
  t.start();  //  update degrees

  // update degree and low degree and lowNum,
  par_for(0, vtxNew.size(), DBTGraph::smallTasksForThreshold,
          [&](const size_t
                  i) {  // update degrees and low degrees from inserts/deletes
            DG->updateDegrees(vtxNew[i]);
          });

  par_for(0, vtxRbl.size(), DBTGraph::smallTasksForThreshold,
          [&](const size_t i) {  // update low degrees from rebalancing
            DG->updateDegrees(vtxRbl[i]);
          });
  // #ifdef DBT_TOMB_MERGE
  //     par_for(0, vtxNew.size(), [&] (const size_t i) { // pack table to
  //     arrays if new degree is low enough, delete tables and change status
  //       DG->downSizeTablesBoth(vtxNew[i]);
  //     });
  // #else
  parallel_for(
      0, vtxNew.size(),
      [&](const size_t i) {  // pack table to arrays if new degree is low
                             // enough, called before downSizeTablesDeletes
        DG->downSizeTables(vtxNew[i]);
      },
      1);
  par_for(0, vtxNew.size(), DBTGraph::smallTasksForThreshold,
          [&](const size_t i) {  // delete tables packed and change status
            DG->downSizeTablesDeletes(vtxNew[i]);
          });
  // #endif
  DG->set_vertices_low(newLowNum);

  t.stop();
  t.reportTotal("8. update degrees");

  //  ============================= Update Wedge Table, L to H
  //  =============================
  // insert HLH where u was low
  parallel_for(0, numLtoH,
               [&](const size_t i) {
                 VtxUpdate u = vtxChangeLH[i];
                 if (DG->use_block_v(u.id)) return;
                 DG->minorRblInsertWedge(u, vtxNew, vtxMap);
               },
               1);

  // insert HLH where center w was high
  // called after minorRblInsertWedge
  // only add to wegdes(u,v) where both (u,v) are not changing from L to H
  // otherwise the count already includes new low centers
  parallel_for(numLtoH, vtxChangeLH.size(),
               [&](const size_t i) {  // called after inserting HLH
                 DG->minorRblInsertWedgeCenter(vtxChangeLH[i], vtxNew, vtxMap);
               },
               1);

  // flag.clear();
  // vtxChange.clear();

  vtxChangeLH.clear();
  rblEdges.clear();
  vtxRbl.clear();
  vtxRblMap.clear();

  return newLowNum;
}
}
}
