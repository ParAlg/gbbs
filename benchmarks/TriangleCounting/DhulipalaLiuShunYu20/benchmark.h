#ifndef SHARED_H
#define SHARED_H

// *************************************************************
//   Misc
// *************************************************************
#include <string>
#include <iostream>
#include <tuple>
#include "gbbs/gbbs.h"

#include "shared.h"
#include "rebalancing.h"

namespace gbbs{

namespace DBTInternal {

  // template <class EdgeT>
  // pbbs::sequence<pair<EdgeT, bool>> generateEdgeUpdates(size_t nv, size_t numEdges, commandLine& P);

  // ///// PRINTING UTILITIES //////
  // inline void PrintCaption(std::string t_in);

  // inline void PrintSubcaption(std::string t_in);

  // template <class T>
  // inline void PrintFunctionTitle(std::string t_func, T t_suffix);

  // template <class T>
  // inline void PrintFunctionItem(std::string t_func, std::string t_item, T t_suffix);

  inline void PrintBreak() {
    std::cout << std::endl;
  }

  // template <typename A>
  // inline void PrintVec(A *t_vec, int t_len);

  // template <typename A>
  // inline void PrintVec(A t_vec, int t_len);

  // template <typename A>
  // inline void PrintPairSeq(A t_vec);

  // template <typename A>
  // inline void PrintVec2(A *t_vec, int t_len);


  template <class T>
  inline void PrintFunctionItem(std::string t_func, std::string t_item, T t_suffix) {
    std::cout << "[" << t_func << "] " << t_item << " = " << t_suffix << std::endl;
  }

  //edges already randomly shuffled
  template <class UT>
  inline size_t staticCount(const std::vector<UT>& edges, size_t batch_size, commandLine& P, size_t n, size_t batch_offset) {
//   size_t block_size = 0; 
  // size_t batch_size = edges.size()/num_batch;
  size_t num_batch =  (edges.size() + batch_size - 1) / batch_size;
  std::cout << "batch_size " << batch_size << std::endl;
  size_t count = 0;
  DBTGraph::DyGraph<DBTGraph::SymGraph> *DGnew;
  for(size_t i = batch_offset; i <= num_batch; ++i){
    size_t batch_end = min(batch_size  * (i+1), edges.size());
    if(batch_end == batch_size * i)  break;
    timer t; t.start();
    tie(count, DGnew)= DBTGraph::majorRebalancing(edges, 0,  batch_end, n, 0, P, false);
    std::cout << "batch " << i << " [" << 0 << " " << batch_end << "]" << std::endl;
    t.stop();t.reportTotal("");
    PrintBreak();
  }
  return count;
}


  //edges already randomly shuffled
  template <class UT>
  inline size_t staticCountMixed(const std::vector<UT>& edges, size_t batch_size, commandLine& P, size_t n, size_t batch_offset) {
//   size_t block_size = 0; 
  // size_t batch_size = edges.size()/num_batch;

  std::vector<UT> edges_new;

  size_t num_batch =  edges.size() / batch_size;
  size_t end = num_batch * batch_size;
  for (size_t i = end/2; i< end; ++i) {
    edges_new.push_back(edges[i]);
  }
  for (size_t i = 0; i< end/2; ++i) {
    edges_new.push_back(edges[i]);
  }
  std::cout << "batch_size " << batch_size << std::endl;
  size_t count = 0;
  DBTGraph::DyGraph<DBTGraph::SymGraph> *DGnew;
  for(size_t i = batch_offset; i <= num_batch; ++i){
    size_t batch_end = i* batch_size/2 + end/2;
    // if(batch_end == batch_size * i)  break;
    timer t; t.start();
    tie(count, DGnew)= DBTGraph::majorRebalancing(edges_new, i* batch_size/2,  batch_end, n, 0, P, false);
    std::cout << "batch " << i << " [" <<  i* batch_size/2 << " " << batch_end << "]" << std::endl;
    t.stop();t.reportTotal("");
    PrintBreak();
  }
  edges_new.clear();
  return count;
}

} // namespace DBTInternal

} // namespace gbbs
#endif
