#ifndef SHARED_H
#define SHARED_H

// *************************************************************
//   Misc
// *************************************************************
#include <string>
#include <iostream>
#include <tuple>
#include "gbbs/gbbs.h"

namespace gbbs{

namespace DBTInternal {

  template <class EdgeT>
  pbbs::sequence<pair<EdgeT, bool>> generateEdgeUpdates(size_t nv, size_t numEdges, commandLine& P);

  ///// PRINTING UTILITIES //////
  inline void PrintCaption(std::string t_in);

  inline void PrintSubcaption(std::string t_in);

  template <class T>
  inline void PrintFunctionTitle(std::string t_func, T t_suffix);

  template <class T>
  inline void PrintFunctionItem(std::string t_func, std::string t_item, T t_suffix);

  template <class T>
  inline void PrintSubtimer(std::string t_item, T t_suffix);

  inline void PrintBreak();

  template <typename A>
  inline void PrintVec(A *t_vec, int t_len);

  template <typename A>
  inline void PrintVec(A t_vec, int t_len);

  template <typename A>
  inline void PrintPairSeq(A t_vec);

  template <typename A>
  inline void PrintVec2(A *t_vec, int t_len);

template <class UT>
inline void staticCount(std::vector<UT>& edges, size_t num_batch);

} // namespace DBTInternal

} // namespace gbbs
#endif
