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
namespace UTIL {

  // struct UpdateTag{
  //   bool is_insert;
  //   double t; // timestamp
  //   UpdateTag(bool a, double b):is_insert(a), t(b){}
  // };

  template <class EdgeT>
  pbbs::sequence<pair<EdgeT, bool>> generateEdgeUpdates(size_t nv, size_t numEdges){
    return pbbs::sequence<pair<EdgeT, bool>>(numEdges, [&](size_t i) { 
      uintE u,v = 0;
      do{u = rand()%nv;v = rand()%nv;}while (u == v);

      return make_pair(EdgeT(u,v), true);
     });

  }

  ///// PRINTING UTILITIES //////
  inline void PrintCaption(std::string t_in) {
    std::cout << "========= " << t_in << " =========" << std::endl;
  }

  inline void PrintSubcaption(std::string t_in) {
    std::cout << "========= " << t_in << std::endl;
  }

  template <class T>
  inline void PrintFunctionTitle(std::string t_func, T t_suffix) {
    std::cout << "[" << t_func << "] " << t_suffix << std::endl;
  }

  template <class T>
  inline void PrintFunctionItem(std::string t_func, std::string t_item, T t_suffix) {
    std::cout << "[" << t_func << "] " << t_item << " = " << t_suffix << std::endl;
  }

  template <class T>
  inline void PrintSubtimer(std::string t_item, T t_suffix) {
    std::cout << "::" << t_item << ": " << t_suffix << std::endl;
  }

  inline void PrintBreak() {
    std::cout << std::endl;
  }

  template <typename A>
  inline void PrintVec(A *t_vec, int t_len) {
      for (int i = 0; i < t_len; ++ i) {
	      std::cout << t_vec[i] << ' ';
      }
      std::cout << std::endl;
    }

  template <typename A>
  inline void PrintVec(A t_vec, int t_len) {
      for (int i = 0; i < t_len; ++ i) {
	      std::cout << t_vec[i] << ' ';
      }
      std::cout << std::endl;
    }

  template <typename A>
  inline void PrintPairSeq(A t_vec) {
      for (size_t i = 0; i < t_vec.size(); ++ i) {
	      std::cout << "(" << t_vec[i].first << "," << t_vec[i].second << ") ";
      }
      std::cout << std::endl;
  }

  template <typename A>
  inline void PrintVec2(A *t_vec, int t_len) {
      for (int i = 0; i < t_len; ++ i) {
	     t_vec[i].print();
      }
      std::cout << std::endl;
    }

}
}
#endif
