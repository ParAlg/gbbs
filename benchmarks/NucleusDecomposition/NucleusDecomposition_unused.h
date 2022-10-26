#pragma once

#include <math.h>
#include <limits>

// Library dependencies
#include "gbbs/bucket.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"
#include "gbbs/helpers/dyn_arr.h"
#include "gbbs/helpers/sparse_table.h"
#include "gbbs/helpers/sparse_additive_map.h"

// Ordering files
#include "benchmarks/DegeneracyOrder/BarenboimElkin08/DegeneracyOrder.h"
#include "benchmarks/DegeneracyOrder/GoodrichPszona11/DegeneracyOrder.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"
#include "benchmarks/CliqueCounting/Clique.h"

// Clique files
#include "benchmarks/CliqueCounting/intersect.h"
#include "benchmarks/CliqueCounting/induced_intersection.h"
#include "benchmarks/CliqueCounting/induced_neighborhood.h"
#include "benchmarks/CliqueCounting/induced_hybrid.h"
#include "benchmarks/CliqueCounting/induced_split.h"
#include "benchmarks/CliqueCounting/relabel.h"

#include "twotable.h"
#include "twotable_nosearch.h"
#include "onetable.h"
#include "commontable.h"
#include "list_buffer.h"

#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"

namespace gbbs {


  template<class Obj>
  class ThreadLocalObj {
    public:
      sequence<unsigned int> table_mark;
      sequence<Obj*> table_obj;
      int nw;
      ThreadLocalObj(){
        nw = num_workers() * 1.5;
        table_mark = sequence<unsigned int>::from_function(nw, [](std::size_t i) {return 0;});
        table_obj = sequence<Obj*>::from_function(nw, [](std::size_t i){return nullptr;});
      }

      Obj* init_idx(unsigned int idx) {
        auto obj = table_obj[idx];
        if (obj != nullptr) return obj;
        table_obj[idx] = new Obj();
        return table_obj[idx];
      }

      void unreserve(unsigned int idx) {
        while(true) {
         if (gbbs::atomic_compare_and_swap(&(table_mark[idx]), (unsigned int) 1, (unsigned int) 0)) { return; }
        }
      }

      std::pair<unsigned int, Obj*> reserve(){
        auto wid = worker_id();
        auto idx = parlay::hash32(wid) % nw;
        while(true) {
          if (table_mark[idx] == 0) {
            if (gbbs::atomic_compare_and_swap(&(table_mark[idx]), (unsigned int) 0, (unsigned int) 1)) {
              return std::make_pair(idx, init_idx(idx));
            }
          }
          idx++;
          idx = idx % nw;
        }
      }
  };


} // namespace gbbs