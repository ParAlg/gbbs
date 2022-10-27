// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

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
#include "NucleusDecomposition_common.h"

// list buffer is to collate indices of r-cliques with changed s-clique counts while peeling
// we use per processor arrays to collate updated clique counts for r-cliques
// this helps with contention, without having to allocate large hash tables

namespace gbbs {

  // For approximate nucleus decomp take log of core value (log_{1+eps})

inline std::vector<uintE> CompressConnect(std::vector<uintE>& connect1, size_t num) {
  std::cout << "Begin compress" << std::endl; fflush(stdout);
  std::vector<uintE> compress(connect1.size(), UINT_E_MAX);
  std::vector<uintE> current_parent(num);
  parallel_for(0, num, [&](size_t i){current_parent[i] = i;});

  auto empty = std::make_tuple<uintE, uintE>(UINT_E_MAX, 0);
  auto duplicate_table = gbbs::sparse_additive_map<uintE, uintE>(num, empty);
  parallel_for(0, connect1.size(), [&](size_t i){
    if (connect1[i] != UINT_E_MAX) duplicate_table.insert(std::make_tuple(connect1[i], 1));
  });
  /*for (size_t i = 0; i < num; i++) {
    if (connect1[i] == UINT_E_MAX) continue;
    bool has_duplicate = false;
    for (size_t j = 0; j < num; j++) {
      if (i == j) continue;
      if (connect1[i] == connect1[j]) {
        has_duplicate = true;
        break;
      }
    }
    if (has_duplicate || connect1[i] >= connect1.size()) {
      compress[i] = connect1[i];
    } else {
      compress[i] = connect1[connect1[i]];
    }
    current_parent[i] = compress[i];
  }*/
  bool has_changed = true;
  while (has_changed) {
    has_changed = false;
    for (size_t i = 0; i < num; i++) {
      auto idx = current_parent[i];
      if (idx == UINT_E_MAX || idx >= connect1.size() || connect1[idx] == UINT_E_MAX) continue;
      if (compress[idx] == UINT_E_MAX) has_changed = true;

      bool has_duplicate = false;
      if (duplicate_table.find(idx) > 1) has_duplicate = true;
      /*for (size_t j = 0; j < connect1.size(); j++) {
        if (idx == j) continue;
        if (connect1[j] == connect1[idx]) {
          has_duplicate = true;
          break;
        }
      }*/

      if (has_duplicate) {
        compress[idx] = connect1[idx];
      } else if (connect1[idx] >= connect1.size()){
        //compress[idx] = connect1[idx];
      } else {
        compress[idx] = connect1[connect1[idx]];
      }
      current_parent[i] = compress[idx];
    }
  }
  std::cout << "End compress" << std::endl; fflush(stdout);
  return compress;
}

inline void CheckConnect(std::vector<uintE>& connect1, std::vector<uintE>& connect2, size_t num) {
  std::vector<uintE> changed1;
  uintE max_val = 0;
  for (size_t i = 0; i < connect1.size(); i++) {
    if (connect1[i] != UINT_E_MAX && connect1[i] > max_val) max_val = connect1[i];
  }
  for (size_t i = 0; i < connect2.size(); i++) {
    if (connect2[i] != UINT_E_MAX && connect2[i] > max_val) max_val = connect2[i];
  }
  max_val++;
  std::vector<uintE> map(max_val);
  std::cout << "start check" << std::endl;
  parallel_for(0, map.size(), [&](size_t i){map[i] = UINT_E_MAX;});
  for (size_t i = 0; i < num; i++) {
    if (connect1[i] == UINT_E_MAX ) {
      assert(connect2[i] == UINT_E_MAX );
      if (connect2[i] != UINT_E_MAX ) {
        std::cout << "Err first round: connect2 i " << i << ", " << connect2[i] << std::endl;
        exit(0);
      }
    } else if (map[connect1[i]] == UINT_E_MAX) {
      map[connect1[i]] = connect2[i];
      changed1.push_back(connect1[i]);
    } else {
      assert(map[connect1[i]] == connect2[i]);
      if (map[connect1[i]] != connect2[i]) {
        std::cout << "Why are these different: map: " << map[connect1[i]] << ", connect2: " << connect2[i] << std::endl;
        exit(0);
      }
    }
  }
  //std::cout << "finish round 1" << std::endl;
  int x = 2;
  while(!changed1.empty()) {
    if (changed1[0] >= connect1.size()) break;
    std::vector<uintE> newchanged;
    for (size_t i = 0; i < changed1.size(); i++) {
      size_t i1 = changed1[i];
      size_t i2 = map[i1];
      if (connect1[i1] == UINT_E_MAX ) {
        assert(connect2[i2] == UINT_E_MAX);
        if (connect2[i2] != UINT_E_MAX) {
          std::cout << "Err second round: connect2 i " << i2 << ", " << connect2[i2] << std::endl;
          exit(0);
        }
      } else if (map[connect1[i1]] == UINT_E_MAX) {
        map[connect1[i1]] = connect2[i2];
        newchanged.push_back(connect1[i1]);
      } else {
        assert(map[connect1[i1]] == connect2[i2]);
        if (map[connect1[i1]] != connect2[i2]) {
          std::cout << "Why are these different second: map: " << map[connect1[i1]] << ", connect2: " << connect2[i2] << std::endl;
          exit(0);
        }
      }
    }
    changed1 = newchanged;
    //std::cout << "finish round " << x << std::endl;
    x++;
  }
}

template <class iden_t, class bucket_t, class Graph, class DirectedGraph, class Table>
inline sequence<bucket_t> NucleusDecompositionRunner(Graph& GA, DirectedGraph& DG,
  size_t r, size_t s, Table& table, 
  size_t max_deg, sequence<uintE>& rank, size_t efficient, bool relabel,
  bool use_compress, bool inline_hierarchy, bool efficient_inline_hierarchy, bool verify) {
  if (efficient_inline_hierarchy) inline_hierarchy = true;

  //std::cout << "Start count" << std::endl;
  timer t; t.start();
  size_t count;
  // efficient = 3 is fake; it means efficient = 1, but run PND clique counting
  if (efficient == 3) {
    count = CountCliquesNucPND(DG, s, r, max_deg, &table);
  }
  else {
    count = CountCliquesNuc(DG, s, r, max_deg, &table);
  }
  double tt = t.stop();
  //std::cout << "End count" << std::endl;

  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << s << " cliques = " << count << "\n";

  timer t2; t2.start();
  sequence<bucket_t> peel;
  EfficientConnectWhilePeeling ecwp;
  ConnectWhilePeeling connect_with_peeling;
  auto num_entries = table.return_total();
  if (inline_hierarchy && !efficient_inline_hierarchy) connect_with_peeling = ConnectWhilePeeling(num_entries);
  else if (efficient_inline_hierarchy) ecwp = EfficientConnectWhilePeeling(num_entries);
  if (use_compress) {
    if (!efficient_inline_hierarchy) peel = Peel_space_efficient<bucket_t, iden_t>(GA, DG, r, s, &table, rank, efficient, relabel, use_compress, inline_hierarchy, connect_with_peeling);
    else peel = Peel_space_efficient<bucket_t, iden_t>(GA, DG, r, s, &table, rank, efficient, relabel, use_compress, inline_hierarchy, ecwp);
  } else {
    if (!efficient_inline_hierarchy) peel = Peel<bucket_t>(GA, DG, r, s, &table, rank, efficient, relabel, use_compress, inline_hierarchy, connect_with_peeling);
    else peel = Peel<bucket_t>(GA, DG, r, s, &table, rank, efficient, relabel, use_compress, inline_hierarchy, ecwp);
  }
  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  std::vector<uintE> connect;
  if (!inline_hierarchy) {
    std::cout << "Running Connectivity" << std::endl;
    timer t3; t3.start();
    connect = construct_nd_connectivity(peel, GA, DG, r-1, s-1, table, rank, relabel);
    double tt3 = t3.stop();
    std::cout << "### Connectivity Running Time: " << tt3 << std::endl;
  } else {
    std::cout << "Constructing tree" << std::endl;
    timer t3; t3.start();
    if (!efficient_inline_hierarchy) connect = construct_nd_connectivity_from_connect(connect_with_peeling, peel, GA, DG, r-1, s-1, table, rank, relabel);
    else connect = construct_nd_connectivity_from_connect(ecwp, peel, GA, DG, r-1, s-1, table, rank, relabel);
    double tt3 = t3.stop();
    std::cout << "### Connectivity Tree Running Time: " << tt3 << std::endl;
    if (verify) {
      std::cout << "Running Connectivity Verify" << std::endl;
      t3.start();
      auto connect2 = construct_nd_connectivity(peel, GA, DG, r-1, s-1, table, rank, relabel);
      tt3 = t3.stop();
      std::cout << "### Connectivity Running Time: " << tt3 << std::endl;
    /*std::cout << "Printing tree 1: " << std::endl;
    for (std::size_t i = 0; i < connect.size(); i++) {
      std::cout << i << ": " << connect[i] << std::endl;
    }*/
    
      if (efficient_inline_hierarchy) connect2 = CompressConnect(connect2, table.return_total());

     /* std::cout << "Printing tree 2: " << std::endl;
    for (std::size_t i = 0; i < connect2.size(); i++) {
      std::cout << i << ": " << connect2[i] << std::endl;
    }*/

      CheckConnect(connect, connect2, table.return_total());
      CheckConnect(connect2, connect, table.return_total());
    }
  }

  /*std::cout << "Printing tree: " << std::endl;
  for (std::size_t i = 0; i < connect.size(); i++) {
    std::cout << i << ": " << connect[i] << std::endl;
  }*/

  return peel;
}

template<class T>
T round_up(T dividend, T divisor)
{
    return (dividend + (divisor - 1)) / divisor;
}

template <class bucket_t, class T, class H, class Graph, class Graph2>
inline sequence<bucket_t> runner(Graph& GA, Graph2& DG, size_t r, size_t s, long table_type, long num_levels,
  bool relabel, bool contiguous_space, size_t max_deg, sequence<uintE>& rank, int shift_factor,
  size_t efficient, bool use_compress, bool output_size, bool inline_hierarchy, bool efficient_inline_hierarchy,
  bool verify) {
  timer t; 
  //sequence<size_t> count;
  nd_global_shift_factor = shift_factor;

  //if (table_type == 2) {
  //  t.start();
  //  twotable::TwolevelHash<T, H, bucket_t> table(r, DG, max_deg, contiguous_space, relabel, shift_factor);
  //  double tt = t.stop();
  //  std::cout << "### Table Running Time: " << tt << std::endl;
  //  return NucleusDecompositionRunner<T, bucket_t>(GA, DG, r, s, table, max_deg, rank, efficient, relabel, use_compress, inline_hierarchy);
  //} else 
  /*if (table_type == 1) {
    t.start();
    onetable::OnelevelHash<T, H, bucket_t> table(r, DG, max_deg, shift_factor);
    double tt = t.stop();
    std::cout << "### Table Running Time: " << tt << std::endl;
    return NucleusDecompositionRunner<T, bucket_t>(GA, DG, r, s, table, max_deg, rank, efficient, relabel, use_compress, inline_hierarchy);
  } */ //else if (table_type == 5) {
    t.start();
    twotable_nosearch::TwolevelHash<T, H, bucket_t> table(r, DG, max_deg, relabel, shift_factor);
    double tt = t.stop();
    std::cout << "### Table Running Time: " << tt << std::endl;
    return NucleusDecompositionRunner<T, bucket_t>(GA, DG, r, s, table, max_deg, rank, efficient, relabel, use_compress, inline_hierarchy, efficient_inline_hierarchy, verify);
  //} 
  //return count;
}

template <class Graph>
inline void NucleusDecomposition(Graph& GA, size_t r, size_t s, long table_type, long num_levels,
  bool relabel, bool contiguous_space, bool verify, size_t efficient, bool use_compress,
  bool output_size, bool inline_hierarchy, bool efficient_inline_hierarchy) {
  // TODO: if r = 2
  using W = typename Graph::weight_type;

  // Obtain vertex ordering
  timer t_rank; t_rank.start();
  sequence<uintE> rank = get_ordering(GA, 3, 0.1); // in clique counting
  double tt_rank = t_rank.stop();
  std::cout << "### Rank Running Time: " << tt_rank << std::endl;

  // Direct the graph based on ordering
  timer t_filter; t_filter.start();
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return (rank[u] < rank[v]);// && GA.get_vertex(u).getOutDegree() >= r-1 && GA.get_vertex(v).getOutDegree() >= r-1;
  };
  // Note: If relabeling, core #s must be translated back
  auto DG = relabel ? relabel_graph(GA, rank.begin(), pack_predicate) : filterGraph(GA, pack_predicate);
  double tt_filter = t_filter.stop();
  std::cout << "### Filter Graph Running Time: " << tt_filter << std::endl;

  auto max_deg = get_max_deg3(DG);


  sequence<size_t> count;

  // unsigned __int128 is 16 bytes
  // unsigned __int32 is 4 bytes
  // unsigned __int64 is 8 bytes

  // Let X be the number of bits needed to express max vertex
  // We have (r, s)
  // round_up<int>(((max(1, (r - (num_levels - 1))) * X) + 1), 8)
  // use whichever type is >= bytes than this

  if (table_type == 1) num_levels = 1;
  else if (table_type == 2 || table_type == 5) num_levels = 2;

  int num_bits_in_n = 1 + parlay::log2_up(DG.n + 1); //32
  int num_bytes_needed = round_up<int>(((std::max(static_cast<int>(1), 
    static_cast<int>(r - (num_levels - 1))) * num_bits_in_n) + 1), 8);
  int shift_factor = num_bits_in_n; //32

  std::cout << "Num bytes needed: " << num_bytes_needed << std::endl;
  std::cout << "Num bits in n: " << shift_factor << std::endl;
  fflush(stdout);

  if (r == 2 && s == 3) {
    using bucket_t = uintE;

  if (num_bytes_needed <= 4 && table_type != 5 && table_type != 4) {
    // unsigned __int32
    auto peel = runner<bucket_t, unsigned int, nhash32>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
      max_deg, rank, shift_factor, efficient, use_compress, output_size, inline_hierarchy, efficient_inline_hierarchy, verify);
  } else if (num_bytes_needed <= 8) {
    // unsigned __int64
    auto peel = runner<bucket_t, unsigned long long, nhash64>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
      max_deg, rank, shift_factor, efficient, use_compress, output_size, inline_hierarchy, efficient_inline_hierarchy, verify);
  } else {
    // unsigned__int128
      auto peel = runner<bucket_t, unsigned __int128, hash128>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
        max_deg, rank, shift_factor, efficient, use_compress, output_size, inline_hierarchy, efficient_inline_hierarchy, verify);
  }


  } else {

  using bucket_t = size_t;

  if (num_bytes_needed <= 4 && table_type != 5 && table_type != 4) {
    // unsigned __int32
    auto peel = runner<bucket_t, unsigned int, nhash32>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
      max_deg, rank, shift_factor, efficient, use_compress, output_size, inline_hierarchy, efficient_inline_hierarchy, verify);
  } else if (num_bytes_needed <= 8) {
    // unsigned __int64
    auto peel = runner<bucket_t, unsigned long long, nhash64>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
      max_deg, rank, shift_factor, efficient, use_compress, output_size, inline_hierarchy, efficient_inline_hierarchy, verify);
  } else {
    // unsigned__int128
    auto peel = runner<bucket_t, unsigned __int128, hash128>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
        max_deg, rank, shift_factor, efficient, use_compress, output_size, inline_hierarchy, efficient_inline_hierarchy, verify);
  }

  }

   

  //table.del();
  //DG.del();

  //return count;
}

}