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
#include "gbbs/pbbslib/dyn_arr.h"
#include "gbbs/pbbslib/sparse_table.h"
#include "gbbs/pbbslib/sparse_additive_map.h"
#include "pbbslib/assert.h"
#include "pbbslib/list_allocator.h"
#include "pbbslib/integer_sort.h"

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

#include "multitable.h"
#include "twotable.h"
#include "twotable_nosearch.h"
#include "onetable.h"
#include "commontable.h"
#include "multitable_nosearch.h"
#include "NucleusDecomposition_common.h"

// list buffer is to collate indices of r-cliques with changed s-clique counts while peeling
// we use per processor arrays to collate updated clique counts for r-cliques
// this helps with contention, without having to allocate large hash tables

namespace gbbs {

template <class Graph, class DirectedGraph, class Table, class Table2>
inline sequence<size_t> NucleusDecompositionVerificationRunner(Graph& GA, DirectedGraph& DG,
  size_t r, size_t s, Table& table, Table2& table2,
  size_t max_deg, sequence<uintE>& rank) {

  //std::cout << "Start count" << std::endl;
  timer t; t.start();
  size_t count = CountCliquesNuc(DG, s, r, max_deg, &table);
  size_t count2 = CountCliquesNuc(DG, s, r, max_deg, &table2);
  assert(count == count2);
  double tt = t.stop();
  //std::cout << "End count" << std::endl;

  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << s << " cliques = " << count << "\n";

  timer t2; t2.start();
  auto peel = Peel_verify<std::size_t>(GA, DG, r, s, &table, &table2, rank);
  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  return peel;
}

template <class Graph, class DirectedGraph, class Table>
inline sequence<size_t> NucleusDecompositionRunner(Graph& GA, DirectedGraph& DG,
  size_t r, size_t s, Table& table, 
  size_t max_deg, sequence<uintE>& rank, size_t efficient, bool relabel) {

  //std::cout << "Start count" << std::endl;
  timer t; t.start();
  size_t count = CountCliquesNuc(DG, s, r, max_deg, &table);
  double tt = t.stop();
  //std::cout << "End count" << std::endl;

  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << s << " cliques = " << count << "\n";

  timer t2; t2.start();
  auto peel = Peel<std::size_t>(GA, DG, r, s, &table, rank, efficient, relabel);
  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  return peel;
}

template<class T>
T round_up(T dividend, T divisor)
{
    return (dividend + (divisor - 1)) / divisor;
}

template <class T, class H, class Graph, class Graph2>
inline sequence<size_t> runner_verify(Graph& GA, Graph2& DG, size_t r, size_t s, long table_type, long num_levels,
  bool relabel, bool contiguous_space, size_t max_deg, sequence<uintE>& rank, int shift_factor) {
  nd_global_shift_factor = shift_factor;
  onetable::OnelevelHash<T, H> table(r, DG, max_deg, shift_factor);
  twotable_nosearch::TwolevelHash<T, H> table2(r, DG, max_deg, relabel, shift_factor);
  return NucleusDecompositionVerificationRunner(GA, DG, r, s, table, table2, max_deg, rank);
}

template <class T, class H, class Graph, class Graph2>
inline sequence<size_t> runner(Graph& GA, Graph2& DG, size_t r, size_t s, long table_type, long num_levels,
  bool relabel, bool contiguous_space, size_t max_deg, sequence<uintE>& rank, int shift_factor,
  size_t efficient = 1) {
  timer t; 
  sequence<size_t> count;
  nd_global_shift_factor = shift_factor;

  if (table_type == 3) {
    t.start();
    // Num levels matches, e.g., 2 for two level
    num_levels -= 1;
    if (!relabel) {
      auto rank_func = [&](uintE a, uintE b){ return rank[a] < rank[b]; };
      multitable::MHash<T, H, decltype(rank_func)> table(r, DG, max_deg, num_levels, contiguous_space, rank_func);
      double tt = t.stop();
      std::cout << "### Table Running Time: " << tt << std::endl;
      count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank, efficient, relabel);
    } else {
      auto rank_func = std::less<uintE>();
      multitable::MHash<T, H, decltype(rank_func)> table(r, DG, max_deg, num_levels, contiguous_space, rank_func);
      double tt = t.stop();
      std::cout << "### Table Running Time: " << tt << std::endl;
      count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank, efficient, relabel);
    }
  } else if (table_type == 2) {
    t.start();
    twotable::TwolevelHash<T, H> table(r, DG, max_deg, contiguous_space, relabel, shift_factor);
    double tt = t.stop();
    std::cout << "### Table Running Time: " << tt << std::endl;
    count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank, efficient, relabel);
  } else if (table_type == 1) {
    t.start();
    onetable::OnelevelHash<T, H> table(r, DG, max_deg, shift_factor);
    double tt = t.stop();
    std::cout << "### Table Running Time: " << tt << std::endl;
    count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank, efficient, relabel);
  } else if (table_type == 4) {
    // Num levels matches, e.g., 2 for two level
    num_levels -= 1;
    if (!relabel) {
      auto rank_func = [&](uintE a, uintE b){ return rank[a] < rank[b]; };
      multitable_nosearch::MHash<T, H, decltype(rank_func)> table(r, DG, max_deg, num_levels, rank_func);
      double tt = t.stop();
      std::cout << "### Table Running Time: " << tt << std::endl;
      count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank, efficient, relabel);
    } else {
      auto rank_func = std::less<uintE>();
      multitable_nosearch::MHash<T, H, decltype(rank_func)> table(r, DG, max_deg, num_levels, rank_func);
      double tt = t.stop();
      std::cout << "### Table Running Time: " << tt << std::endl;
      count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank, efficient, relabel);
    }
  } else if (table_type == 5) {
    t.start();
    twotable_nosearch::TwolevelHash<T, H> table(r, DG, max_deg, relabel, shift_factor);
    double tt = t.stop();
    std::cout << "### Table Running Time: " << tt << std::endl;
    count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank, efficient, relabel);
  } 
  return count;
}

template <class Graph>
inline sequence<size_t> NucleusDecomposition(Graph& GA, size_t r, size_t s, long table_type, long num_levels,
  bool relabel, bool contiguous_space, bool verify, size_t efficient) {
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

  int num_bits_in_n = 1 + pbbslib::log2_up(DG.n + 1); //32
  int num_bytes_needed = round_up<int>(((std::max(static_cast<int>(1), 
    static_cast<int>(r - (num_levels - 1))) * num_bits_in_n) + 1), 8);
  int shift_factor = num_bits_in_n; //32

  std::cout << "Num bytes needed: " << num_bytes_needed << std::endl;
  std::cout << "Num bits in n: " << shift_factor << std::endl;
  fflush(stdout);

  if (num_bytes_needed <= 4 && table_type != 5 && table_type != 4) {
    // unsigned __int32
    count = runner<unsigned int, nhash32>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
      max_deg, rank, shift_factor);
  } else if (num_bytes_needed <= 8) {
    // unsigned __int64
    count = runner<unsigned long long, nhash64>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
      max_deg, rank, shift_factor);
  } else {
    // unsigned__int128
    if (!verify)
      count = runner<unsigned __int128, hash128>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
        max_deg, rank, shift_factor, efficient);
    else
      count = runner_verify<unsigned __int128, hash128>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
        max_deg, rank, shift_factor);
  }

   

  //table.del();
  DG.del();

  return count;
}

}