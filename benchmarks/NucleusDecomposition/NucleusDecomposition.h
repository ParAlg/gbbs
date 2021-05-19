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

namespace gbbs {

  template <class Graph, class T>
  inline size_t CountCliquesNuc(Graph& DG, size_t k, size_t r, size_t max_deg, T* table) {
    k--; r--;
    timer t2; t2.start();

    auto add_f = [&] (long* ct, const std::tuple<unsigned __int128, long>& tup) {
      pbbs::fetch_and_add(ct, (long)1);
    };
// Nested hash tables where the first level you hash the first vertex, second level you hash the second, etc.
// Saves space if cliques share vertices
// Renumber vert by which is in the most r-cliques
    auto base_f = [&](sequence<uintE>& base){
      table->insert(base, r, k);
    };
    auto tots = sequence<size_t>(DG.n, size_t{0});

    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, DG.n, true, true); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del();
    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
        if (DG.get_vertex(i).getOutDegree() != 0) {
          induced->setup(DG, k, i);
          auto base = sequence<uintE>(k + 1);
          base[0] = i;
          //auto base_f2 = [&](uintE vtx, size_t _count) {};
          //tots[i] = induced_hybrid::KCliqueDir_fast_hybrid_rec(DG, 1, k, induced, base_f2, 0);
          tots[i] = NKCliqueDir_fast_hybrid_rec(DG, 1, k, induced, base_f, base);
        } else tots[i] = 0;
    }, 1, false);
    double tt2 = t2.stop();
    std::cout << "##### Actual counting: " << tt2 << std::endl;

    return pbbslib::reduce_add(tots);
  }

class list_buffer {
  public:
    int buffer;
    sequence<size_t> list;
    sequence<size_t> starts;
    sequence<bool> to_pack;
    size_t next;
    size_t num_workers2;
    list_buffer(size_t s){
      num_workers2 = num_workers();
      buffer = 64;
      int buffer2 = 64;
      list = sequence<size_t>(s + buffer2 * num_workers2, static_cast<size_t>(0));
      starts = sequence<size_t>(num_workers2, [&](size_t i){return i * buffer2;});
      next = num_workers2 * buffer2;
      to_pack = sequence<bool>(s + buffer2 * num_workers2, true);
    }
    void add(size_t index) {
      size_t worker = worker_id();
      list[starts[worker]] = index;
      starts[worker]++;
      if (starts[worker] % buffer == 0) {
        size_t use_next = pbbs::fetch_and_add(&next, buffer);
        starts[worker] = use_next;
      }
    }

    template <class I>
    size_t filter(I update_changed, sequence<size_t>& per_processor_counts) {
      parallel_for(0, num_workers2, [&](size_t worker) {
        size_t divide = starts[worker] / buffer;
        for (size_t j = starts[worker]; j < (divide + 1) * buffer; j++) {
          to_pack[j] = false;
        }
      });
      // Pack out 0 to next of list into pack
      parallel_for(0, next, [&] (size_t i) {
        if (to_pack[i])
          update_changed(per_processor_counts, i, list[i]);
        else
          update_changed(per_processor_counts, i, UINT_E_MAX);
      });
      parallel_for(0, num_workers2, [&](size_t worker) {
        size_t divide = starts[worker] / buffer;
        for (size_t j = starts[worker]; j < (divide + 1) * buffer; j++) {
          to_pack[j] = true;
        }
      });
      return next;
    }

    void reset() {
      parallel_for (0, num_workers2, [&] (size_t j) {
        starts[j] = j * buffer;
      });
      next = num_workers2 * buffer;
    }
};

template <class Graph, class Graph2, class F, class H, class I, class T>
inline size_t cliqueUpdate(Graph& G, Graph2& DG, size_t r, 
size_t k, size_t max_deg, bool label, F get_active, size_t active_size,
  size_t granularity, char* still_active, sequence<uintE> &rank, 
  sequence<size_t>& per_processor_counts, H update,
  bool do_update_changed, I update_changed,
  T* cliques, size_t n, list_buffer& count_idxs, timer& t1) {

  // Set up space for clique counting
  auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k-r+1, G.n, true, true); };
  auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };

  // Mark every vertex in the active set
  parallel_for (0, active_size, [&] (size_t j) {
    auto index = get_active(j); //cliques->find_index(get_active(j));
    still_active[index] = 1;
    }, 2048);

  // Hash table to contain clique count updates
  //size_t edge_table_size = (size_t) (n);
  //auto edge_table = pbbslib::sparse_table<uintE, bool, hashtup>(edge_table_size, std::make_tuple(0, false), hashtup());

  // Function that dictates which edges to consider in first level of recursion
  //auto ignore_f = [&](const uintE& u, const uintE& v) {
  //  return true;
    //auto index_u = cliques->find_index(u);
    //auto index_v = cliques->find_index(v);
    //auto status_u = still_active[index_u]; auto status_v = still_active[index_v];
    //if (status_u == 2 || status_v == 2) return false; // deleted edge
    //if (status_v == 0) return true; // non-deleted, non-active edge
    //return rank[u] < rank[v]; // orient edges if in active set
  //};

  // Collate clique counts by processor
  //count_idxs[0] = 0;
  auto update_d = [&](sequence<uintE>& base){
    cliques->extract_indices(base, [&](std::size_t index){
      size_t ct = pbbs::fetch_and_add(&(per_processor_counts[index]), 1);
      if (ct == 0) count_idxs.add(index);
    }, r, k);
  };

t1.start();
  // Clique count updates
  parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active_size,
                                     [&](size_t i, HybridSpace_lw* induced) {
    // TODO: THIS PART IS WRONG
    // you wanna start from the clique given by vert
    auto x = get_active(i);
    auto base = sequence<uintE>(k + 1);
    cliques->extract_clique(x, base, G, k);
      // Fill base[k] ... base[k-r+1] and base[0]
      induced->setup_nucleus(G, DG, k, base, r);
      // Need to fix so that k_idx is 1, but ends as if it was r
      NKCliqueDir_fast_hybrid_rec(G, 1, k-r, induced, update_d, base);
  }, granularity, false);
t1.stop();

  // Extract all vertices with changed clique counts
  //auto changed_vtxs = edge_table.entries();
  //edge_table.del();

  // Aggregate the updated counts across all worker's local arrays, as specified by update
  /*parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
    size_t nthreads = num_workers();
    uintE v = std::get<0>(changed_vtxs[i]);
    auto index = cliques->find_index(v);
    for (size_t j=0; j<nthreads; j++) {
      update(per_processor_counts, j, index);
    }
  }, 128);*/

  // Perform update_changed on each vertex with changed clique counts
  std::size_t num_count_idxs = 0;
  if (do_update_changed) {
    /*parallel_for(0, changed_vtxs.size(), [&] (size_t i) {
      auto index = cliques->find_index(std::get<0>(changed_vtxs[i]));
      update_changed(per_processor_counts, i, index);
    });*/
      num_count_idxs = count_idxs.filter(update_changed, per_processor_counts);
      count_idxs.reset();
    /*
    parallel_for(0, num_count_idxs, [&] (size_t i) {//count_idxs[0]
      //assert(count_idxs[i+1] < n);
      //assert(per_processor_counts[count_idxs[i+1]] > 0);
      update_changed(per_processor_counts, i, count_idxs.pack[i]);//count_idxs[i + 1]
    });*/
    
  }

  // Mark every vertex in the active set as deleted
  parallel_for (0, active_size, [&] (size_t j) {
    auto index = get_active(j); //cliques->find_index(get_active(j));
    still_active[index] = 2;}, 2048);

  return num_count_idxs; //count_idxs[0];
}

template <typename bucket_t, class Graph, class Graph2, class T>
sequence<bucket_t> Peel(Graph& G, Graph2& DG, size_t r, size_t k, 
  T* cliques, sequence<uintE> &rank,
  size_t num_buckets=16) {
    k--; r--;
  timer t2; t2.start();
// Here's the mistake: You're thinking that get_active should return the key,
// which is the concatenation of vertices. That's wrong. get_active should
// return an index into cliques->table, which we can then convert to the
// appropriate key. When we bucket, we only maintain indices into
// cliques->table, and we transfer to the key as necessary.

  //sequence<std::tuple<__int128 unsigned, long int> > entries = cliques->entries;
  size_t num_entries = cliques->return_total();
  auto D = sequence<bucket_t>(num_entries, [&](size_t i) -> bucket_t { 
    return cliques->get_count(i);
  });

  auto D_filter = sequence<std::tuple<uintE, bucket_t>>(num_entries);

  auto b = make_vertex_custom_buckets<bucket_t>(num_entries, D, increasing, num_buckets);

  auto per_processor_counts = sequence<size_t>(num_entries , static_cast<size_t>(0));
  
  list_buffer count_idxs(num_entries);
  //auto count_idxs = sequence<size_t>(num_entries + 1 + 
  //  count_buffer * num_workers(), static_cast<size_t>(0));

  char* still_active = (char*) calloc(num_entries, sizeof(char));
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active?
  //auto update_idxs = sequence<uintE>(max_deg);
  size_t n = num_entries;

  timer t_extract;
  timer t_count;
  timer t_update;
  timer t_x;

  size_t rounds = 0;
  size_t finished = 0;
  bucket_t cur_bkt = 0;
  bucket_t max_bkt = 0;
  double max_density = 0;
  bool use_max_density = false;
  // Peel each bucket
  auto update_clique = [&](sequence<size_t>& ppc, size_t j, uintE v) {
    if (j == 0) return;
    ppc[v] += ppc[j*n + v];
    ppc[j*n + v] = 0;
  };
  while (finished != n) {
    t_extract.start();
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    //auto active = vertexSubset(n, bkt.identifiers);
    auto active_size = (bkt.identifiers).size();
    cur_bkt = bkt.id;
    t_extract.stop();

    finished += active_size;
    max_bkt = std::max(cur_bkt, max_bkt);

    auto get_active = [&](size_t j) -> unsigned __int128 { return (bkt.identifiers)[j]; };
      //return active.vtx(j); };

    if (active_size == 0 || D[get_active(0)] == 0) continue;
    //std::cout << "PEEL" << std::endl;
    //fflush(stdout);

    // TESTING EXTRACTION OF KEy
    /*
    for (int i = 0; i < static_cast<int>(active_size); i++){
    auto x = get_active(i);
    auto key = std::get<0>((cliques->table)[x]);
       for (int j = 0; j < static_cast<int>(r)+1; ++j) {
          int extract = (int) key;
          if (static_cast<uintE>(extract) >= DG.n) {
            std::cout << "Extract: " << static_cast<uintE>(extract) << ", n: " << DG.n << std::endl;
          }
          assert(static_cast<uintE>(extract) < DG.n);
        key = key >> 32;
      }
    }*/

    size_t granularity = (cur_bkt * active_size < 10000) ? 1024 : 1;

    size_t filter_size = 0;

      auto update_changed = [&](sequence<size_t>& ppc, size_t i, uintE v){
        /* Update the clique count for v, and zero out first worker's count */
        //auto index = cliques->find_index(v);
        //auto val = std::get<1>((cliques->table)[index]) - ppc[v];
    //(cliques->table)[index] = std::make_tuple(std::get<0>((cliques->table)[index]), val);
        //cliques[v] -= ppc[v];
        if (v == UINT_E_MAX) {
          D_filter[i] = std::make_tuple(num_entries + 1, 0);
          return;
        }
        if (ppc[v] == 0) D_filter[i] = std::make_tuple(num_entries + 1, 0);
        else {
          auto val = cliques->update_count(v, ppc[v]);
        ppc[v] = 0;
        bucket_t deg = D[v];
        if (deg > cur_bkt) {
          bucket_t new_deg = std::max((bucket_t) val, (bucket_t) cur_bkt);
          D[v] = new_deg;
          // store (v, bkt) in an array now, pass it to apply_f below instead of what's there right now -- maybe just store it in D_filter?
          D_filter[i] = std::make_tuple(v, b.get_bucket(deg, new_deg));
        } else D_filter[i] = std::make_tuple(num_entries + 1, 0);
        }
      };
    t_count.start();
     filter_size = cliqueUpdate(G, DG, r, k, max_deg, true, get_active, active_size, 
     granularity, still_active, rank, per_processor_counts,
      update_clique, true, update_changed, cliques, num_entries, count_idxs, t_x);
      t_count.stop();

    auto apply_f = [&](size_t i) -> std::optional<std::tuple<unsigned __int128, bucket_t>> {
      auto v = std::get<0>(D_filter[i]);
      bucket_t bucket = std::get<1>(D_filter[i]);
      if (v != num_entries + 1) {
        if (still_active[v] != 2) return wrap(v, bucket);
      }
      return std::nullopt;
    };

t_update.start();
    b.update_buckets(apply_f, filter_size);

    parallel_for (0, active_size, [&] (size_t j) {
      auto index = get_active(j);
      //auto index = cliques->find_index(v);
      cliques->clear_count(index);
      //cliques[active.vtx(j)] = 0;
      }, 2048);
      t_update.stop();

    //active.del();

    rounds++;
  }

  t_extract.reportTotal("### Peel Extract time: ");
  t_count.reportTotal("### Peel Count time: ");
  t_update.reportTotal("### Peel Update time: ");
  t_x.reportTotal("Inner counting: ");

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  std::cout.precision(17);
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "clique core: " << static_cast<uintE>(max_bkt) << std::endl;
  if (use_max_density) std::cout << "max density: " << max_density << std::endl;

  b.del();
  free(still_active);

  return D;
}

template <class Graph, class DirectedGraph, class Table>
inline sequence<size_t> NucleusDecompositionRunner(Graph& GA, DirectedGraph& DG, size_t r, size_t s, Table& table, 
  size_t max_deg, sequence<uintE>& rank) {

  std::cout << "Start count" << std::endl;
  timer t; t.start();
  size_t count = CountCliquesNuc(DG, s, r, max_deg, &table);
  double tt = t.stop();
  std::cout << "End count" << std::endl;

  std::cout << "### Count Running Time: " << tt << std::endl;
  std::cout << "### Num " << s << " cliques = " << count << "\n";

  timer t2; t2.start();
  auto peel = Peel<std::size_t>(GA, DG, r, s, &table, rank);
  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  return peel;
}

template<class T>
T round_up(T dividend, T divisor)
{
    return (dividend + (divisor - 1)) / divisor;
}

template <class T, class Graph, class Graph2>
inline sequence<size_t> runner(Graph& GA, Graph2& DG, size_t r, size_t s, long table_type, long num_levels,
  bool relabel, bool contiguous_space, size_t max_deg, sequence<uintE>& rank, int shift_factor) {
    sequence<size_t> count;

  if (table_type == 3) {
    t.start();
    // Num levels matches, e.g., 2 for two level
    num_levels -= 1;
    if (!relabel) {
      auto rank_func = [&](uintE a, uintE b){ return rank[a] < rank[b]; };
      multitable::MHash table(r, DG, max_deg, num_levels, contiguous_space, rank_func);
      double tt = t.stop();
      std::cout << "### Table Running Time: " << tt << std::endl;
      count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank);
    } else {
      multitable::MHash table(r, DG, max_deg, num_levels, contiguous_space, std::less<uintE>());
      double tt = t.stop();
      std::cout << "### Table Running Time: " << tt << std::endl;
      count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank);
    }
  } else if (table_type == 2) {
    t.start();
    twotable::TwolevelHash table(r, DG, max_deg, contiguous_space, relabel);
    double tt = t.stop();
    std::cout << "### Table Running Time: " << tt << std::endl;
    count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank);
  } else if (table_type == 1) {
    t.start();
    onetable::OnelevelHash table(r, DG, max_deg);
    double tt = t.stop();
    std::cout << "### Table Running Time: " << tt << std::endl;
    count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank);
  } else if (table_type == 4) {
    // Num levels matches, e.g., 2 for two level
    num_levels -= 1;
    if (!relabel) {
      auto rank_func = [&](uintE a, uintE b){ return rank[a] < rank[b]; };
      multitable_nosearch::MHash table(r, DG, max_deg, num_levels, rank_func);
      double tt = t.stop();
      std::cout << "### Table Running Time: " << tt << std::endl;
      count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank);
    } else {
      multitable_nosearch::MHash table(r, DG, max_deg, num_levels, std::less<uintE>());
      double tt = t.stop();
      std::cout << "### Table Running Time: " << tt << std::endl;
      count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank);
    }
  } else if (table_type == 5) {
    t.start();
    twotable_nosearch::TwolevelHash<T> table(r, DG, max_deg, relabel, shift_factor);
    double tt = t.stop();
    std::cout << "### Table Running Time: " << tt << std::endl;
    count = NucleusDecompositionRunner(GA, DG, r, s, table, max_deg, rank);
  } 
  return count;
}

template <class Graph>
inline sequence<size_t> NucleusDecomposition(Graph& GA, size_t r, size_t s, long table_type, long num_levels,
  bool relabel, bool contiguous_space) {
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
    return (rank[u] < rank[v]) && GA.get_vertex(u).getOutDegree() >= r-1 && GA.get_vertex(v).getOutDegree() >= r-1;
  };
  // Note: If relabeling, core #s must be translated back
  auto DG = relabel ? relabel_graph(GA, rank.begin(), pack_predicate) : filterGraph(GA, pack_predicate);
  double tt_filter = t_filter.stop();
  std::cout << "### Filter Graph Running Time: " << tt_filter << std::endl;

  auto max_deg = get_max_deg3(DG);

  timer t; 
  sequence<size_t> count;

  // unsigned __int128 is 16 bytes
  // unsigned __int32 is 4 bytes
  // unsigned __int64 is 8 bytes

  // Let X be the number of bits needed to express max vertex
  // We have (r, s)
  // round_up<int>(((max(1, (r - (num_levels - 1))) * X) + 1), 8)
  // use whichever type is >= bytes than this

  int num_bits_in_n = pbbslib::log2_up(DG.n + 1);
  int num_bytes_needed = round_up<int>(((std::max(static_cast<int>(1), 
    static_cast<int>(r - (num_levels - 1))) * num_bites_in_n) + 1), 8);
  int shift_factor = num_bits_in_n; //32

  std::cout << "Num bytes needed: " << num_bytes_needed << std::endl;
  std::cout << "Num bits in n: " << shift_factor << std::endl;
  fflush(stdout);

  if (num_bytes_needed <= 4) {
    // unsigned __int32
    count = runner<unsigned __int32>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
      max_deg, rank, shift_factor);
  } else if (num_bytes_needed <= 8) {
    // unsigned __int64
    count = runner<unsigned __int64>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
      max_deg, rank, shift_factor);
  } else {
    // unsigned__int128
    count = runner<unsigned __int128>(GA, DG, r, s, table_type, num_levels, relabel, contiguous_space,
      max_deg, rank, shift_factor);
  }

   

  //table.del();
  DG.del();

  return count;
}

}