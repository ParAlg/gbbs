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

//#include "twotable.h"
#include "twotable_nosearch.h"
//#include "onetable.h"
#include "commontable.h"
#include "list_buffer.h"
#include "NucleusDecomposition_pnd.h"
#include "NucleusDecomposition_structs.h"

#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"

namespace gbbs {

  template <class Graph, class T>
  inline size_t CountCliquesNuc(Graph& DG, size_t k, size_t r, size_t max_deg, T* table, bool tmp_verify = false) {
    k--; r--;
    //timer t2; t2.start();

    if (k == 2) {
      auto counts = sequence<size_t>::from_function(DG.n, [&](size_t i){ return 0; });
      auto base_f = [&](uintE v1, uintE v2, uintE v3) {
        table->insert_twothree(v1, v2, v3, r, k);
      };
      return CountDirectedBalanced(DG, counts.begin(), base_f);
    }

    auto base_f = [&](sequence<uintE>& base){
      table->insert(base, r, k);
    };
    auto tots = sequence<size_t>(DG.n, size_t{0});

    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, DG.n, true, true); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del();
    parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
    //for(size_t i =0; i < DG.n; i++) {
    //  HybridSpace_lw* induced = new HybridSpace_lw();
    //  init_induced(induced);
        if (DG.get_vertex(i).out_degree() != 0) {
          induced->setup(DG, k, i);
          auto base = sequence<uintE>(k + 1);
          base[0] = i;
          tots[i] = NKCliqueDir_fast_hybrid_rec(DG, 1, k, induced, base_f, base);
        } else tots[i] = 0;
    //  finish_induced(induced);
    }, 1, false);
    //double tt2 = t2.stop();

    auto total_count = parlay::reduce(tots);

    if (tmp_verify && k != 2) {
      uint64_t verify_total = 0;
      auto tmp_base_f = [&](sequence<uintE>& base){
        auto idx = table->extract_indices_check(base.data(), r);
        auto count = table->get_count(idx);
        gbbs::fetch_and_add(&verify_total, count);
      };

      auto verify_init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, r, DG.n, true, true); };
      auto verify_finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del();
      parallel_for_alloc<HybridSpace_lw>(verify_init_induced, verify_finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
        if (DG.get_vertex(i).out_degree() != 0) {
          induced->setup(DG, r, i);
          auto base = sequence<uintE>(r + 1);
          base[0] = i;
          NKCliqueDir_fast_hybrid_rec(DG, 1, r, induced, tmp_base_f, base);
        }
      }, 1, false);
      verify_total = verify_total / BinomialCoefficient(k+1, r+1);
      if (verify_total != total_count){
        std::cout << "verify: " << verify_total << std::endl;
        std::cout << "total count: " << total_count << std::endl;
        fflush(stdout); exit(0);
      }
      assert(verify_total == total_count);
    }

    return total_count;
  }

template <class Graph, class Graph2, class F, class I, class T, class D, class CWP, class bucket_t>
inline size_t cliqueUpdate(Graph& G, Graph2& DG, size_t r, 
size_t k, size_t max_deg, bool label, F get_active, size_t active_size,
  size_t granularity, char* still_active, sequence<uintE> &rank, 
  sequence<double>& per_processor_counts, 
  bool do_update_changed, I& update_changed,
  T* cliques, size_t n, list_buffer& count_idxs, timer& t1,
  sequence<uintE>& inverse_rank, bool relabel, timer& t_update_d, D& cores,
  bool use_ppc,
  bool inline_hierarchy, CWP& cwp, bucket_t cur_bkt) {
  
  using W = typename Graph::weight_type;

    // Mark every vertex in the active set
    parallel_for (0, active_size, [&] (size_t j) {
      auto index = get_active(j); //cliques->find_index(get_active(j));
      still_active[index] = 1;
    }, 2048);

  auto is_active = [&](size_t index) {
    return still_active[index] == 1  || still_active[index] == 4;
  };
  auto is_inactive = [&](size_t index) {
    return still_active[index] == 2;
  };
  auto is_inactive_hierarchy = [&](size_t index) {
    return still_active[index] == 0 || still_active[index] == 3;
  };

  auto cores_func = [&](size_t a) -> bucket_t {
    if (is_inactive_hierarchy(a)) return n + 1;
    return is_active(a) ? cur_bkt : cores[a]; // 
  };

  auto update_d = [&](unsigned __int128 x, uintE* base){
    if (use_ppc) {
      cliques->extract_indices(base, is_active, is_inactive, [&](std::size_t index, double val){
      double ct = gbbs::fetch_and_add(&(per_processor_counts[index]), val);
      if (ct == 0 && val != 0) {
        count_idxs.add(index);
      }
      }, r, k, x);
      if (inline_hierarchy) {
      cliques->extract_indices_conn(base, is_inactive_hierarchy, [&](std::size_t index){
        if (index != x) {
          cwp.link(x, index, cores_func);
        }
      }, r, k, x);
    }
    } else {
      if (inline_hierarchy) {
      cliques->extract_indices_conn(base, is_inactive_hierarchy, [&](std::size_t index){
        if (index != x) {
          cwp.link(x, index, cores_func);
        }
      }, r, k, x);
    }
    cliques->extract_indices(base, is_active, is_inactive, [&](std::size_t index, double val){
        if (!is_inactive(index) && !is_active(index)) {
        cliques->update_count_atomic(index, gbbs::uintE{std::round(val)});
        if (gbbs::CAS(&(still_active[index]), char{0}, char{3}) || gbbs::CAS(&(still_active[index]), char{1}, char{4}))
          count_idxs.add(index);
        }
    }, r, k, x);
    }
  };

  // Set up space for clique counting
  auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k-r, G.n, true, true, true); };
  auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };
t1.start();
  // Clique count updates
  //std::cout << "Start setup nucleus" << std::endl; fflush(stdout);
  if (k - r == 1) {
    // For each vert from 0 to active_size, intersect + place the intersection
    // in base[1], and then call update_d
    if (k == 2 && r == 1) { // This is (2, 3)

      parallel_for(0, active_size, [&](size_t i){
        auto x = get_active(i);

        auto update_d_twothree = [&](uintE v1, uintE v2, uintE v3){
          if (inline_hierarchy) {
           std::vector<uintE> base = {v1, v2, v3};
      cliques->extract_indices_conn(base.data(), is_inactive_hierarchy, [&](std::size_t index){
        if (index != x) {
          cwp.link(x, index, cores_func);
        }
      }, r, k, x);
    }
        cliques->extract_indices_twothree(v1, v2, v3, is_active, is_inactive,
          [&](std::size_t index, double val){
            if (use_ppc) {
              double ct = gbbs::fetch_and_add(&(per_processor_counts[index]), val);
          if (ct == 0 && val != 0) count_idxs.add(index);
            } else {
              if (!is_inactive(index) && !is_active(index)) {
              cliques->update_count_atomic(index, gbbs::uintE{std::round(val)});
        if (gbbs::CAS(&(still_active[index]), char{0}, char{3}) || gbbs::CAS(&(still_active[index]), char{1}, char{4}))
          count_idxs.add(index);
              }
            }
        }, r, k);

         
      };

        std::tuple<uintE, uintE> v1v2 = cliques->extract_clique_two(x, k);
        uintE u = relabel ? inverse_rank[std::get<0>(v1v2)] : std::get<0>(v1v2);
        uintE v = relabel ? inverse_rank[std::get<1>(v1v2)] : std::get<1>(v1v2);
        auto v_v = G.get_vertex(v);
        auto process_f = [&](uintE a, uintE b, uintE intersect_w){
          uintE v3 = relabel ? rank[intersect_w] : intersect_w;
          update_d_twothree(std::get<0>(v1v2), std::get<1>(v1v2), v3);
        };
        auto their_neighbors = v_v.out_neighbors();
        G.get_vertex(u).out_neighbors().intersect_f_par(&their_neighbors, process_f);

      });
    } else if (k == 3 && r == 2) { // (3, 4))
      

      /*auto init_intersect = [&](IntersectSpace* arr){ arr->alloc(G.n); };
      auto finish_intersect = [&](IntersectSpace* arr){
        if (arr != nullptr){
          if (arr->labels != nullptr) { free(arr->labels); arr->labels = nullptr; }
          delete arr;
        }
      };
      parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active_size,
                                     [&](size_t i, HybridSpace_lw* induced) {*/
      //ThreadLocalObj<IntersectSpace> thread_local_is = ThreadLocalObj<IntersectSpace>();
      parallel_for(0, active_size, [&](size_t i) {
      //  auto is_pair = thread_local_is.reserve();
      //  IntersectSpace* is = is_pair.second;
        static thread_local IntersectSpace* is = nullptr;
        if (is == nullptr) is = new IntersectSpace();
        is->alloc(G.n);
        auto labels = is->labels; //induced->old_labels;
        auto x = get_active(i);

        auto update_d_threefour = [&](uintE v1, uintE v2, uintE v3, uintE v4){
          if (inline_hierarchy) {
           std::vector<uintE> base = {v1, v2, v3, v4};
      cliques->extract_indices_conn(base.data(), is_inactive_hierarchy, [&](std::size_t index){
        if (index != x) {
         cwp.link(x, index, cores_func);
        }
      }, r, k, x);
    }

        //uintE base[4]; base[0] = v1; base[1] = v2; base[3] = v4; base[2] = v3;
        cliques->extract_indices_threefour(v1, v2, v3, v4, is_active, is_inactive,
          [&](std::size_t index, double val){
            if (use_ppc) {
          double ct = gbbs::fetch_and_add(&(per_processor_counts[index]), val);
          if (ct == 0 && val != 0) count_idxs.add(index);
            } else {
              if (!is_inactive(index) && !is_active(index)) {
              cliques->update_count_atomic(index, gbbs::uintE{std::round(val)});
        if (gbbs::CAS(&(still_active[index]), char{0}, char{3}) || gbbs::CAS(&(still_active[index]), char{1}, char{4}))
          count_idxs.add(index);
              }
            }
        }, r, k);

        
      };


        std::tuple<uintE, uintE, uintE> v1v2v3 = cliques->extract_clique_three(x, k);
        uintE u = relabel ? inverse_rank[std::get<0>(v1v2v3)] : std::get<0>(v1v2v3);
        auto map_label_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
          uintE actual_ngh = relabel ? rank[ngh] : ngh;
          assert(labels[actual_ngh] == 0);
          labels[actual_ngh] = 1;
        };
        G.get_vertex(u).out_neighbors().map(map_label_f, false);
        uintE v = relabel ? inverse_rank[std::get<1>(v1v2v3)] : std::get<1>(v1v2v3);
        auto map_label_inner_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
          uintE actual_ngh = relabel ? rank[ngh] : ngh;
          assert(labels[actual_ngh] <= 1);
          if (labels[actual_ngh] > 0) labels[actual_ngh]++;
        };
        G.get_vertex(v).out_neighbors().map(map_label_inner_f, false);
        v = relabel ? inverse_rank[std::get<2>(v1v2v3)] : std::get<2>(v1v2v3);
        auto map_label_inner_f2 = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
          uintE actual_ngh = relabel ? rank[ngh] : ngh;
          assert(labels[actual_ngh] <= 2);
          if (labels[actual_ngh] > 0) labels[actual_ngh]++;
        };
        G.get_vertex(v).out_neighbors().map(map_label_inner_f2, false);
        // Any vtx with labels[vtx] = k - 1 is in the intersection
        auto map_update_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
          uintE actual_ngh = relabel ? rank[ngh] : ngh;
          if (labels[actual_ngh] == k) {
            assert(is_edge(G, ngh, relabel ? inverse_rank[std::get<1>(v1v2v3)] : std::get<1>(v1v2v3)));
            update_d_threefour(std::get<0>(v1v2v3), std::get<1>(v1v2v3), std::get<2>(v1v2v3), actual_ngh);
          }
          labels[actual_ngh] = 0;
        };
        G.get_vertex(u).out_neighbors().map(map_update_f, false);
        //thread_local_is.unreserve(is_pair.first);
      },1, true);
    } else { // This is not (2, 3)
      /*auto init_intersect = [&](IntersectSpace* arr){
        arr->alloc(G.n);
      };
      auto finish_intersect = [&](IntersectSpace* arr){
        if (arr != nullptr){
          if (arr->labels != nullptr) {
            free(arr->labels);
            arr->labels = nullptr;
          }
          delete arr;
        }
      };*/
      /*parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active_size,
                                     [&](size_t i, HybridSpace_lw* induced) {*/
      //ThreadLocalObj<IntersectSpace> thread_local_is = ThreadLocalObj<IntersectSpace>();
      parallel_for(0, active_size, [&](size_t i) {
        //auto is_pair = thread_local_is.reserve();
        //IntersectSpace* is = is_pair.second;
        static thread_local IntersectSpace* is = nullptr;
        if (is == nullptr) is = new IntersectSpace();
        is->alloc(G.n);
        auto labels = is->labels; //induced->old_labels;
        auto x = get_active(i);
        uintE base[10];
        cliques->extract_clique(x, base, G, k);

        // Sequentially find min deg and swap with 0
        uintE u = relabel ? inverse_rank[base[0]] : base[0];
        auto min_deg = G.get_vertex(u).out_degree();
        size_t u_idx = 0;
        for (size_t j = k; j > 1; j--) {
          uintE v = relabel ? inverse_rank[base[j]] : base[j];
          auto v_deg = G.get_vertex(v).out_degree();
          if (v_deg < min_deg) {
            u_idx = j;
            min_deg = v_deg;
          }
        }
        // Swap u_idx and 0
        if (u_idx != 0) {
          auto tmp = base[u_idx];
          base[u_idx] = base[0];
          base[0] = tmp;
        }
        u = relabel ? inverse_rank[base[0]] : base[0];
        auto map_label_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
          uintE actual_ngh = relabel ? rank[ngh] : ngh;
          labels[actual_ngh] = 1;
        };
        G.get_vertex(u).out_neighbors().map(map_label_f, false);
        for (size_t j = k; j > 1; j--) {
          uintE v = relabel ? inverse_rank[base[j]] : base[j];
          auto map_label_inner_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
            uintE actual_ngh = relabel ? rank[ngh] : ngh;
            if (labels[actual_ngh] > 0) labels[actual_ngh]++;
          };
          G.get_vertex(v).out_neighbors().map(map_label_inner_f, false);
        }
        // Any vtx with labels[vtx] = k - 1 is in the intersection
        auto map_update_f = [&] (const uintE& src, const uintE& ngh, const W& wgh) {
          uintE actual_ngh = relabel ? rank[ngh] : ngh;
          if (labels[actual_ngh] == k) {
            base[1] = actual_ngh;
            update_d(x, base);
          }
          labels[actual_ngh] = 0;
        };
        G.get_vertex(u).out_neighbors().map(map_update_f, false);
        //thread_local_is.unreserve(is_pair.first);
      },1, true);
    }

  } else {

  //parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active_size,
  //                                   [&](size_t i, HybridSpace_lw* induced) {
  //ThreadLocalObj<HybridSpace_lw> thread_local_is = ThreadLocalObj<HybridSpace_lw>();
  parallel_for(0, active_size, [&](size_t i) {
    //auto is_pair = thread_local_is.reserve();
    //HybridSpace_lw* induced = is_pair.second;
    static thread_local HybridSpace_lw* induced = nullptr;
    if (induced == nullptr) induced = new HybridSpace_lw();
    induced->alloc(max_deg, k-r, G.n, true, true, true);
  //for(size_t i =0; i < active_size; i++) {
  //  HybridSpace_lw* induced = new HybridSpace_lw();
    init_induced(induced);

    auto x = get_active(i);

      auto update_d_bind = [&](uintE* base){
         if (inline_hierarchy) {
      cliques->extract_indices_conn(base, is_inactive_hierarchy, [&](std::size_t index){
        if (index != x) {
          cwp.link(x, index, cores_func);
        }
      }, r, k, x);
    }
    
    cliques->extract_indices(base, is_active, is_inactive, [&](std::size_t index, double val){
      if (use_ppc) {
      double ct = gbbs::fetch_and_add(&(per_processor_counts[index]), val);
      if (ct == 0 && val != 0) {
        count_idxs.add(index);
      }
      } else {
        if (!is_inactive(index) && !is_active(index)) {
        cliques->update_count_atomic(index, gbbs::uintE{std::round(val)});
        if (gbbs::CAS(&(still_active[index]), char{0}, char{3}) || gbbs::CAS(&(still_active[index]), char{1}, char{4}))
          count_idxs.add(index);
        }
      }
    }, r, k, x);

   
  };


    //auto base = sequence<uintE>(k + 1);
    uintE base[10];
    cliques->extract_clique(x, base, G, k);
    // Fill base[k] ... base[k-r+1] and base[0]
    if (relabel) {
      auto g_vert_map = [&](uintE vert){return rank[vert];};
      auto inverse_g_vert_map = [&](uintE vert){return inverse_rank[vert];};
      induced->setup_nucleus(G, DG, k, base, r, g_vert_map, inverse_g_vert_map);
    } else {
      auto g_vert_map = [&](uintE vert){return vert;};
      induced->setup_nucleus(G, DG, k, base, r, g_vert_map, g_vert_map);
    }
    NKCliqueDir_fast_hybrid_rec(DG, 1, k-r, induced, update_d_bind, base);
    //thread_local_is.unreserve(is_pair.first);
  //  finish_induced(induced);
  }, 1, true); //granularity
  //std::cout << "End setup nucleus" << std::endl; fflush(stdout);
  }
t1.stop();

  std::size_t num_count_idxs = 0;
if (do_update_changed && use_ppc) {
t_update_d.start();
  // Perform update_changed on each vertex with changed clique counts

  num_count_idxs = count_idxs.filter(update_changed, per_processor_counts);
  count_idxs.reset();
t_update_d.stop();
}

  // Mark every vertex in the active set as deleted
  parallel_for (0, active_size, [&] (size_t j) {
    auto index = get_active(j); 
    still_active[index] = 2;}, 2048);

  return num_count_idxs; //count_idxs[0];
}

// Sort vertices from highest core # to lowest core #
// Take each "bucket" of vertices in each core #
// Maintain connectivity UF throughout all rounds
// Run connectivity considering only vertices previously considered, or vertices
// in the current bucket
// k and r here should be the same as that used in cliqueUpdate
template <class bucket_t, class Graph, class Graph2, class Table>
inline std::vector<uintE> construct_nd_connectivity(sequence<bucket_t>& cores, Graph& GA, Graph2& DG,
  size_t r, size_t k, Table& table, sequence<uintE>& rank, bool relabel){
  auto n = table.return_total();
  // TODO(jeshi): don't repeat work from peel
  size_t max_deg = induced_hybrid::get_max_deg(GA);
  sequence<uintE> inverse_rank;
  if (relabel) {
    // This maps a DG vertex to a G vertex
    inverse_rank = sequence<uintE>(GA.n);
    parallel_for(0, GA.n, [&](size_t i){
      inverse_rank[rank[i]] = i;
    });
  }

  // Sort vertices from highest core # to lowest core #
  auto get_core = [&](uintE p, uintE q){ return cores[p] > cores[q]; };
  auto sorted_vert = sequence<uintE>::from_function(n, [&](size_t i) { return i; });
  parlay::sample_sort_inplace(make_slice(sorted_vert), get_core);


  sequence<uintE> mark_keys(n + 1);
  parallel_for(0, n, [&](std::size_t i) {
    if (i != 0 && cores[sorted_vert[i]] == cores[sorted_vert[i-1]])
      mark_keys[i] = UINT_E_MAX;
    else
      mark_keys[i] = i;
  });
  mark_keys[n] = n;
  auto vert_buckets = 
      parlay::filter(mark_keys, [&](uintE x) -> bool { return x != UINT_E_MAX; });

  //auto vert_buckets = GetBoundaryIndices<uintE>(n, [&](size_t i, size_t j){
  //  return cores[sorted_vert[i]] == cores[sorted_vert[j]];
  //});

  auto uf = gbbs::simple_union_find::SimpleUnionAsyncStruct(n);
  // TODO(jeshi): This isn't parallel, but I'd just like easy resizing atm
  std::vector<uintE> connectivity_tree(n);
  sequence<uintE> prev_parent = sequence<uintE>::from_function(n, [&](size_t i){ return i; });
  uintE prev_max_parent = n;
  //uintE prev_offset = 0;
  parallel_for(0, n, [&](std::size_t i){connectivity_tree[i] = UINT_E_MAX;});

  auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k-r, GA.n, true, true, true); };
  auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };

  for (size_t i = 0; i < vert_buckets.size()-1; i++) {
    size_t start_index = vert_buckets[i];
    size_t end_index = vert_buckets[i + 1];

    auto first_x = sorted_vert[start_index];
    auto first_current_core = cores[first_x];
    // TODO: to get rid of this, gotta check for valid indices
    if (first_current_core != 0) {

      auto is_active = [&](size_t index) {
        return cores[index] == first_current_core;
      };
      auto is_inactive = [&](size_t index) {
        return cores[index] < first_current_core;
      };
      //std::cout << "Core: " << first_current_core << std::endl; fflush(stdout);

    parallel_for(start_index, end_index, [&](size_t j){
      // Vertices are those given by sorted_vert[j]
      auto x = sorted_vert[j];
      if (table.is_valid(x)) {
      auto current_core = cores[x];
      assert(current_core == first_current_core);

      // Steps:
      // 1. Extract vertices given by index -- these will be the r-clique X
      // 2. Find all s-clique containing the r-clique X
      // 3. For each s-clique containing X, iterate over all combinations of
      // r vertices in that s-clique; these form r-clique X'
      // 4. Union X and X'
    static thread_local HybridSpace_lw* induced = nullptr;
    if (induced == nullptr) induced = new HybridSpace_lw();
    induced->alloc(max_deg, k-r, GA.n, true, true, true);
    init_induced(induced);

    auto update_d_bind = [&](uintE* base){
      table.extract_indices_conn(base, is_inactive, [&](std::size_t index){
        if (index != x) {
          uf.unite(x, index);
        }
      }, r, k, x);
    };
    uintE base[10];
    table.extract_clique(x, base, GA, k);
    // Fill base[k] ... base[k-r+1] and base[0]
    if (relabel) {
      auto g_vert_map = [&](uintE vert){return rank[vert];};
      auto inverse_g_vert_map = [&](uintE vert){return inverse_rank[vert];};
      induced->setup_nucleus(GA, DG, k, base, r, g_vert_map, inverse_g_vert_map);
    } else {
      auto g_vert_map = [&](uintE vert){return vert;};
      induced->setup_nucleus(GA, DG, k, base, r, g_vert_map, g_vert_map);
    }
    NKCliqueDir_fast_hybrid_rec(DG, 1, k-r, induced, update_d_bind, base);
      }
    }, 1, true); // granularity
    }

    connectivity_tree.resize(prev_max_parent, UINT_E_MAX);

        // TODO(jeshi): Here we gotta update the connectivity tree
    parallel_for(0, n, [&](size_t l) { gbbs::simple_union_find::find_compress(l, uf.parents); });
    //for (std::size_t yyy = 0 ; yyy < uf.parents.size(); yyy++) {
    //  std::cout << "yyy: " << yyy << ", " << uf.parents[yyy] << std::endl;
    //}
    // It looks like at least for simple_union_find, uf.parents is guaranteed to be in the range 0
    // to n, which we could inefficiently use...but maybe we have enough time here to compress
    //auto old_max_parent = 1 + parlay::reduce(make_slice(uf.parents), parlay::maxm<uintE>());
    // This is to compress the space of parents to at least be contiguous
    sequence<uintE> map_parents = sequence<uintE>::from_function(n, [&](std::size_t l){return 0;});
    parallel_for(0, n, [&](size_t l){
      if (table.is_valid(l) && cores[l] >= first_current_core) map_parents[uf.parents[l]] = 1;
    });
    auto max_parent = parlay::scan_inplace(make_slice(map_parents));
    //std::cout << "Max parent: " << max_parent << std::endl; fflush(stdout);
    //***: auto max_parent = n;

    //std::cout << "Max parent: " << max_parent << std::endl;
    //std::cout << "Prev max parent: " << prev_max_parent << std::endl; fflush(stdout);
    parallel_for(0, n, [&](size_t l){
      // One check is that once one person changes this value, anyone else trying to change
      // this value should not succeed (they should be trying to change this to be the same val)
      if (table.is_valid(l) && cores[l] >= first_current_core) { //***: if (table.is_valid(l)) {
        //assert(uf.parents[l] < old_max_parent);
        assert(prev_parent[l] < prev_max_parent);
        connectivity_tree[prev_parent[l]] = prev_max_parent + map_parents[uf.parents[l]];  //***: uf.parents[l];
        // Update previous parent
        prev_parent[l] = connectivity_tree[prev_parent[l]];
      }
    });
    // Update previous max parent
    //prev_offset = prev_max_parent;
    prev_max_parent += max_parent;
    //std::cout << "Prev max parent 2: " << prev_max_parent << std::endl; fflush(stdout);
  }
  //for (std::size_t i = 0; i < n; i++) {
  //  std::cout << "core: " << i << ", " << cores[i] << std::endl;
  //}
  //for (std::size_t i = 0; i < connectivity_tree.size(); i++) {
  //  std::cout << "i: " << i << ", parent: " << connectivity_tree[i] << std::endl;
  //}
  return connectivity_tree; //uf.finish();
}

template <typename bucket_t, class Graph, class Graph2, class T, class CWP>
sequence<bucket_t> Peel(Graph& G, Graph2& DG, size_t r, size_t k, 
  T* cliques, sequence<uintE> &rank, size_t fake_efficient, bool relabel, 
  bool use_compress, bool inline_hierarchy, CWP& connect_while_peeling,
  size_t num_buckets=16) {
    size_t efficient = fake_efficient;
    if (fake_efficient == 3) efficient = 1;
  sequence<uintE> inverse_rank;
  if (relabel) {
    // This maps a DG vertex to a G vertex
    inverse_rank = sequence<uintE>(G.n);
    parallel_for(0, G.n, [&](size_t i){
      inverse_rank[rank[i]] = i;
    });
  }
    k--; r--;

  timer t2; t2.start();

  size_t num_entries = cliques->return_total();
  std::cout << "num entries: " << num_entries << std::endl;
  auto D = sequence<bucket_t>::from_function(num_entries, [&](size_t i) -> bucket_t { 
    return cliques->get_count(i);
  });

  auto num_entries_filter = num_entries;
  if (efficient == 1) num_entries_filter += num_workers() * 1024;
  else if (efficient == 4) num_entries_filter = 1 + 10000 * ((1 + (num_entries / 10000) / 1024) * 1024  + 1024* num_workers());
  auto D_filter = sequence<std::tuple<uintE, bucket_t>>(num_entries_filter);

  auto b = make_vertex_custom_buckets<bucket_t>(num_entries, D, increasing, num_buckets);

  auto per_processor_counts = sequence<double>(num_entries , static_cast<double>(0));
  
  list_buffer count_idxs(num_entries, efficient);

  char* still_active = (char*) calloc(num_entries, sizeof(char));
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active?

  /*timer t_compress;
  compress_utils compress_util;
  if (r == 1 && k == 2 && use_compress) {
    compress_util = compress_utils(G);
  }*/

  timer t_extract;
  timer t_count;
  timer t_update;
  timer t_x;
  timer t_update_d;

  size_t rounds = 0;
  size_t finished = 0;
  bucket_t cur_bkt = 0;
  bucket_t max_bkt = 0;
  bucket_t prev_bkt = 0;
  double max_density = 0;
  bool use_max_density = false;
  size_t iter = 0;
  size_t total_active_size = 0;
  //std::vector<unsigned __int128> total_active;

  while (finished < num_entries) {
    t_extract.start();
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers; //vertexSubset(num_entries, bkt.identifiers);
    auto active_size = active.size();
    cur_bkt = bkt.id;
    t_extract.stop();

    auto get_active = [&](size_t j) -> unsigned __int128 { 
      return active[j]; }; //active.vtx(j); };
    
    if (active_size == 0) continue;

    finished += active_size;

    if (inline_hierarchy && prev_bkt != cur_bkt && cur_bkt != 0) {
      connect_while_peeling.init(cur_bkt);
    }

    /*if (inline_hierarchy && prev_bkt != cur_bkt && prev_bkt != 0) {
      auto hierarchy_is_inactive = [&](std::size_t j){ return still_active[j] == 0 || still_active[j] == 1; };
      auto hierarchy_is_active = [&](std::size_t j){
        return still_active[j] == 2; // && D[j] == prev_bkt;
      };
      auto get_total_active = [&](size_t j) -> unsigned __int128 { 
        return total_active[j]; }; //active.vtx(j); };
      // here, is_active is same, but is_inactive is never been active before
      connect_while_peeling.update_cores(prev_bkt, get_total_active, total_active_size, G, 
      DG, r, k, *cliques, rank, relabel, max_deg, inverse_rank, 
      hierarchy_is_active, hierarchy_is_inactive, D);
      total_active_size = 0;
      total_active = std::vector<unsigned __int128>(0);
    }*/

    //if (cur_bkt >= std::numeric_limits<bucket_t>::max() - 1) continue;

    max_bkt = std::max(cur_bkt, max_bkt);
    if (cur_bkt == 0) { // || finished == num_entries) { 
      parallel_for (0, active_size, [&] (size_t j) {
        auto index = get_active(j);
        still_active[index] = 2;
        cliques->set_count(index, 0); // UINT_E_MAX
        D[index] = std::max((bucket_t) 0, (bucket_t) cur_bkt);
      }, 2048);
      /*if (cur_bkt != 0 && inline_hierarchy && finished == num_entries) {
        auto hierarchy_is_inactive = [&](std::size_t j){ return still_active[j] == 0 || still_active[j] == 1; };
        auto hierarchy_is_active = [&](std::size_t j){
          return still_active[j] == 2; //&& D[j] == cur_bkt;
        };
        total_active_size += active_size;
        for (size_t i = 0; i < active_size; i++) {
          total_active.push_back(active[i]);
        }
        auto get_total_active = [&](size_t j) -> unsigned __int128 { 
          return total_active[j]; }; //active.vtx(j); };
        // here, is_active is same, but is_inactive is never been active before
        connect_while_peeling.update_cores(cur_bkt, get_total_active, total_active_size, G, 
        DG, r, k, *cliques, rank, relabel, max_deg, inverse_rank, 
        hierarchy_is_active, hierarchy_is_inactive, D);
        total_active_size = 0;
        total_active = std::vector<unsigned __int128>(0);
    }*/
      continue;
    }

    //std::cout << "k = " << cur_bkt << " iter = " << iter << " #edges = " << active_size << std::endl;
    //std::cout << "Finished: " << finished << ", num_entries: " << num_entries << std::endl;
    iter++;

    size_t granularity = (cur_bkt * active_size < 10000) ? 1024 : 1;

    size_t filter_size = 0;

      auto update_changed = [&](sequence<double>& ppc, size_t i, uintE v){
        if (v == UINT_E_MAX) {
          D_filter[i] = std::make_tuple(num_entries + 1, 0);
          return;
        }
        //assert(ppc[v] != 0);
        if (ppc[v] == 0) D_filter[i] = std::make_tuple(num_entries + 1, 0);
        else {
          bucket_t deg = D[v];
          assert(deg > cur_bkt);
          auto clique_count = cliques->get_count(v);
          //if (std::round(ppc[v]) > clique_count){
            //std::cout << "PPC: " << std::round(ppc[v]) << ", count: " << clique_count << ", v: " << v << std::endl;
            //fflush(stdout);
            //exit(0);
          //}
          //assert(std::round(ppc[v]) <= clique_count);
          
          auto val = clique_count - std::round(ppc[v]);
          cliques->set_count(v, val);
          if (deg > cur_bkt) {
            bucket_t new_deg = std::max((bucket_t) val, (bucket_t) cur_bkt);
            D[v] = new_deg;
            // store (v, bkt) in an array now, pass it to apply_f below instead of what's there right now -- maybe just store it in D_filter?
            //cliques->set_count(v, (size_t) new_deg);
            D_filter[i] = std::make_tuple(v, b.get_bucket(deg, new_deg));
          } else D_filter[i] = std::make_tuple(num_entries + 1, 0);
        }
        ppc[v] = 0;
      };
    t_count.start();
    count_idxs.resize(active_size, k, r, cur_bkt);
    /*if (fake_efficient == 3) {
      filter_size = cliqueUpdatePND(G, DG, r, k, max_deg, true, get_active, active_size, 
     granularity, still_active, rank, per_processor_counts,
      true, update_changed, cliques, num_entries, count_idxs, t_x, inverse_rank, relabel,
      t_update_d);
    } else {*/
     filter_size = cliqueUpdate(G, DG, r, k, max_deg, true, get_active, active_size, 
     granularity, still_active, rank, per_processor_counts,
      true, update_changed, cliques, num_entries, count_idxs, t_x, inverse_rank, relabel,
      t_update_d, D, true, inline_hierarchy, connect_while_peeling, cur_bkt);
    //}
      t_count.stop();

    auto apply_f = [&](size_t i) -> std::optional<std::tuple<uintE, bucket_t>> {
      auto v = std::get<0>(D_filter[i]);
      bucket_t bucket = std::get<1>(D_filter[i]);
      if (v != num_entries + 1) {
        if (still_active[v] != 2 && still_active[v] != 1) return wrap(v, bucket);
      }
      return std::nullopt;
    };

    t_update.start();
    b.update_buckets(apply_f, filter_size);

    /*parallel_for (0, active_size, [&] (size_t j) {
      auto index = get_active(j);
      //auto index = cliques->find_index(v);
      cliques->clear_count(index);
      //cliques[active.vtx(j)] = 0;
      }, 2048);*/
    t_update.stop();

    rounds++;
    prev_bkt = cur_bkt;
    /*total_active_size += active_size;
    for (size_t i = 0; i < active_size; i++) {
      total_active.push_back(active[i]);
    }*/

    /*if (r == 1 && k == 2 && use_compress) {
      bool to_compress = compress_util.update_del_edges(active);
      if (to_compress) {
        //std::cout << "Start compress" << std::endl;
        if (!relabel) {
          auto map_g_dg = [](uintE x) {return x;};
          CompressOut(compress_util, G, map_g_dg, map_g_dg,still_active, t_compress, cliques);
        } else {
          auto g_to_dg = [&](uintE x) {return rank[x];};
          auto dg_to_g = [&](uintE x) {return inverse_rank[x];};
          CompressOut(compress_util, G, g_to_dg, dg_to_g,still_active, t_compress, cliques);
        }
        //std::cout << "End compress" << std::endl;
      }
    }*/
  
  }

  

  t_extract.next("### Peel Extract time: ");
  t_count.next("### Peel Count time: ");
  t_update.next("### Peel Update time: ");
  //t_compress.reportTotal("### Compress time: ");
  t_x.next("Inner counting: ");
  t_update_d.next("Update d time: ");

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  std::cout.precision(17);
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "clique core: " << max_bkt << std::endl;
  if (use_max_density) std::cout << "max density: " << max_density << std::endl;

  //b.del();
  free(still_active);

  return D;
}


//*************************************************SPACE EFFICIENT CODE******************


template <typename bucket_t, typename iden_t, class Graph, class Graph2, class T, class CWP>
sequence<bucket_t> Peel_space_efficient(Graph& G, Graph2& DG, size_t r, size_t k, 
  T* cliques, sequence<uintE> &rank, size_t fake_efficient, bool relabel, 
  bool use_compress, bool inline_hierarchy, CWP& connect_while_peeling,
  size_t num_buckets=16) {
    size_t efficient = fake_efficient;
    if (fake_efficient == 3) efficient = 5;
  sequence<uintE> inverse_rank;
  if (relabel) {
    // This maps a DG vertex to a G vertex
    inverse_rank = sequence<uintE>(G.n);
    parallel_for(0, G.n, [&](size_t i){
      inverse_rank[rank[i]] = i;
    });
  }
    k--; r--;
  timer t2; t2.start();

  size_t num_entries = cliques->return_total();
  std::cout << "num entries: " << num_entries << std::endl;
  auto D = sequence<bucket_t>::from_function(num_entries, [&](size_t i) -> bucket_t { 
    return cliques->get_count(i);
  });

  //auto num_entries_filter = num_entries;
  //if (efficient == 1) num_entries_filter += num_workers() * 1024;
  //else if (efficient == 4) num_entries_filter = 1 + 10000 * ((1 + (num_entries / 10000) / 1024) * 1024  + 1024* num_workers());
  //std::cout << "created 1 " << std::endl; fflush(stdout);

  auto b = buckets<sequence<bucket_t>, iden_t, bucket_t>(num_entries, D, increasing, num_buckets);
  //make_vertex_custom_buckets<bucket_t>(num_entries, D, increasing, num_buckets);
  //std::cout << "created 3 " << std::endl; fflush(stdout);

  auto per_processor_counts = sequence<double>(0);
  //std::cout << "created 4 " << std::endl; fflush(stdout);
  
  list_buffer count_idxs(num_entries, efficient);
  //std::cout << "created 5 " << std::endl; fflush(stdout);

  char* still_active = (char*) calloc(num_entries, sizeof(char));
  //std::cout << "created 6 " << std::endl; fflush(stdout);
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active?

  timer t_extract;
  timer t_count;
  timer t_update;
  timer t_x;
  timer t_update_d;

  size_t rounds = 0;
  size_t finished = 0;
  bucket_t cur_bkt = 0;
  bucket_t max_bkt = 0;
  double max_density = 0;
  bool use_max_density = false;
  size_t iter = 0;
  bucket_t prev_bkt = 0;

  while (finished < num_entries) {
    t_extract.start();
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers; //vertexSubset(num_entries, bkt.identifiers);
    auto active_size = active.size();
    cur_bkt = bkt.id;
    t_extract.stop();

    auto get_active = [&](size_t j) -> unsigned __int128 { 
      return active[j]; }; //active.vtx(j); };

    if (active_size == 0) continue;

    finished += active_size;

    if (inline_hierarchy && prev_bkt != cur_bkt && cur_bkt != 0) {
      connect_while_peeling.init(cur_bkt);
    }

    max_bkt = std::max(cur_bkt, max_bkt);
    if (cur_bkt == 0) { //|| finished == num_entries
      parallel_for (0, active_size, [&] (size_t j) {
        auto index = get_active(j);
        still_active[index] = 2;
        cliques->set_count(index, UINT_E_MAX);
        D[index] = std::max((bucket_t) 0, (bucket_t) cur_bkt);
      }, 2048);
      continue;
    }

    //std::cout << "k = " << cur_bkt << " iter = " << iter << " #edges = " << active_size << std::endl;
    //std::cout << "Finished: " << finished << ", num_entries: " << num_entries << std::endl;
    iter++;

    size_t granularity = (cur_bkt * active_size < 10000) ? 1024 : 1;

    size_t filter_size = 0;

      auto update_changed = [&](sequence<double>& ppc, size_t i, uintE v){

      };
    t_count.start();
    count_idxs.resize(active_size, k, r, cur_bkt);
    /*if (fake_efficient == 3) {
      filter_size = cliqueUpdatePND(G, DG, r, k, max_deg, true, get_active, active_size, 
     granularity, still_active, rank, per_processor_counts,
      false, update_changed, cliques, num_entries, count_idxs, t_x, inverse_rank, relabel,
      t_update_d, false);
    } else {*/
     filter_size = cliqueUpdate(G, DG, r, k, max_deg, true, get_active, active_size, 
     granularity, still_active, rank, per_processor_counts,
      false, update_changed, cliques, num_entries, count_idxs, t_x, inverse_rank, relabel,
      t_update_d, D, false, inline_hierarchy, connect_while_peeling, cur_bkt);
    //}
      t_count.stop();

      //std::cout << "FLAG 1" << std::endl; fflush(stdout);

  // Perform update_changed on each vertex with changed clique counts
  //std::size_t num_count_idxs = 0;
  //num_count_idxs = count_idxs.filter(update_changed, per_processor_counts);
  //count_idxs.reset();
  t_update.start();
    auto num_count_idxs = count_idxs.num_entries();
    auto D_filter = sequence<bucket_t>(num_count_idxs);
    
    parallel_for(0, num_count_idxs, [&](size_t i){
      auto v = count_idxs.get_v(i);
      
        if (v == UINT_E_MAX) {
          //v = num_entries + 1;
          count_idxs.void_v(i, v);
          D_filter[i] = num_entries + 1;
        } 
        else {
          bucket_t deg = D[v];
          //assert(deg > cur_bkt);
          auto val = cliques->get_count(v);
          if (deg > cur_bkt) {
            bucket_t new_deg = std::max((bucket_t) val, (bucket_t) cur_bkt);
            D[v] = new_deg;
            D_filter[i] = b.get_bucket(deg, new_deg); //std::make_tuple(v, 
          } else {
            //v = num_entries + 1;
            count_idxs.void_v(i, v);
            D_filter[i] = num_entries + 1; //std::make_tuple(num_entries + 1, 0);
          }
        }

        if (v != UINT_E_MAX) {
          gbbs::CAS(&(still_active[v]), char{3}, char{0});
          gbbs::CAS(&(still_active[v]), char{4}, char{1});
        }
    });
    //std::cout << "FLAG 2" << std::endl; fflush(stdout);
    auto apply_f = [&](size_t i) -> std::optional<std::tuple<uintE, bucket_t>> {
      auto v = count_idxs.get_v(i);
      if (v != UINT_E_MAX) {
        //if (v >= D.size()) {std::cout << "v: " << v << ", size: " << D.size() << std::endl; fflush(stdout);}
        //assert(v < D.size());
      bucket_t bucket = D_filter[i];
      //assert(bucket != num_entries + 1);
        if (still_active[v] != 2 && still_active[v] != 1 && still_active[v] != 4) return wrap(v, bucket);
      }
      return std::nullopt;
    };

    
    b.update_buckets(apply_f, num_count_idxs);
    //std::cout << "FLAG 3" << std::endl; fflush(stdout);
    count_idxs.reset();
    //std::cout << "FLAG 4" << std::endl; fflush(stdout);
    t_update.stop();

    rounds++;
    prev_bkt = cur_bkt;
  
  }

  t_extract.next("### Peel Extract time: ");
  t_count.next("### Peel Count time: ");
  t_update.next("### Peel Update time: ");
  //t_compress.reportTotal("### Compress time: ");
  t_x.next("Inner counting: ");
  t_update_d.next("Update d time: ");

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  std::cout.precision(17);
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "clique core: " << max_bkt << std::endl;
  if (use_max_density) std::cout << "max density: " << max_density << std::endl;

  //b.del();
  free(still_active);

  return D;
}


//******************************APPROX CODE**********************************

template <typename bucket_t, typename iden_t, class Graph, class Graph2, class T, class CWP>
sequence<bucket_t> ApproxPeel_space_efficient(Graph& G, Graph2& DG, size_t r, size_t k, 
  T* cliques, sequence<uintE> &rank, size_t fake_efficient, bool relabel, 
  bool use_compress, bool inline_hierarchy, CWP& connect_while_peeling,
  size_t num_buckets=16, double eps = 0.2, double delta = 0.1, bool use_pow = false) {
    std::cout << "Eps: " << eps << ", delta: " << delta << std::endl;

    size_t efficient = fake_efficient;
    if (fake_efficient == 3) efficient = 5;
  sequence<uintE> inverse_rank;
  if (relabel) {
    // This maps a DG vertex to a G vertex
    inverse_rank = sequence<uintE>(G.n);
    parallel_for(0, G.n, [&](size_t i){
      inverse_rank[rank[i]] = i;
    });
  }
    k--; r--;
  timer t2; t2.start();

  size_t num_entries = cliques->return_total();
  std::cout << "num entries: " << num_entries << std::endl;
  double one_plus_delta = log(1 + delta);
  auto get_bucket = [&](size_t deg) -> uintE {
    return ceil(log(1 + deg) / one_plus_delta);
  };
  auto D = sequence<bucket_t>::from_function(num_entries, [&](size_t i) -> bucket_t {
    auto deg = cliques->get_count(i);
    if (deg == 0) return 0; //return deg;
    if (deg == UINT_E_MAX) return 0;
    std::cout << "Deg: " << get_bucket(deg) << std::endl;
    return get_bucket(deg);
  });

  auto D_capped = sequence<bucket_t>::from_function(num_entries, [&](size_t i) -> bucket_t {
    return cliques->get_count(i);
  });

std::cout << "Start b" << std::endl; fflush(stdout);
  auto b = buckets<sequence<bucket_t>, iden_t, bucket_t>(num_entries, D, increasing, num_buckets);
std::cout << "End b" << std::endl; fflush(stdout);
  auto per_processor_counts = sequence<double>(0);
  
  list_buffer count_idxs(num_entries, efficient);

  char* still_active = (char*) calloc(num_entries, sizeof(char));
  size_t max_deg = induced_hybrid::get_max_deg(G);

  timer t_extract;
  timer t_count;
  timer t_update;
  timer t_x;
  timer t_update_d;

  size_t rounds = 0;
  size_t finished = 0;
  bucket_t cur_bkt = 0;
  bucket_t max_bkt = 0;
  double max_density = 0;
  bool use_max_density = false;
  size_t iter = 0;
  bucket_t prev_bkt = 0;

  size_t cur_inner_rounds = 0;
  size_t max_inner_rounds = log(num_entries) / log(1.0 + eps);

  while (finished < num_entries) {
    t_extract.start();
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers; //vertexSubset(num_entries, bkt.identifiers);
    auto active_size = active.size();
    cur_bkt = bkt.id;
    t_extract.stop();

    auto get_active = [&](size_t j) -> unsigned __int128 { 
      return active[j]; }; //active.vtx(j); };

    if (active_size == 0) continue;

    finished += active_size;

    if (cur_bkt == UINT_E_MAX) continue;

    if (prev_bkt != cur_bkt) {
      cur_inner_rounds = 0;
    }

     // Check if we hit the threshold for inner peeling rounds.
    if (cur_inner_rounds == max_inner_rounds) {
      // new re-insertions will go to at least bucket k (one greater than before).
      cur_bkt++;
      cur_inner_rounds = 0;
    }
    uintE lower_bound = ceil(pow((1 + delta), cur_bkt-1));

    if (inline_hierarchy && prev_bkt != cur_bkt && cur_bkt != 0) {
      connect_while_peeling.init(cur_bkt);
    }

    max_bkt = std::max(cur_bkt, max_bkt);
    if (cur_bkt == 0) { //|| finished == num_entries
      parallel_for (0, active_size, [&] (size_t j) {
        auto index = get_active(j);
        still_active[index] = 2;
        cliques->set_count(index, UINT_E_MAX);
        D[index] = std::max((bucket_t) 0, (bucket_t) cur_bkt);
      }, 2048);
      continue;
    }

    iter++;

    size_t granularity = (cur_bkt * active_size < 10000) ? 1024 : 1;

    size_t filter_size = 0;

      auto update_changed = [&](sequence<double>& ppc, size_t i, uintE v){

      };
    t_count.start();
    count_idxs.resize(active_size, k, r, cur_bkt);
    filter_size = cliqueUpdate(G, DG, r, k, max_deg, true, get_active, active_size, 
      granularity, still_active, rank, per_processor_counts,
      false, update_changed, cliques, num_entries, count_idxs, t_x, inverse_rank, relabel,
      t_update_d, D, false, inline_hierarchy, connect_while_peeling, cur_bkt);

    t_count.stop();

    t_update.start();
    auto num_count_idxs = count_idxs.num_entries();
    auto D_filter = sequence<bucket_t>(num_count_idxs);
    
    parallel_for(0, num_count_idxs, [&](size_t i){
      auto v = count_idxs.get_v(i);
        if (v == UINT_E_MAX) {
          count_idxs.void_v(i, v);
          D_filter[i] = num_entries + 1;
        } 
        else {
          bucket_t deg = D[v];
          auto val = cliques->get_count(v);
          bucket_t new_deg = std::max((bucket_t) val, (bucket_t) lower_bound);
          D_capped[v] = new_deg;
          uintE new_bkt = std::max((bucket_t) get_bucket(new_deg),(bucket_t) cur_bkt);
          if (deg != new_bkt) {
            D[v] = new_bkt;
            D_filter[i] = b.get_bucket(deg, new_bkt);
          } else {
            count_idxs.void_v(i, v);
            D_filter[i] = num_entries + 1;
          }
        }

        if (v != UINT_E_MAX) {
          gbbs::CAS(&(still_active[v]), char{3}, char{0});
          gbbs::CAS(&(still_active[v]), char{4}, char{1});
        }
    });
    //std::cout << "FLAG 2" << std::endl; fflush(stdout);
    auto apply_f = [&](size_t i) -> std::optional<std::tuple<uintE, bucket_t>> {
      auto v = count_idxs.get_v(i);
      if (v != UINT_E_MAX) {
      bucket_t bucket = D_filter[i];
        if (still_active[v] != 2 && still_active[v] != 1 && still_active[v] != 4) return wrap(v, bucket);
      }
      return std::nullopt;
    };

    
    b.update_buckets(apply_f, num_count_idxs);
    count_idxs.reset();
    t_update.stop();

    rounds++;
    prev_bkt = cur_bkt;
    cur_inner_rounds++;
  
  }

  t_extract.next("### Peel Extract time: ");
  t_count.next("### Peel Count time: ");
  t_update.next("### Peel Update time: ");
  t_x.next("Inner counting: ");
  t_update_d.next("Update d time: ");

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  std::cout.precision(17);
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "clique core: " << max_bkt << std::endl;
  if (use_max_density) std::cout << "max density: " << max_density << std::endl;

  free(still_active);

  parallel_for(0, num_entries, [&] (size_t i) {
    if (use_pow) {  // use 2^{peeled_bkt} as the coreness estimate
      D[i] = (D[i] == 0) ? 0 : 1 << D[i];
    } else {
      D[i] = D_capped[i];  // use capped induced degree when peeled as the coreness estimate
    }
  });

  return D;
}

} // end namespace gbbs