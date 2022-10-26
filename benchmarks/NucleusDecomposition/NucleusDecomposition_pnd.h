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


template<class Graph>
void intersectionPND(Graph& G, uintE v, uintE u, std::vector<uintE>& intersection){
  auto vert_v = G.get_vertex(v);
  auto vert_u = G.get_vertex(u);
  uintE idx_v = 0;
  uintE idx_u = 0;
  auto iter_v = vert_v.out_neighbors().get_iter();
  auto iter_u = vert_u.out_neighbors().get_iter();
  auto deg_v = vert_v.out_degree();
  auto deg_u = vert_u.out_degree();
  while(idx_v < deg_v && idx_u < deg_u) {
    uintE v_nbhr = std::get<0>(iter_v.cur());
    uintE u_nbhr = std::get<0>(iter_u.cur());
    if (v_nbhr < u_nbhr) {
      idx_v++;
      if (iter_v.has_next()) iter_v.next();
    }
    else if (u_nbhr < v_nbhr){
      idx_u++;
      if (iter_u.has_next()) iter_u.next();
    } 
    else {
      intersection.push_back(v_nbhr);
      idx_u++; idx_v++;
      if (iter_v.has_next()) iter_v.next();
      if (iter_u.has_next()) iter_u.next();
    }
  }
}

template<class Graph>
void intersectionPND(Graph& G, uintE v, uintE u, uintE w, std::vector<uintE>& intersection){
  auto vert_v = G.get_vertex(v);
  auto vert_u = G.get_vertex(u);
  auto vert_w = G.get_vertex(w);
  uintE idx_v = 0;
  uintE idx_u = 0;
  uintE idx_w = 0;
  auto iter_v = vert_v.out_neighbors().get_iter();
  auto iter_u = vert_u.out_neighbors().get_iter();
  auto iter_w = vert_w.out_neighbors().get_iter();
  auto deg_v = vert_v.out_degree();
  auto deg_u = vert_u.out_degree();
  auto deg_w = vert_w.out_degree();
  while(idx_v < deg_v && idx_u < deg_u && idx_w < deg_w) {
    uintE v_nbhr = std::get<0>(iter_v.cur());
    uintE u_nbhr = std::get<0>(iter_u.cur());
    uintE w_nbhr = std::get<0>(iter_w.cur());
    if (u_nbhr == v_nbhr && u_nbhr == w_nbhr) {
      intersection.push_back(v_nbhr);
      idx_u++; idx_v++; idx_w++;
      if (iter_v.has_next()) iter_v.next();
      if (iter_u.has_next()) iter_u.next();
      if (iter_w.has_next()) iter_w.next();
    } else {
      auto max_nbhr = std::max(std::max(v_nbhr, u_nbhr), w_nbhr);
      if (v_nbhr < max_nbhr) {
        idx_v++;
        if (iter_v.has_next()) iter_v.next();
      }
      if (u_nbhr < max_nbhr) {
        idx_u++;
        if (iter_u.has_next()) iter_u.next();
      }
      if (w_nbhr < max_nbhr) {
        idx_w++;
        if (iter_w.has_next()) iter_w.next();
      }
    }
  }
}


  template<class Graph, class T>
  inline size_t CountCliquesNucPND(Graph& DG, size_t k, size_t r, size_t max_deg, T* table) {
    k--; r--;

    if (k == 2) {
      auto tots = sequence<size_t>(DG.n, size_t{0});

      auto base_f = [&](uintE v1, uintE v2, uintE v3) {
        table->insert_twothree(v1, v2, v3, r, k);
      };
      parallel_for(0, DG.n, [&](size_t i) {
        auto vert_i = DG.get_vertex(i);
        auto iter_i = vert_i.out_neighbors().get_iter();
        for (std::size_t j = 0; j < vert_i.out_degree(); j++) {
          auto x = std::get<0>(iter_i.cur());
          if (iter_i.has_next()) iter_i.next();
          std::vector<uintE> inter;
          intersectionPND(DG, i, x, inter);
          tots[i] += inter.size();
          for (std::size_t l = 0; l < inter.size(); l++) {
            table->insert_twothree(i, x, inter[l], r, k);
          }
        }
      });
      return parlay::reduce(tots);
    }
    // Now k must be 3 (we'll count 4-cliques)
    assert(k == 3);
    auto tots = sequence<size_t>(DG.n, size_t{0});

    auto base_f = [&](sequence<uintE>& base){
      table->insert(base, r, k);
    };
    parallel_for(0, DG.n, [&](size_t i) {
      auto base = sequence<uintE>(k + 1);
      base[0] = i;
      auto vert_i = DG.get_vertex(i);
      auto iter_i = vert_i.out_neighbors().get_iter();
      for (std::size_t j = 0; j < vert_i.out_degree(); j++) {
        auto x = std::get<0>(iter_i.cur());
        if (iter_i.has_next()) iter_i.next();
        base[1] = x;
        std::vector<uintE> inter;
        intersectionPND(DG, i, x, inter);
        // Now any two vert in inter that are in each other's nbhr list...
        for (std::size_t l = 0; l < inter.size(); l++) {
          for (std::size_t p = l + 1; p < inter.size(); p++) {
            auto v1 = inter[l];
            auto v2 = inter[p];
            if (is_edge(DG, v1, v2)) {
              base[2] = v1; base[3] = v2; base_f(base); tots[i]++;
            } else if (is_edge(DG, v2, v1)) {
              base[2] = v1; base[3] = v2; base_f(base); tots[i]++;
            }
          }
        }
      }
    });

    return parlay::reduce(tots);
  }

// This should only be used for (2, 3) and (3, 4) nucleus decomposition
template <class Graph, class Graph2, class F, class I, class T>
inline size_t cliqueUpdatePND(Graph& G, Graph2& DG, size_t r, 
size_t k, size_t max_deg, bool label, F get_active, size_t active_size,
  size_t granularity, char* still_active, sequence<uintE> &rank, 
  sequence<double>& per_processor_counts, 
  bool do_update_changed, I& update_changed,
  T* cliques, size_t n, list_buffer& count_idxs, timer& t1,
  sequence<uintE>& inverse_rank, bool relabel, timer& t_update_d,
  bool use_ppc = true) {
    using W = typename Graph::weight_type;

    // Mark every vertex in the active set
    parallel_for (0, active_size, [&] (size_t j) {
      auto index = get_active(j); //cliques->find_index(get_active(j));
      still_active[index] = 1;
    }, 2048);
    
    auto is_active = [&](size_t index) {
    return still_active[index] == 1 || still_active[index] == 4;
  };
  auto is_inactive = [&](size_t index) {
    return still_active[index] == 2;
  };

  t1.start();
  if (k == 2 && r == 1) { // This is (2, 3)
    auto update_d_twothree = [&](uintE v1, uintE v2, uintE v3){
    cliques->extract_indices_twothree(v1, v2, v3, is_active, is_inactive,
      [&](std::size_t index, double val){
        if (use_ppc) {
          double ct = gbbs::fetch_and_add(&(per_processor_counts[index]), val);
          if (ct == 0 && val != 0) count_idxs.add(index);
        } else {
          if (!is_inactive(index)) {
          cliques->update_count_atomic(index, gbbs::uintE{std::round(val)});
          if (gbbs::CAS(&(still_active[index]), char{0}, char{3}) || gbbs::CAS(&(still_active[index]), char{1}, char{4}))
            count_idxs.add(index);
          }
        }
      }, r, k);
    };

      parallel_for(0, active_size, [&](size_t i){
        auto x = get_active(i);
        std::tuple<uintE, uintE> v1v2 = cliques->extract_clique_two(x, k);
        uintE u = relabel ? inverse_rank[std::get<0>(v1v2)] : std::get<0>(v1v2);
        uintE v = relabel ? inverse_rank[std::get<1>(v1v2)] : std::get<1>(v1v2);
        std::vector<uintE> u_v_intersection;
        intersectionPND(G, u, v, u_v_intersection);
        for (std::size_t p = 0; p < u_v_intersection.size(); p++) {
          uintE v3 = relabel ? rank[u_v_intersection[p]] : u_v_intersection[p];
          update_d_twothree(std::get<0>(v1v2), std::get<1>(v1v2), v3);
        }
      });
  } else if (k == 3 && r == 2) { // (3, 4))
  auto update_d_threefour = [&](uintE v1, uintE v2, uintE v3, uintE v4){
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

      parallel_for(0, active_size, [&](size_t i) {
        auto x = get_active(i);
        std::tuple<uintE, uintE, uintE> v1v2v3 = cliques->extract_clique_three(x, k);
        uintE u = relabel ? inverse_rank[std::get<0>(v1v2v3)] : std::get<0>(v1v2v3);
        uintE v = relabel ? inverse_rank[std::get<1>(v1v2v3)] : std::get<1>(v1v2v3);
        uintE w = relabel ? inverse_rank[std::get<2>(v1v2v3)] : std::get<2>(v1v2v3);
        std::vector<uintE> u_v_w_intersection;
        intersectionPND(G, u, v, w, u_v_w_intersection);
        for (std::size_t p = 0; p < u_v_w_intersection.size(); p++) {
          uintE actual_ngh = relabel ? rank[u_v_w_intersection[p]] : u_v_w_intersection[p];
          update_d_threefour(std::get<0>(v1v2v3), std::get<1>(v1v2v3), std::get<2>(v1v2v3), actual_ngh);
        }
      },1, true);
  }
  else {std::cout << "ERROR" << std::endl; assert(false); exit(0); } // should never happen
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

} // namespace gbbs