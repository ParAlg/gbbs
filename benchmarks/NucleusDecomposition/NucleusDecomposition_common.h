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

#include "multitable.h"
#include "twotable.h"
#include "twotable_nosearch.h"
#include "onetable.h"
#include "commontable.h"
#include "multitable_nosearch.h"
#include "list_buffer.h"

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
/*
  template <class CU, class Graph, class F, class FF, class Table>
  void CompressOut(CU& compress_utils, Graph& GA, F& g_to_dg, FF& dg_to_g,
    char* still_active, timer& ct, Table* cliques) {
      using bucket_t = uintE;
    using W = typename Graph::weight_type;
      ct.start();
    auto del_edges = compress_utils.del_edges;
    auto em = compress_utils.em;
    auto actual_degree = compress_utils.actual_degree;
      // compact
      // map over both endpoints, update counts using histogram
      // this is really a uintE seq, but edge_t >= uintE, and this way we can
      // re-use the same histogram structure.
      auto decr_seq = pbbs::sequence<uintE>(2*del_edges.size);
      parallel_for(0, del_edges.size, [&] (size_t i) {
        size_t fst = 2*i; size_t snd = fst+1;
        uintE id = del_edges.A[i];
        std::tuple<uintE, uintE> uv = cliques->extract_clique_two(id, 2);
        decr_seq[fst] = dg_to_g(std::get<0>(uv));
        decr_seq[snd] = dg_to_g(std::get<1>(uv));
      });

      // returns only those vertices that have enough degree lost to warrant
      // packing them out. Again note that edge_t >= uintE
      auto apply_vtx_f = [&](const std::tuple<uintE, bucket_t>& p)
        -> const std::optional<std::tuple<uintE, bucket_t> > {
        uintE id = std::get<0>(p);
        uintE degree_lost = std::get<1>(p);
        actual_degree[id] -= degree_lost;
        // compare with GA.V[id]. this is the current space used for this vtx.
        return std::nullopt;
      };

      auto vs = vertexSubset(GA.n, decr_seq.size(), decr_seq.begin());
      auto cond_f = [&] (const uintE& u) { return true; };
      nghCount(GA, vs, cond_f, apply_vtx_f, em);

      auto all_vertices = pbbs::delayed_seq<uintE>(GA.n, [&] (size_t i) { return i; });
      auto to_pack_seq = pbbs::filter(all_vertices, [&] (uintE u) {
        return 4*actual_degree[u] >= GA.get_vertex(u).getOutDegree();
      });
      auto to_pack = vertexSubset(GA.n, std::move(to_pack_seq));

      auto pack_predicate = [&](const uintE& u, const uintE& ngh, const W& wgh) {
        // return true iff edge is still alive
        auto index = cliques->extract_indices_two(g_to_dg(u), g_to_dg(ngh));
        return still_active[index] == 0;
      };
      edgeMapFilter(GA, to_pack, pack_predicate, pack_edges | no_output);

      ct.stop();
  }

  class compress_utils {
    public:
    using bucket_t = uintE;
    pbbslib::dyn_arr<uintE> del_edges;
    hist_table<uintE, bucket_t> em;
    pbbs::sequence<uintE> actual_degree;
    size_t n;

    compress_utils(){}

    template <class Graph>
    compress_utils(Graph& GA){
      n = GA.n;
      del_edges = pbbslib::dyn_arr<uintE>(6 * GA.n);
      std::tuple<uintE, bucket_t> histogram_empty = std::make_tuple(std::numeric_limits<uintE>::max(), 0);
      em = hist_table<uintE, bucket_t>(histogram_empty, GA.m/50);
      actual_degree = pbbs::sequence<uintE>(GA.n, [&] (size_t i) {
        return GA.get_vertex(i).getOutDegree();
      });
    }

    template <class A>
    bool update_del_edges(A& active) {
      del_edges.copyIn(active, active.size());
      return del_edges.size > 2*n;
    }
  };*/

  struct IntersectSpace {
  // Label each vertex for fast intersect for first recursive level
  uintE* labels = nullptr;

  IntersectSpace () {}

  void alloc(size_t n) {
    if (labels == nullptr) {
      labels = (uintE*) calloc(n, sizeof(uintE));
      assert(labels[0] == 0);
    }
  }
  };

  template <class Graph, class T>
  inline size_t CountCliquesNuc(Graph& DG, size_t k, size_t r, size_t max_deg, T* table) {
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

    return parlay::reduce(tots);
  }



template <class Graph>
bool is_edge(Graph& DG, uintE v, uintE u) {
  using W = typename Graph::weight_type;
  bool is = false;
  auto map_f = [&] (const uintE& src, const uintE& vv, const W& wgh) {
    if (vv == u) is = true;
    };
    DG.get_vertex(v).out_neighbors().map(map_f, false);
    return is;
}

/*template <
    template <class W> class vertex, class W, 
    typename std::enable_if<
        std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value, int>::type = 0>
void intersectionPND(
    symmetric_graph<vertex, W>& GA, uintE v, uintE u, std::vector<uintE>& intersection) {  // -> decltype(GA)
    std::cout << "IntersectionPND not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
    }

template <
    template <class W> class vertex, class W, 
    typename std::enable_if<std::is_same<vertex<W>, symmetric_vertex<W>>::value,
                            int>::type = 0>
void intersectionPND(symmetric_graph<vertex, W>& G, */
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

/*
template <
    template <class W> class vertex, class W, 
    typename std::enable_if<
        std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value, int>::type = 0>
void intersectionPND(
    symmetric_graph<vertex, W>& GA, uintE v, uintE u, uintE w, std::vector<uintE>& intersection) {  // -> decltype(GA)
    std::cout << "IntersectionPND not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
    }

template <
    template <class W> class vertex, class W, 
    typename std::enable_if<std::is_same<vertex<W>, symmetric_vertex<W>>::value,
                            int>::type = 0>
void intersectionPND(symmetric_graph<vertex, W>& G, */
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
/*
template <
    template <class W> class vertex, class W, class T,
    typename std::enable_if<
        std::is_same<vertex<W>, asymmetric_vertex<W>>::value, int>::type = 0>
inline size_t CountCliquesNucPND(asymmetric_graph<vertex, W>& G, size_t k, size_t r, size_t max_deg, T* table){
  std::cout << "Filter graph not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return 0;
}

template <
    template <class W> class vertex, class W, class T,
    typename std::enable_if<
        std::is_same<vertex<W>, cav_bytepd_amortized<W>>::value, int>::type = 0>
inline size_t CountCliquesNucPND(asymmetric_graph<vertex, W>& G, size_t k, size_t r, size_t max_deg, T* table){
  std::cout << "Filter graph not implemented for directed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return 0;
}

template <
    template <class W> class vertex, class W, class T,
    typename std::enable_if<
        std::is_same<vertex<W>, csv_bytepd_amortized<W>>::value, int>::type = 0>
inline size_t CountCliquesNucPND(
    symmetric_graph<vertex, W>& GA, size_t k, size_t r, size_t max_deg, T* table) {  // -> decltype(GA)
    std::cout << "CountCliquesNucPND not implemented for compressed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return 0;
    }

    template <
    template <class W> class vertex, class W, class T,
    typename std::enable_if<
        std::is_same<vertex<W>, csv_byte<W>>::value, int>::type = 0>
inline size_t CountCliquesNucPND(
    symmetric_graph<vertex, W>& GA, size_t k, size_t r, size_t max_deg, T* table) {  // -> decltype(GA)
    std::cout << "CountCliquesNucPND not implemented for compressed graphs" << std::endl;
  assert(false);  // Not implemented for directed graphs
  return 0;
    }

template <
    template <class W> class vertex, class W, class T,
    typename std::enable_if<std::is_same<vertex<W>, symmetric_vertex<W>>::value,
                            int>::type = 0>
  inline size_t CountCliquesNucPND( symmetric_graph<vertex, W>& DG, */
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
    auto update_d = [&](uintE* base){
    cliques->extract_indices(base, is_active, is_inactive, [&](std::size_t index, double val){
      if (use_ppc) {
        double ct = gbbs::fetch_and_add(&(per_processor_counts[index]), val);
        if (ct == 0 && val != 0) {
          count_idxs.add(index);
        }
      } else {
        if (!is_inactive(index)) {
        cliques->update_count_atomic(index, gbbs::uintE{std::round(val)});
        if (gbbs::CAS(&(still_active[index]), char{0}, char{3}) || gbbs::CAS(&(still_active[index]), char{1}, char{4}))
          count_idxs.add(index);
        }
      }
    }, r, k);
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
              if (!is_inactive(index)) {
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

template <class Graph, class Graph2, class F, class I, class T>
inline size_t cliqueUpdate(Graph& G, Graph2& DG, size_t r, 
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
    return still_active[index] == 1  || still_active[index] == 4;
  };
  auto is_inactive = [&](size_t index) {
    return still_active[index] == 2;
  };

  auto update_d = [&](uintE* base){
    cliques->extract_indices(base, is_active, is_inactive, [&](std::size_t index, double val){
      if (use_ppc) {
      double ct = gbbs::fetch_and_add(&(per_processor_counts[index]), val);
      if (ct == 0 && val != 0) {
        count_idxs.add(index);
      }
      } else {
        if (!is_inactive(index)) {
        cliques->update_count_atomic(index, gbbs::uintE{std::round(val)});
        if (gbbs::CAS(&(still_active[index]), char{0}, char{3}) || gbbs::CAS(&(still_active[index]), char{1}, char{4}))
          count_idxs.add(index);
        }
      }
    }, r, k);
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
        auto v_v = G.get_vertex(v);
        auto process_f = [&](uintE a, uintE b, uintE intersect_w){
          uintE v3 = relabel ? rank[intersect_w] : intersect_w;
          update_d_twothree(std::get<0>(v1v2), std::get<1>(v1v2), v3);
        };
        auto their_neighbors = v_v.out_neighbors();
        G.get_vertex(u).out_neighbors().intersect_f_par(&their_neighbors, process_f);

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
              if (!is_inactive(index)) {
              cliques->update_count_atomic(index, gbbs::uintE{std::round(val)});
        if (gbbs::CAS(&(still_active[index]), char{0}, char{3}) || gbbs::CAS(&(still_active[index]), char{1}, char{4}))
          count_idxs.add(index);
              }
            }
        }, r, k);
      };

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
            update_d(base);
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
    NKCliqueDir_fast_hybrid_rec(DG, 1, k-r, induced, update_d, base);
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

template <typename bucket_t, class Graph, class Graph2, class T>
sequence<bucket_t> Peel(Graph& G, Graph2& DG, size_t r, size_t k, 
  T* cliques, sequence<uintE> &rank, size_t fake_efficient, bool relabel, 
  bool use_compress,
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
  double max_density = 0;
  bool use_max_density = false;
  size_t iter = 0;

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

    //if (cur_bkt >= std::numeric_limits<bucket_t>::max() - 1) continue;

    max_bkt = std::max(cur_bkt, max_bkt);
    if (cur_bkt == 0 || finished == num_entries) {
      parallel_for (0, active_size, [&] (size_t j) {
        auto index = get_active(j);
        still_active[index] = 2;
        cliques->set_count(index, UINT_E_MAX);
      }, 2048);
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
          if (std::round(ppc[v]) > clique_count){
            std::cout << "PPC: " << std::round(ppc[v]) << ", count: " << clique_count << ", v: " << v << std::endl;
            fflush(stdout);
            exit(0);
          }
          assert(std::round(ppc[v]) <= clique_count);
          
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
    if (fake_efficient == 3) {
      filter_size = cliqueUpdatePND(G, DG, r, k, max_deg, true, get_active, active_size, 
     granularity, still_active, rank, per_processor_counts,
      true, update_changed, cliques, num_entries, count_idxs, t_x, inverse_rank, relabel,
      t_update_d);
    } else {
     filter_size = cliqueUpdate(G, DG, r, k, max_deg, true, get_active, active_size, 
     granularity, still_active, rank, per_processor_counts,
      true, update_changed, cliques, num_entries, count_idxs, t_x, inverse_rank, relabel,
      t_update_d);
    }
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


template <typename bucket_t, class Graph, class Graph2, class T>
sequence<bucket_t> Peel_space_efficient(Graph& G, Graph2& DG, size_t r, size_t k, 
  T* cliques, sequence<uintE> &rank, size_t fake_efficient, bool relabel, 
  bool use_compress,
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
  std::cout << "created 1 " << std::endl; fflush(stdout);

  auto b = make_vertex_custom_buckets<bucket_t>(num_entries, D, increasing, num_buckets);
  std::cout << "created 3 " << std::endl; fflush(stdout);

  auto per_processor_counts = sequence<double>(0);
  std::cout << "created 4 " << std::endl; fflush(stdout);
  
  list_buffer count_idxs(num_entries, efficient);
  std::cout << "created 5 " << std::endl; fflush(stdout);

  char* still_active = (char*) calloc(num_entries, sizeof(char));
  std::cout << "created 6 " << std::endl; fflush(stdout);
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

    max_bkt = std::max(cur_bkt, max_bkt);
    if (cur_bkt == 0 || finished == num_entries) {
      parallel_for (0, active_size, [&] (size_t j) {
        auto index = get_active(j);
        still_active[index] = 2;
        cliques->set_count(index, UINT_E_MAX);
      }, 2048);
      continue;
    }

    std::cout << "k = " << cur_bkt << " iter = " << iter << " #edges = " << active_size << std::endl;
    std::cout << "Finished: " << finished << ", num_entries: " << num_entries << std::endl;
    iter++;

    size_t granularity = (cur_bkt * active_size < 10000) ? 1024 : 1;

    size_t filter_size = 0;

      auto update_changed = [&](sequence<double>& ppc, size_t i, uintE v){

      };
    t_count.start();
    count_idxs.resize(active_size, k, r, cur_bkt);
    if (fake_efficient == 3) {
      filter_size = cliqueUpdatePND(G, DG, r, k, max_deg, true, get_active, active_size, 
     granularity, still_active, rank, per_processor_counts,
      false, update_changed, cliques, num_entries, count_idxs, t_x, inverse_rank, relabel,
      t_update_d, false);
    } else {
     filter_size = cliqueUpdate(G, DG, r, k, max_deg, true, get_active, active_size, 
     granularity, still_active, rank, per_processor_counts,
      false, update_changed, cliques, num_entries, count_idxs, t_x, inverse_rank, relabel,
      t_update_d, false);
    }
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
        if (v >= D.size()) {std::cout << "v: " << v << ", size: " << D.size() << std::endl; fflush(stdout);}
        assert(v < D.size());
      bucket_t bucket = D_filter[i];
      assert(bucket != num_entries + 1);
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

//*************************************************************VERIFICATION CODE**********


template <typename bucket_t, class Graph, class Graph2, class T, class T2>
sequence<bucket_t> Peel_verify(Graph& G, Graph2& DG, size_t r, size_t k, 
  T* cliques, T2* cliques2, sequence<uintE> &rank, size_t num_rounds_verify=2,
  size_t num_buckets=16) {
  bool relabel = false;
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
  size_t num_entries2 = cliques2->return_total();

  auto D = sequence<bucket_t>::from_function(num_entries, [&](size_t i) -> bucket_t { 
    return cliques->get_count(i);
  });
  auto D2 = sequence<bucket_t>::from_function(num_entries2, [&](size_t i) -> bucket_t { 
    return cliques2->get_count(i);
  });

  auto D_filter = sequence<std::tuple<uintE, bucket_t>>(num_entries);
  auto D_filter2 = sequence<std::tuple<uintE, bucket_t>>(num_entries2);

  auto b = make_vertex_custom_buckets<bucket_t>(num_entries, D, increasing, num_buckets);
  auto b2 = make_vertex_custom_buckets<bucket_t>(num_entries2, D2, increasing, num_buckets);

  auto per_processor_counts = sequence<double>(num_entries , static_cast<double>(0));
  auto per_processor_counts2 = sequence<double>(num_entries2 , static_cast<double>(0));
  
  list_buffer count_idxs(num_entries);
  list_buffer count_idxs2(num_entries2);

  char* still_active = (char*) calloc(num_entries, sizeof(char));
  char* still_active2 = (char*) calloc(num_entries2, sizeof(char));

  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active?

  timer t_extract;
  timer t_count;
  timer t_update;
  timer t_x;
  timer t_update_d;

  size_t rounds = 0;
  size_t finished = 0;
  bucket_t cur_bkt = 0;
  bucket_t cur_bkt2 = 0;
  bucket_t max_bkt = 0;
  size_t iter = 0;

  while (finished != num_entries) {
    t_extract.start();
    // Retrieve next bucket
    auto bkt = b.next_bucket();
    auto active = bkt.identifiers; //vertexSubset(num_entries, bkt.identifiers);
    auto active_size = active.size();
    cur_bkt = bkt.id;

    auto bkt2 = b2.next_bucket();
    auto active2 = bkt2.identifiers;
    auto active_size2 = active2.size();
    cur_bkt2 = bkt2.id;
    t_extract.stop();

    assert(cur_bkt == cur_bkt2);
    assert(active_size == active_size2);

    auto get_active = [&](size_t j) -> unsigned __int128 {  return active[j]; };
    auto get_active2 = [&](size_t j) -> unsigned __int128 {  return active2[j]; };

    // Assume active_size == active_size2
    if (active_size == 0) continue;
    finished += active_size;

    // Assume cur_bkt == cur_bkt2
    if (cur_bkt == UINT_E_MAX) continue;
    max_bkt = std::max(cur_bkt, max_bkt);
  
    if (cur_bkt == 0 || finished == num_entries) {
      parallel_for (0, active_size, [&] (size_t j) {
        auto index = get_active(j);
        still_active[index] = 2;
        cliques->set_count(index, UINT_E_MAX);

        auto index2 = get_active2(j);
        still_active2[index2] = 2;
        cliques2->set_count(index2, UINT_E_MAX);
      }, 2048);
      continue;
    }

    std::cout << "k = " << cur_bkt << " iter = " << iter << " #edges = " << active_size << std::endl;
    std::cout << "k2 = " << cur_bkt2 << " iter = " << iter << " #edges2 = " << active_size2 << std::endl;
    if (iter >= num_rounds_verify) exit(0);
    iter++;

    // Verify active set
    /*for (size_t i = 0; i < active_size; i++) {
      auto vtx = get_active(i);
      // Check that vtx exists in cliques2
      sequence<uintE> base(k + 1);
      // Vertices will be in inclusive k, ..., k - r + 1, 0
      cliques->extract_clique(vtx, base, G, k);
      sequence<uintE> actual_base(r + 1);
      for (size_t j = 0; j <= r; j++) {
        auto base_idx = k - j;
        if (j == r) base_idx = 0;
        actual_base[j] = base[base_idx];
      }

      // Check that extracting a clique from cliques and turning it back to an index works
      auto check_vtx = cliques->extract_indices_check(actual_base, r);
      assert(check_vtx == vtx);

      auto vtx2 = cliques2->extract_indices_check(actual_base, r);
      sequence<uintE> base2(k + 1);
      cliques2->extract_clique(vtx2, base2, G, k);
      sequence<uintE> actual_base2(r + 1);
      for (size_t j = 0; j <= r; j++) {
        auto base_idx = k - j;
        if (j == r) base_idx = 0;
        actual_base2[j] = base2[base_idx];
      }

      // Check that extracting a clique from cliques2 and turning it back to an index works
      auto check_vtx2 = cliques2->extract_indices_check(actual_base2, r);
      assert(check_vtx2 == vtx2);

      // Check that the clique counts match up
      assert(cliques->get_count(vtx) == cliques2->get_count(vtx2));
      // On first round, degree should be 1
      if (iter == 1) {
        assert(cliques->get_count(vtx) == 1);
      }

      // Check that cliques and cliques2 are thinking about the same clique
      pbbslib::sample_sort_inplace (actual_base.slice(), [&](const uintE& u, const uintE&  v) {
          return u < v;
      });
      pbbslib::sample_sort_inplace (actual_base2.slice(), [&](const uintE& u, const uintE&  v) {
          return u < v;
      });
      for (size_t j = 0; j < r + 1; j++) {
        assert(actual_base[j] == actual_base2[j]);
      }
    }*/

    size_t granularity = (cur_bkt * active_size < 10000) ? 1024 : 1;

    auto update_changed = [&](sequence<double>& ppc, size_t i, uintE v){
      if (v == UINT_E_MAX) {
        D_filter[i] = std::make_tuple(num_entries + 1, 0);
        return;
      }
      assert(ppc[v] != 0);
      if (ppc[v] == 0) D_filter[i] = std::make_tuple(num_entries + 1, 0);
      else {
        bucket_t deg = D[v];
        assert(deg > cur_bkt);
        auto val = cliques->get_count(v) - std::round(ppc[v]);
        if (v == 2297107) {
          std::cout << "Val: " << val << ", prev count: " << cliques->get_count(v) << std::endl;
        }
        cliques->set_count(v, val);
        if (v == 2297107) {
          std::cout << "new count: " << cliques->get_count(v) << std::endl; fflush(stdout);
        }
        //auto val = cliques->update_count(v, (size_t) ppc[v]);
        if (deg > cur_bkt) {
            bucket_t new_deg = std::max((bucket_t) val, (bucket_t) cur_bkt);
            D[v] = new_deg;
            D_filter[i] = std::make_tuple(v, b.get_bucket(deg, new_deg));
        } else D_filter[i] = std::make_tuple(num_entries + 1, 0);
      }
      //ppc[v] = 0;
    };
    auto update_changed2 = [&](sequence<double>& ppc, size_t i, uintE v){
      if (v == UINT_E_MAX) {
        D_filter2[i] = std::make_tuple(num_entries2 + 1, 0);
        return;
      }
      assert(ppc[v] != 0);
      if (ppc[v] == 0) D_filter2[i] = std::make_tuple(num_entries2 + 1, 0);
      else {
        bucket_t deg = D2[v];
        assert(deg > cur_bkt2);
        auto val = cliques2->get_count(v) - std::round(ppc[v]);
        if (v == 1509298) {
          std::cout << "Val: " << val << ", prev count: " << cliques2->get_count(v) << std::endl;
        }
        cliques2->set_count(v, val);
        if (v == 1509298) {
          std::cout << "new count: " << cliques2->get_count(v) << std::endl; fflush(stdout);
        }
        //auto val = cliques2->update_count(v, (size_t) ppc[v]);
        if (deg > cur_bkt2) {
            bucket_t new_deg = std::max((bucket_t) val, (bucket_t) cur_bkt2);
            D2[v] = new_deg;
            D_filter2[i] = std::make_tuple(v, b2.get_bucket(deg, new_deg));
        } else D_filter2[i] = std::make_tuple(num_entries2 + 1, 0);
      }
      //ppc[v] = 0;
    };
    t_count.start();
    size_t filter_size = cliqueUpdate(G, DG, r, k, max_deg, true, get_active, active_size, 
      granularity, still_active, rank, per_processor_counts,
      true, update_changed, cliques, num_entries, count_idxs, t_x, inverse_rank, relabel,
      t_update_d);
    
    size_t filter_size2 = cliqueUpdate(G, DG, r, k, max_deg, true, get_active2, active_size2, 
      granularity, still_active2, rank, per_processor_counts2,
      true, update_changed2, cliques2, num_entries2, count_idxs2, t_x, inverse_rank, relabel,
      t_update_d);
    t_count.stop();

    std::cout << "Finished verifying active set; starting verifying updates" << std::endl;
    fflush(stdout);

    // Verify updated vertices from D_filter
    for (size_t i = 0; i < filter_size; i++) {
      auto vtx = std::get<0>(D_filter[i]);
      if (vtx == num_entries + 1) continue;
      // Check that vtx exists in cliques2
      sequence<uintE> base(k + 1);
      // Vertices will be in inclusive k, ..., k - r + 1, 0
      cliques->extract_clique(vtx, base.data(), G, k);
      sequence<uintE> actual_base(r + 1);
      for (size_t j = 0; j <= r; j++) {
        auto base_idx = k - j;
        if (j == r) base_idx = 0;
        actual_base[j] = base[base_idx];
      }

      // Check that extracting a clique from cliques and turning it back to an index works
      auto check_vtx = cliques->extract_indices_check(actual_base.data(), r);
      assert(check_vtx == vtx);

      auto vtx2 = cliques2->extract_indices_check(actual_base.data(), r);
      sequence<uintE> base2(k + 1);
      cliques2->extract_clique(vtx2, base2.data(), G, k);
      sequence<uintE> actual_base2(r + 1);
      for (size_t j = 0; j <= r; j++) {
        auto base_idx = k - j;
        if (j == r) base_idx = 0;
        actual_base2[j] = base2[base_idx];
      }

      // Check that extracting a clique from cliques2 and turning it back to an index works
      auto check_vtx2 = cliques2->extract_indices_check(actual_base2.data(), r);
      assert(check_vtx2 == vtx2);

      // Check that the clique counts match up
      if (cliques->get_count(vtx) != cliques2->get_count(vtx2)) {
        std::cout << "Vtx1: " << static_cast<long>(vtx) << ", vtx2: " << static_cast<long>(vtx2) << std::endl;
        std::cout << "count1: " << cliques->get_count(vtx) << ", count2: " << cliques2->get_count(vtx2) << std::endl;
        std::cout << "Vert in 1: " << std::endl;
        for (size_t j = 0; j < r + 1; j++) {
          std::cout << actual_base[j] << ", " << std::endl;
        }
        std::cout << "Vert in 2: " << std::endl;
        for (size_t j = 0; j < r + 1; j++) {
          std::cout << actual_base2[j] << ", " << std::endl;
        }
        std::cout << "ppc1: " << per_processor_counts[vtx] << std::endl;
        std::cout << "ppc2: " << per_processor_counts2[vtx2] << std::endl;
        fflush(stdout);
      }
      assert(cliques->get_count(vtx) == cliques2->get_count(vtx2));

      // Check that cliques and cliques2 are thinking about the same clique
      parlay::sample_sort_inplace (make_slice(actual_base), [&](const uintE& u, const uintE&  v) {
          return u < v;
      });
      parlay::sample_sort_inplace (make_slice(actual_base2), [&](const uintE& u, const uintE&  v) {
          return u < v;
      });
      for (size_t j = 0; j < r + 1; j++) {
        assert(actual_base[j] == actual_base2[j]);
      }
    }

    parallel_for(0, num_entries, [&](size_t i) {
      per_processor_counts[i] = 0;
    });
    parallel_for( 0, num_entries2, [&](size_t i) {
      per_processor_counts2[i] = 0;
    });

    auto apply_f = [&](size_t i) -> std::optional<std::tuple<unsigned __int128, bucket_t>> {
      auto v = std::get<0>(D_filter[i]);
      bucket_t bucket = std::get<1>(D_filter[i]);
      if (v != num_entries + 1) {
        if (still_active[v] != 2 && still_active[v] != 1) return wrap(v, bucket);
      }
      return std::nullopt;
    };
    auto apply_f2 = [&](size_t i) -> std::optional<std::tuple<unsigned __int128, bucket_t>> {
      auto v = std::get<0>(D_filter2[i]);
      bucket_t bucket = std::get<1>(D_filter2[i]);
      if (v != num_entries2 + 1) {
        if (still_active2[v] != 2 && still_active2[v] != 1) return wrap(v, bucket);
      }
      return std::nullopt;
    };

    t_update.start();
    b.update_buckets(apply_f, filter_size);
    b2.update_buckets(apply_f2, filter_size2);
    t_update.stop();

    rounds++;
  }

  t_extract.next("### Peel Extract time: ");
  t_count.next("### Peel Count time: ");
  t_update.next("### Peel Update time: ");
  t_x.next("Inner counting: ");

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  std::cout.precision(17);
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "clique core: " << max_bkt << std::endl;

  //b.del();
  free(still_active);

  return D;
}

} // end namespace gbbs