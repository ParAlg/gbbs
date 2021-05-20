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

    auto base_f = [&](sequence<uintE>& base){
      table->insert(base, r, k);
    };
    auto tots = sequence<size_t>(DG.n, size_t{0});

    auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k, DG.n, true, true); };
    auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } }; //induced->del();
    //parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, DG.n, [&](size_t i, HybridSpace_lw* induced) {
    for(size_t i =0; i < DG.n; i++) {
      HybridSpace_lw* induced = new HybridSpace_lw();
      init_induced(induced);
        if (DG.get_vertex(i).getOutDegree() != 0) {
          induced->setup(DG, k, i);
          auto base = sequence<uintE>(k + 1);
          base[0] = i;
          tots[i] = NKCliqueDir_fast_hybrid_rec(DG, 1, k, induced, base_f, base);
        } else tots[i] = 0;
      finish_induced(induced);
    }//, 1, false);
    double tt2 = t2.stop();

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
    size_t ss;
    list_buffer(size_t s){
      ss = s;
      num_workers2 = num_workers();
      buffer = 64;
      int buffer2 = 64;
      //list = sequence<size_t>(s + buffer2 * num_workers2, static_cast<size_t>(UINT_E_MAX));
      list = sequence<size_t>(s, static_cast<size_t>(UINT_E_MAX));
      std::cout << "list size: " << list.size() << std::endl;
      //starts = sequence<size_t>(num_workers2, [&](size_t i){return i * buffer2;});
      //next = num_workers2 * buffer2;
      //to_pack = sequence<bool>(s + buffer2 * num_workers2, true);
      next = 0;
    }
    void add(size_t index) {
      //std::cout << "Add: " << index << std::endl; fflush(stdout);
      size_t use_next = pbbs::fetch_and_add(&next, 1);
      list[use_next] = index;
      /*size_t worker = worker_id();
      list[starts[worker]] = index;
      starts[worker]++;
      if (starts[worker] % buffer == 0) {
        size_t use_next = pbbs::fetch_and_add(&next, buffer);
        starts[worker] = use_next;
      }*/
    }

    template <class I>
    size_t filter(I update_changed, sequence<double>& per_processor_counts) {
      //std::cout << "Next: "<< next << std::endl;
      parallel_for(0, next, [&](size_t worker) {
        assert(list[worker] != UINT_E_MAX);
        assert(per_processor_counts[list[worker]] != 0);
        update_changed(per_processor_counts, worker, list[worker]);
      });
      return next;
/*
      parallel_for(0, num_workers2, [&](size_t worker) {
        size_t divide = starts[worker] / buffer;
        for (size_t j = starts[worker]; j < (divide + 1) * buffer; j++) {
          to_pack[j] = false;
        }
      });
      // Pack out 0 to next of list into pack
      parallel_for(0, next, [&] (size_t i) {
        if (list[i] != UINT_E_MAX)//(to_pack[i])
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
      return next;*/
    }

    void reset() {
      /*parallel_for (0, num_workers2, [&] (size_t j) {
        starts[j] = j * buffer;
      });
      parallel_for (0, ss + buffer * num_workers2, [&] (size_t j) {
        list[j] = UINT_E_MAX;
      });
      next = num_workers2 * buffer;*/
      next = 0;
    }
};

template <class Graph>
bool is_edge(Graph& DG, uintE v, uintE u) {
  using W = typename Graph::weight_type;
  bool is = false;
  auto map_f = [&] (const uintE& src, const uintE& vv, const W& wgh) {
    if (vv == u) is = true;
    };
    DG.get_vertex(v).mapOutNgh(v, map_f, false);
    return is;
}

template <class Graph, class Graph2, class F, class I, class T>
inline size_t cliqueUpdate(Graph& G, Graph2& DG, size_t r, 
size_t k, size_t max_deg, bool label, F get_active, size_t active_size,
  size_t granularity, char* still_active, sequence<uintE> &rank, 
  sequence<double>& per_processor_counts, 
  bool do_update_changed, I update_changed,
  T* cliques, size_t n, list_buffer& count_idxs, timer& t1) {

  // Mark every vertex in the active set
  parallel_for (0, active_size, [&] (size_t j) {
    auto index = get_active(j); //cliques->find_index(get_active(j));
    still_active[index] = 1;
    }, 2048);

  auto is_active = [&](size_t index) {
    return still_active[index] == 1;
  };
  auto is_inactive = [&](size_t index) {
    return still_active[index] == 2;
  };
  auto update_d = [&](sequence<uintE>& base){
    cliques->extract_indices(base, is_active, is_inactive, [&](std::size_t index, double val){
      double ct = pbbs::fetch_and_add(&(per_processor_counts[index]), val);
      if (ct == 0 && val != 0) {
        count_idxs.add(index);
      }
    }, r, k);
  };

  // Set up space for clique counting
  auto init_induced = [&](HybridSpace_lw* induced) { induced->alloc(max_deg, k-r, G.n, true, true, true); };
  auto finish_induced = [&](HybridSpace_lw* induced) { if (induced != nullptr) { delete induced; } };

t1.start();
  // Clique count updates
  //std::cout << "Start setup nucleus" << std::endl; fflush(stdout);

  //parallel_for_alloc<HybridSpace_lw>(init_induced, finish_induced, 0, active_size,
  //                                   [&](size_t i, HybridSpace_lw* induced) {
  for(size_t i =0; i < active_size; i++) {
    HybridSpace_lw* induced = new HybridSpace_lw();
    init_induced(induced);

    auto x = get_active(i);
    auto base = sequence<uintE>(k + 1, [](size_t j){return UINT_E_MAX;});
    cliques->extract_clique(x, base, G, k);
    // Fill base[k] ... base[k-r+1] and base[0]
    induced->setup_nucleus(G, DG, k, base, r);
    NKCliqueDir_fast_hybrid_rec(DG, 1, k-r, induced, update_d, base);
    finish_induced(induced);
  }//, 1, true); //granularity
  //std::cout << "End setup nucleus" << std::endl; fflush(stdout);
t1.stop();

  // Perform update_changed on each vertex with changed clique counts
  std::size_t num_count_idxs = 0;
  num_count_idxs = count_idxs.filter(update_changed, per_processor_counts);
  count_idxs.reset();

  // Mark every vertex in the active set as deleted
  parallel_for (0, active_size, [&] (size_t j) {
    auto index = get_active(j); 
    still_active[index] = 2;}, 2048);

  return num_count_idxs; //count_idxs[0];
}

template <typename bucket_t, class Graph, class Graph2, class T>
sequence<bucket_t> Peel(Graph& G, Graph2& DG, size_t r, size_t k, 
  T* cliques, sequence<uintE> &rank,
  size_t num_buckets=16) {
    k--; r--;
  timer t2; t2.start();

  size_t num_entries = cliques->return_total();
  auto D = sequence<bucket_t>(num_entries, [&](size_t i) -> bucket_t { 
    return cliques->get_count(i);
  });

  auto D_filter = sequence<std::tuple<uintE, bucket_t>>(num_entries);

  auto b = make_vertex_custom_buckets<bucket_t>(num_entries, D, increasing, num_buckets);

  auto per_processor_counts = sequence<double>(num_entries , static_cast<double>(0));
  
  list_buffer count_idxs(num_entries);

  char* still_active = (char*) calloc(num_entries, sizeof(char));
  size_t max_deg = induced_hybrid::get_max_deg(G); // could instead do max_deg of active?

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
  size_t iter = 0;

  while (finished != num_entries) {
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
    iter++;

    size_t granularity = (cur_bkt * active_size < 10000) ? 1024 : 1;

    size_t filter_size = 0;

      auto update_changed = [&](sequence<double>& ppc, size_t i, uintE v){

        if (v == UINT_E_MAX) {
          D_filter[i] = std::make_tuple(num_entries + 1, 0);
          return;
        }
    
        assert(ppc[v] != 0);
        if (ppc[v] == 0) D_filter[i] = std::make_tuple(num_entries + 1, 0);
        else {
          bucket_t deg = D[v];
          //bucket_t deg = cliques->get_count(v);
          assert(deg > cur_bkt);
          auto val = cliques->update_count(v, (size_t) ppc[v]);
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
     filter_size = cliqueUpdate(G, DG, r, k, max_deg, true, get_active, active_size, 
     granularity, still_active, rank, per_processor_counts,
      true, update_changed, cliques, num_entries, count_idxs, t_x);
      t_count.stop();

    auto apply_f = [&](size_t i) -> std::optional<std::tuple<unsigned __int128, bucket_t>> {
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
  }

  t_extract.reportTotal("### Peel Extract time: ");
  t_count.reportTotal("### Peel Count time: ");
  t_update.reportTotal("### Peel Update time: ");
  t_x.reportTotal("Inner counting: ");

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  std::cout.precision(17);
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "clique core: " << max_bkt << std::endl;
  if (use_max_density) std::cout << "max density: " << max_density << std::endl;

  b.del();
  free(still_active);

  return D;
}

//*************************************************************VERIFICATION CODE**********


template <typename bucket_t, class Graph, class Graph2, class T, class T2>
sequence<bucket_t> Peel_verify(Graph& G, Graph2& DG, size_t r, size_t k, 
  T* cliques, T2* cliques2, sequence<uintE> &rank, size_t num_rounds_verify=2,
  size_t num_buckets=16) {
    k--; r--;
  timer t2; t2.start();

  size_t num_entries = cliques->return_total();
  size_t num_entries2 = cliques2->return_total();

  auto D = sequence<bucket_t>(num_entries, [&](size_t i) -> bucket_t { 
    return cliques->get_count(i);
  });
  auto D2 = sequence<bucket_t>(num_entries2, [&](size_t i) -> bucket_t { 
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
      true, update_changed, cliques, num_entries, count_idxs, t_x);
    
    size_t filter_size2 = cliqueUpdate(G, DG, r, k, max_deg, true, get_active2, active_size2, 
      granularity, still_active2, rank, per_processor_counts2,
      true, update_changed2, cliques2, num_entries2, count_idxs2, t_x);
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
      pbbslib::sample_sort_inplace (actual_base.slice(), [&](const uintE& u, const uintE&  v) {
          return u < v;
      });
      pbbslib::sample_sort_inplace (actual_base2.slice(), [&](const uintE& u, const uintE&  v) {
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

  t_extract.reportTotal("### Peel Extract time: ");
  t_count.reportTotal("### Peel Count time: ");
  t_update.reportTotal("### Peel Update time: ");
  t_x.reportTotal("Inner counting: ");

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  std::cout.precision(17);
  std::cout << "rho: " << rounds << std::endl;
  std::cout << "clique core: " << max_bkt << std::endl;

  b.del();
  free(still_active);

  return D;
}

} // end namespace gbbs