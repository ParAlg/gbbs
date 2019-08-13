// Based on AsyncST/CC.C
#pragma once

#include <algorithm>
#include "ligra.h"
#include "utils/stats.h"
#include "../benchmark/MIS.h"
#include "pbbslib/random_shuffle.h"
#include <unordered_map>

template <class G, class S>
inline symmetric_graph<symmetricVertex, pbbs::empty> renumber(G& GA, S& perm) {
  using W = typename G::weight_type;
  uintE n = GA.n;

  // build the new graph
  auto offsets = sequence<size_t>(n);
  parallel_for(0, n, [&] (size_t i) {
    uintE our_id = perm[i];
    auto count_f = [&](uintE src, uintE ngh, const W& wgh) {
      uintE ngh_id = perm[ngh];
      return ngh_id < our_id;
    };
    offsets[our_id] = GA.get_vertex(i).countOutNgh(i, count_f);
  }, 1);

  size_t new_m = pbbslib::scan_add_inplace(offsets.slice());
  cout << "new_m = " << new_m << endl;

  auto new_edges = pbbs::sequence<uintE>(new_m);
  parallel_for(0, n, [&] (size_t i) {
    uintE new_id = perm[i];
    size_t offset = offsets[new_id];
    size_t directed_degree = ((new_id == n-1) ? new_m : offsets[new_id+1]) - offset;

    size_t k = 0;
    auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      uintE id_v = perm[v];
      if (id_v < new_id) {
        assert(offset + k <= new_m);
        new_edges[offset + k] = id_v;
        k++;
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f, false); // run seq
    assert(directed_degree == k);
    auto seq = pbbslib::make_sequence(new_edges.begin() + offset, k);
    // sort by id
    pbbs::sample_sort_inplace(seq, std::less<uintE>());
  });

  auto ss = pbbs::sequence<uintE>(n);
  auto v = pbbs::sequence<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    uintE new_id = perm[i];
    uintE o = offsets[new_id];
    uintE degree = ((new_id == n-1) ? new_m : offsets[new_id+1]) - o;
    v[new_id].degree = degree;
    v[new_id].offset = o;
  });

  auto vv = v.to_array();
  auto edge_arr = (std::tuple<uintE, pbbs::empty>*)new_edges.to_array();
  return symmetric_graph<symmetricVertex, pbbs::empty>(vv, n, new_m, get_deletion_fn(vv, edge_arr), edge_arr);
}

template <class G>
std::tuple<bool, size_t> YoshidaCache(G& GA, uintE u, std::unordered_map<uintE, bool>& cache) {
  using W = typename G::weight_type;
  bool in_mis = true;
  size_t total_work = 1;
  bool quit = false;
  auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) -> bool {
    bool ngh_in = false;
    size_t ngh_total_work = 0;
    auto lookup = cache.find(v);
    if (lookup == cache.end()) {
      std::tie(ngh_in, ngh_total_work) = YoshidaCache(GA, v, cache);
      total_work += ngh_total_work + 1; // recursive work done below, plus one for this neighbor
      cache.insert(std::make_pair(v, ngh_in));
    } else {
      ngh_in = lookup->second;
      total_work += 1; // 1 for this neighbor (pruned)
    }
    if (ngh_in) {
      in_mis = false;
      quit = true;
      return false; // break
    } else {
      return true; // keep going, ngh was not in MIS
    }
  };
  if (GA.get_vertex(u).getOutDegree() > 0) {
    GA.get_vertex(u).decodeOutNghSparseSeqCond(u, map_f);
  } // otherwise base case; return that we are in mis, and total_work is 1 (to check this vertex)
  return std::make_tuple(in_mis, total_work);
}

template <class G>
std::tuple<bool, size_t> Yoshida(G& GA, uintE u) {
  using W = typename G::weight_type;
  bool in_mis = true;
  size_t total_work = 1;
  auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) -> bool {
    if (v >= u) {
      return false;
    }
    bool ngh_in = false;
    size_t ngh_total_work = 0;
    std::tie(ngh_in, ngh_total_work) = Yoshida(GA, v);
    total_work += ngh_total_work + 1; // recursive work done below, plus one for this neighbor
    if (ngh_in) {
      in_mis = false;
      return false; // break
    } else {
      return true; // keep going, ngh was not in MIS
    }
  };
  if (GA.get_vertex(u).getOutDegree() > 0) {
    GA.get_vertex(u).decodeOutNghSparseSeqCond(u, map_f);
  } // otherwise base case; return that we are in mis, and total_work is 1 (to check this vertex)
  return std::make_tuple(in_mis, total_work);
}

template <class G>
auto YoshidaMIS(G& GA) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;

  auto in_mis = pbbs::sequence<bool>(n, false);

  auto perm = pbbslib::random_permutation<uintE>(n);

  auto GR = renumber(GA, perm);

  timer search_t; search_t.start();
  cout << "GR.n = " << GR.n << " GR.m = " << GR.m << endl;
  auto search_lengths = pbbs::sequence<size_t>(n, static_cast<size_t>(0));
  parallel_for(0, n, [&] (size_t i) {
    bool is_in_mis = false;
    size_t search_length = 0;
    uintE new_id = perm[i];
    // auto mm = std::unordered_map<uintE, bool>();
    std::tie(is_in_mis, search_length) = Yoshida(GR, new_id);
    in_mis[i] = is_in_mis;
    search_lengths[i] = search_length;
  }, 512);
  search_t.stop(); search_t.reportTotal("yoshida search time");

//  auto mm = std::unordered_map<uintE, bool>();
//  YoshidaCache(GR, 478623, mm);
//  cout << mm.size() << endl;

//  size_t mis_size = 0;
//  for (size_t i=0; i<n; i++) {
//    bool is_in_mis = false;
//    size_t total_work = 0;
//    uintE new_id = perm[i];
//    std::tie(is_in_mis, total_work) = Yoshida(GR, new_id);
//    cout << "id = " << i << " is_in_mis = " << is_in_mis << " total work = " << total_work << endl;
//    search_lengths[i] = total_work;
//    mis_size += static_cast<size_t>(is_in_mis);
//  }

  size_t total_work = pbbslib::reduce_add(search_lengths);
  cout << "total work = " << total_work << endl;

  size_t max_work = pbbslib::reduce_max(search_lengths);
  cout << "max work = " << max_work << endl;

//  MIS_rootset::verify_mis(GA, in_mis);
//  cout << "mis size = " << mis_size << endl;
}

