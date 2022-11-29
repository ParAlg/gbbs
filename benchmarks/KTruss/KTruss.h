// Copyright (c) 2019 Laxman Dhulipala, Guy Blelloch, and Julian Shun
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

#include <cassert>

#include "gbbs/bridge.h"
#include "gbbs/bucket.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"
#include "gbbs/helpers/dyn_arr.h"
#include "gbbs/helpers/sparse_table.h"

#include "benchmarks/Connectivity/SimpleUnionAsync/Connectivity.h"

#include "truss_utils.h"

namespace gbbs {


sequence<uintE> GetBoundaryIndices(
    std::size_t num_keys,
    const std::function<bool(std::size_t, std::size_t)>& key_eq_func) {
      uintE null_key = UINT_E_MAX;
  sequence<uintE> mark_keys(num_keys + 1);
  parlay::parallel_for(0, num_keys, [&](std::size_t i) {
    if (i != 0 && key_eq_func(i, i - 1))
      mark_keys[i] = null_key;
    else
      mark_keys[i] = i;
  });
  mark_keys[num_keys] = num_keys;
  auto vert_buckets = 
      parlay::filter(mark_keys, [&](uintE x) -> bool { return x != null_key; });
  return vert_buckets;
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

class EfficientConnectWhilePeeling {
  public:
    EfficientConnectWhilePeeling() {}
    EfficientConnectWhilePeeling(size_t _n) {
      n = _n;
      uf = gbbs::simple_union_find::SimpleUnionAsyncStruct(n);
      links = sequence<uintE>::from_function(n, [&](size_t s) { return UINT_E_MAX; });
    }
    
    void initialize(size_t _n);

    template<class X, class Y, class F>
    void link(X a, Y b, F& cores);

    template<class X, class Y, class F>
    void check_equal_for_merge(X a, Y b, F& cores);

    template<class bucket_t>
    void init(bucket_t cur_bkt);

    gbbs::simple_union_find::SimpleUnionAsyncStruct uf =  gbbs::simple_union_find::SimpleUnionAsyncStruct(0);
    sequence<uintE> links;
    size_t n; // table size
};

void EfficientConnectWhilePeeling::initialize(size_t _n)  {
  this->n = _n;
  this->uf = gbbs::simple_union_find::SimpleUnionAsyncStruct(this->n);
  this->links = sequence<uintE>::from_function(this->n, [&](size_t s) { return UINT_E_MAX; });
}

template<class X, class Y, class F>
void EfficientConnectWhilePeeling::check_equal_for_merge(X a, Y b, F& cores) {
  if (cores(a) == cores(b)) {
    this->uf.unite(a, b);
  } else {
    auto link_b = links[b];
    if (link_b != UINT_E_MAX && cores(link_b) >= cores(a)) this->check_equal_for_merge(a, link_b, cores);
  }
}

template<class X, class Y, class F>
void EfficientConnectWhilePeeling::link(X a, Y b, F& cores) {
  a = simple_union_find::find_compress(a, this->uf.parents);
  b = simple_union_find::find_compress(b, this->uf.parents);

  if (cores(a) == cores(b)) {
    this->uf.unite(a, b);
    uintE parent = simple_union_find::find_compress(a, this->uf.parents);
    auto link_a = links[a]; auto link_b = links[b];
    if (link_a != UINT_E_MAX && parent != a) this->link(link_a, parent, cores);
    if (link_b != UINT_E_MAX && parent != b) this->link(link_b, parent, cores);
  }
  else if (cores(a) < cores(b)) {
      gbbs::uintE c = links[b];
      while (true) {
        c = links[b];
        if (c == UINT_E_MAX) {
          if (gbbs::atomic_compare_and_swap<uintE>(&(links[b]), UINT_E_MAX, a)) break;
        } else if (cores(c) < cores(a)) { // || (cores(c) == cores(a) && a < c)
          if (gbbs::atomic_compare_and_swap<uintE>(&(links[b]), c, a)) {
            auto parent_b = simple_union_find::find_compress(b, this->uf.parents);
            if (b != parent_b) this->link(a, parent_b, cores);
            this->link(a, c, cores);
            break;
          }
        } else {
          this->link(a, c, cores);
          break;
        }
      }
  }
  else {
    this->link(b, a, cores);
  }
}

template <class bucket_t>
void EfficientConnectWhilePeeling::init(bucket_t cur_bkt) {
}


class ConnectWhilePeeling {
  public:
    ConnectWhilePeeling() {}
    ConnectWhilePeeling(size_t _n){
      n = _n;
    }

    void initialize(size_t _n);

    template<class X, class Y, class F>
    void link(X x, Y index, F& cores);

    template<class bucket_t>
    void init(bucket_t cur_bkt);
    size_t n; // table size
    std::vector<gbbs::simple_union_find::SimpleUnionAsyncStruct> set_uf;
    std::vector<uintE> set_core;
};

void ConnectWhilePeeling::initialize(size_t _n) { this->n = _n; }

template<class X, class Y, class F>
void ConnectWhilePeeling::link(X x, Y index, F& cores) {
  parallel_for(0, set_uf.size(), [&](size_t idx){
    if (cores(index) >= set_core[idx]) set_uf[idx].unite(x, index);
  });
}

template <class bucket_t>
void ConnectWhilePeeling::init(bucket_t cur_bkt) {
  set_uf.push_back(gbbs::simple_union_find::SimpleUnionAsyncStruct(n));
  set_core.push_back(cur_bkt);
}

template <class Table>
inline std::vector<uintE> construct_nd_connectivity_from_connect(EfficientConnectWhilePeeling& cwp, Table& trussness_multi){
  auto parents = cwp.uf.finish();

  auto n = trussness_multi.size();
  auto is_valid = [&](uintE id) { return std::get<1>(trussness_multi.big_table[id]) != UINT_E_MAX; };

  // Sort vertices from highest core # to lowest
  auto sort_by_parent = [&](uintE p, uintE q) {
      return parents[p] < parents[q];
  };
  auto sorted_vert = sequence<uintE>::from_function(n, [&](size_t i) { return i; });
  parlay::sample_sort_inplace(make_slice(sorted_vert), sort_by_parent);

  auto parent_eq_func = [&](size_t i, size_t j) {return parents[sorted_vert[i]] == parents[sorted_vert[j]];};
  auto vert_buckets = GetBoundaryIndices(n, parent_eq_func);


  std::vector<uintE> connectivity_tree(n);
  parallel_for(0, n, [&](std::size_t i){connectivity_tree[i] = UINT_E_MAX;});
  uintE prev_max_parent = n;
  parallel_for (0, vert_buckets.size()-1, [&](size_t i) {
    size_t start_index = vert_buckets[i];
    size_t end_index = vert_buckets[i + 1];

    auto first_x = sorted_vert[start_index];
    if (is_valid(first_x)) {
      parallel_for(0, end_index - start_index, [&](size_t a){
        if (is_valid(sorted_vert[start_index + a])) {
          connectivity_tree[sorted_vert[start_index + a]] = prev_max_parent + i;
        }
      });
    }
  });
  prev_max_parent += vert_buckets.size() - 1;
  //std::cout << "Finish first pass" << std::endl; fflush(stdout);
  connectivity_tree.resize(prev_max_parent);
  parallel_for(n, prev_max_parent, [&](std::size_t i){connectivity_tree[i] = UINT_E_MAX;});

  for (size_t i = 0; i < cwp.links.size(); i++) {
    if (!is_valid(i)) continue;
    if (cwp.links[i] == UINT_E_MAX) continue;
    if (i == parents[i]) {
      connectivity_tree[connectivity_tree[i]] = connectivity_tree[cwp.links[i]];
    }
  }
  return connectivity_tree;
}

template <class Table>
inline std::vector<uintE> construct_nd_connectivity_from_connect(ConnectWhilePeeling& connect_with_peeling, Table& trussness_multi){
  auto n = trussness_multi.size();
  auto is_valid = [&](uintE id) { return std::get<1>(trussness_multi.big_table[id]) != UINT_E_MAX; };

  std::vector<uintE> connectivity_tree(n);
  sequence<uintE> prev_parent = sequence<uintE>::from_function(n, [&](size_t i){ return i; });
  uintE prev_max_parent = n;
  parallel_for(0, n, [&](std::size_t i){connectivity_tree[i] = UINT_E_MAX;});
  for (long idx = connect_with_peeling.set_uf.size() - 1; idx >= 0; idx--) {
    connectivity_tree.resize(prev_max_parent, UINT_E_MAX);
    parallel_for(0, n, [&](size_t l) { gbbs::simple_union_find::find_compress(l, connect_with_peeling.set_uf[idx].parents); });

    parallel_for(0, n, [&](size_t l){
      if (is_valid(l)) {
        connectivity_tree[prev_parent[l]] = prev_max_parent + connect_with_peeling.set_uf[idx].parents[l];
        // Update previous parent
        prev_parent[l] = connectivity_tree[prev_parent[l]];
      }
    });
    // Update previous max parent
    prev_max_parent += n;
  }
  return connectivity_tree;
}

// Sort vertices from highest core # to lowest core #
// Take each "bucket" of vertices in each core #
// Maintain connectivity UF throughout all rounds
// Run connectivity considering only vertices previously considered, or vertices
// in the current bucket
// k and r here should be the same as that used in cliqueUpdate
template <class Graph, class Table>
inline std::vector<uintE> construct_nd_connectivity(Graph& GA, Table& trussness_multi){
  using edge_t = uintE;
  using bucket_t = uintE;
  using trussness_t = uintE;

  auto n = trussness_multi.size();

  auto get_trussness_and_id = [&trussness_multi](uintE u, uintE v) {
    // Precondition: uv is an edge in G.
    edge_t id = trussness_multi.idx(u, v);
    trussness_t truss = std::get<1>(trussness_multi.big_table[id]);
    return std::make_tuple(truss, id);
  };

  auto cores = [&](uintE id) -> uintE {
    auto truss = std::get<1>(trussness_multi.big_table[id]);
    if (truss == std::numeric_limits<int>::max()) return 0;
    if (truss != UINT_E_MAX) return truss - 1;
    return truss;
  };

  // Sort vertices from highest core # to lowest core #
  auto get_core = [&](uintE p, uintE q){
    return cores(p) > cores(q);
  };
  auto sorted_vert = sequence<uintE>::from_function(n, [&](size_t i) { return i; });
  parlay::sample_sort_inplace(make_slice(sorted_vert), get_core);


  sequence<uintE> mark_keys(n + 1);
  parallel_for(0, n, [&](std::size_t i) {
    if (i != 0 && cores(sorted_vert[i]) == cores(sorted_vert[i-1]))
      mark_keys[i] = UINT_E_MAX;
    else
      mark_keys[i] = i;
  });
  mark_keys[n] = n;
  auto vert_buckets = 
      parlay::filter(mark_keys, [&](uintE x) -> bool { return x != UINT_E_MAX; });

  auto uf = gbbs::simple_union_find::SimpleUnionAsyncStruct(n);
  // TODO(jeshi): This isn't parallel, but I'd just like easy resizing atm
  std::vector<uintE> connectivity_tree(n);
  sequence<uintE> prev_parent = sequence<uintE>::from_function(n, [&](size_t i){ return i; });
  uintE prev_max_parent = n;
  //uintE prev_offset = 0;
  parallel_for(0, n, [&](std::size_t i){connectivity_tree[i] = UINT_E_MAX;});

  for (size_t i = 0; i < vert_buckets.size()-1; i++) {
    size_t start_index = vert_buckets[i];
    size_t end_index = vert_buckets[i + 1];

    auto first_x = sorted_vert[start_index];
    auto first_current_core = cores(first_x);
    // TODO: to get rid of this, gotta check for valid indices
    if (first_current_core != UINT_E_MAX && first_current_core != 0) {

      auto is_active = [&](size_t index) {
        return cores(index) == first_current_core;
      };
      auto is_inactive = [&](size_t index) {
        return cores(index) < first_current_core;
      };

    parallel_for(start_index, end_index, [&](size_t j){
      // Vertices are those given by sorted_vert[j]
      auto x = sorted_vert[j];
  
      auto current_core = cores(x);
      assert(current_core == first_current_core);

      // Steps:
      // 1. Extract vertices given by index -- these will be the r-clique X
      // 2. Find all s-clique containing the r-clique X
      // 3. For each s-clique containing X, iterate over all combinations of
      // r vertices in that s-clique; these form r-clique X'
      // 4. Union X and X'
      auto unite_func = [&](edge_t a, edge_t b){ uf.unite(a, b); };
      uintE u = trussness_multi.u_for_id(x);
      uintE v = std::get<0>(trussness_multi.big_table[x]);

      truss_utils::do_union_things<edge_t, trussness_t>(GA, x, u, v, get_trussness_and_id, unite_func, is_inactive);


    }, 1, true); // granularity
    }

    connectivity_tree.resize(prev_max_parent, UINT_E_MAX);

    parallel_for(0, n, [&](size_t l) { gbbs::simple_union_find::find_compress(l, uf.parents); });

    sequence<uintE> map_parents = sequence<uintE>::from_function(n, [&](std::size_t l){return 0;});
    parallel_for(0, n, [&](size_t l){
      if (cores(l) != UINT_E_MAX && cores(l) >= first_current_core) map_parents[uf.parents[l]] = 1;
    });
    auto max_parent = parlay::scan_inplace(make_slice(map_parents));

    parallel_for(0, n, [&](size_t l){
      if (cores(l) != UINT_E_MAX && cores(l) >= first_current_core) {
        assert(prev_parent[l] < prev_max_parent);
        connectivity_tree[prev_parent[l]] = prev_max_parent + map_parents[uf.parents[l]];  //***: uf.parents[l];
        // Update previous parent
        prev_parent[l] = connectivity_tree[prev_parent[l]];
      }
    });
    prev_max_parent += max_parent;
  }

  return connectivity_tree;
}

template <class Graph, class T>
void CountCliquesNucPND(Graph& DG, T apply_func) {
      //auto tots = sequence<size_t>(DG.n, size_t{0});
      parallel_for(0, DG.n, [&](size_t i) {
        auto vert_i = DG.get_vertex(i);
        auto iter_i = vert_i.out_neighbors().get_iter();
        for (std::size_t j = 0; j < vert_i.out_degree(); j++) {
          auto x = std::get<0>(iter_i.cur());
          if (iter_i.has_next()) iter_i.next();
          std::vector<uintE> inter;
          truss_utils::truss_intersectionPND(DG, i, x, inter);
          //tots[i] += inter.size();
          for (std::size_t l = 0; l < inter.size(); l++) {
            apply_func(i, x, inter[l]);
          }
        }
      });
      //return pbbslib::reduce_add(tots);
  }

template <class Graph, class MT>
void initialize_trussness_values(Graph& GA, MT& multi_table, bool use_pnd = false) {
  using W = typename Graph::weight_type;

  timer it;
  it.start();
  GA.mapEdges([&](const uintE& u, const uintE& v, const W& wgh) {
    if (u < v) {
      multi_table.insert(u, std::make_tuple(v, 0));
    }
  });
  it.stop();
  it.next("insertion time");

  // 2. Triangle count, update trussness scores for each edge
  // 2.(a) Rank vertices based on degree
  auto rank = truss_utils::rankNodes(GA);

  // 2.(b) Direct edges to point from lower to higher rank vertices.
  auto pack_predicate = [&](const uintE& u, const uintE& v, const W& wgh) {
    return rank[u] < rank[v];
  };
  auto DG = filterGraph(GA, pack_predicate);
  std::cout << "Filtered graph to construct dirgraph: n = " << DG.n
            << " m = " << DG.m << std::endl;

  // Each triangle only found once---increment all three edges
  auto inc_truss_f = [&](const uintE& u, const uintE& v, const uintE& w) {
    multi_table.increment(u, v);
    multi_table.increment(u, w);
    multi_table.increment(v, w);
  };
  timer tct;
  tct.start();
  if (!use_pnd) truss_utils::TCDirected(DG, inc_truss_f);
  else CountCliquesNucPND(DG, inc_truss_f);
  tct.stop();
  tct.next("TC time");
}

// High-level desc:
// 1. Compute a hash table mapping each edge (u, v), u < v, to its trussness The
// initial trussness values are just the number of triangles each edge
// participates in.
//
// 2. Next, we compute a bucketing where each edge is represented by its index
// in the HT (locs without an edge are in the "infinity" bucket, and never
// leave.
//
// 3. Peel. Each peeling step removes the edges in bucket k, which implicitly
// fixes their trussness numbers.
//   3.a Each edge intersects its two endpoints using the edges in the original,
//   undirected graph. For each neighbor w in intersect(N(u), N(v)), we find the
//   (u, w) and (v, w) edges in the hash-table, check if we should decrement
//   them, and if so, insert them (their ids) into a hashtable.
//
//   3.b Get the entries of the HT, actually decrement their coreness, see if
//   their bucket needs to be updated and if so, update.
template <class Graph, class CWP>
truss_utils::multi_table<uintE, uintE, std::function<size_t(size_t)>> KTruss_ht(Graph& GA, CWP& connect_while_peeling, 
  size_t num_buckets = 16, bool inline_hierarchy = false, bool use_compact = true, bool use_pnd = false) {
  using W = typename Graph::weight_type;
  size_t n_edges = GA.m / 2;

  using edge_t = uintE;
  using bucket_t = uintE;
  using trussness_t = uintE;

  auto deg_lt = parlay::delayed_seq<uintE>(GA.n, [&](size_t i) {
    return GA.get_vertex(i).out_degree() < (1 << 15);
  });
  std::cout << "count = " << parlay::reduce(deg_lt) << std::endl;
  auto deg_lt_ct = parlay::delayed_seq<size_t>(GA.n, [&](size_t i) {
    if (GA.get_vertex(i).out_degree() < (1 << 15)) {
      return GA.get_vertex(i).out_degree();
    }
    return (uintE)0;
  });
  std::cout << "total degree = " << parlay::reduce(deg_lt_ct) << std::endl;

  auto counts = sequence<size_t>(GA.n, (size_t)0);
  parallel_for(0, GA.n, [&](size_t i) {
    auto d_i = GA.get_vertex(i).out_degree();
    bool d_i_lt = d_i <= (1 << 15);
    auto count_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      return u < v ? d_i_lt : 0;
    };
    counts[i] = GA.get_vertex(i).out_neighbors().count(count_f);
  });
  std::cout << "total lt ct = " << parlay::reduce(counts) << std::endl;

  std::tuple<edge_t, bucket_t> histogram_empty =
      std::make_tuple(std::numeric_limits<edge_t>::max(), 0);
  auto em = hist_table<edge_t, bucket_t>(histogram_empty, GA.m / 50);

  // Store the initial trussness of each edge in the trussness table.
  std::function<size_t(size_t)> get_size = [&](size_t vtx) -> size_t {
    auto count_f = [&](uintE u, uintE v, W& wgh) { return vtx < v; };
    return GA.get_vertex(vtx).out_neighbors().count(count_f);
  };
  auto trussness_multi =
      truss_utils::make_multi_table<uintE, uintE>(GA.n, UINT_E_MAX, get_size);

  // Note that this multi-table business is a performance optimization. The
  // previous version is somewhere in git history; we should measure how much
  // using a multi-table helps.
  //
  // Laxman (6/12): experiment with making the multi_table oriented by degree.
  // This requires an extra random access when handling an edge to place it in
  // the proper orientation. The simple ordering is to use ids, but using
  // low-deg --> high-deg has the advantage of reducing the max hash-table size,
  // which could improve locality.
  // * for small enough vertices, use an array instead of a hash table.

  // Initially stores #triangles incident/edge.
  initialize_trussness_values(GA, trussness_multi, use_pnd);

  // Initialize the bucket structure. #ids = trussness table size
  std::cout << "multi_size = " << trussness_multi.size() << std::endl;
  auto multi_size = trussness_multi.size();


  timer t2; t2.start();

  connect_while_peeling.initialize(multi_size);

  auto get_bkt = parlay::delayed_seq<uintE>(multi_size, [&](size_t i) {
    auto table_value =
        std::get<1>(trussness_multi.big_table[i]);  // the trussness.
    return (uintE)table_value;
  });
  auto b = make_buckets<edge_t, bucket_t>(trussness_multi.size(), get_bkt,
                                          increasing, num_buckets);

  // Stores edges idents that lose a triangle, including duplicates (MultiSet)
  auto hash_edge_id = [&](const edge_t& e) { return parlay::hash32(e); };
  auto decr_source_table = gbbs::make_sparse_table<edge_t, uintE>(
      1 << 20, std::make_tuple(std::numeric_limits<edge_t>::max(), (uintE)0),
      hash_edge_id);

  auto del_edges = gbbs::dyn_arr<edge_t>(6 * GA.n);
  auto actual_degree = sequence<uintE>::from_function(
      GA.n, [&](size_t i) { return GA.get_vertex(i).out_degree(); });

  auto get_trussness_and_id = [&trussness_multi](uintE u, uintE v) {
    // Precondition: uv is an edge in G.
    edge_t id = trussness_multi.idx(u, v);
    trussness_t truss = std::get<1>(trussness_multi.big_table[id]);
    return std::make_tuple(truss, id);
  };

  size_t data_structure_size = sizeof(decltype(trussness_multi)) + sizeof(std::tuple<uintE, uintE>) * trussness_multi.big_size;
  std::cout << "Data Structure Size: " << data_structure_size << std::endl; fflush(stdout);

  timer em_t, decrement_t, bt, ct, peeling_t;
  peeling_t.start();
  size_t finished = 0, k_max = 0;
  size_t iter = 0;
  uintE prev_bkt = 0;

  char* still_active = (char*) calloc(multi_size, sizeof(char));

  while (finished != n_edges) {
    bt.start();
    auto bkt = b.next_bucket();
    bt.stop();
    auto rem_edges = bkt.identifiers;
    if (rem_edges.size() == 0) {
      continue;
    }

    uintE k = bkt.id;
    finished += rem_edges.size();
    k_max = std::max(k_max, bkt.id);

    if (inline_hierarchy && prev_bkt != k && k != 0) {
      connect_while_peeling.init(k);
    }

    //std::cout << "k = " << k << " iter = " << iter << " #edges = " << rem_edges.size() << std::endl;

    if (k == 0) { // || finished == n_edges
      // No triangles incident to these edges. We set their trussness to MAX,
      // which is safe since there are no readers until we output.
      parallel_for(0, rem_edges.size(), [&](size_t i) {
        edge_t id = rem_edges[i];
        still_active[rem_edges[i]] = 2;
        std::get<1>(trussness_multi.big_table[id]) =
            std::numeric_limits<int>::max();  // UINT_E_MAX is reserved
      });
      continue;
    }

    parallel_for (0, rem_edges.size(), [&] (size_t j) {
        still_active[rem_edges[j]] = 1;
     });

    size_t e_size = 2 * k * rem_edges.size();
    size_t e_space_required = (size_t)1
                              << parlay::log2_up((size_t)(e_size * 1.2));

    // Resize the table that stores edge updates if necessary.
    decr_source_table.resize_no_copy(e_space_required);
    auto decr_tab = gbbs::make_sparse_table<edge_t, uintE>(
        decr_source_table.backing.begin(), e_space_required,
        std::make_tuple(std::numeric_limits<edge_t>::max(), (uintE)0),
        hash_edge_id, false /* do not clear */);

    auto cores_func = [&](size_t a) -> uintE {
      if (still_active[a] == 0) return multi_size + 1;
      if (still_active[a] == 1) return k;
      auto truss = std::get<1>(trussness_multi.big_table[a]);
      if (truss == std::numeric_limits<int>::max()) return 0;
      if (truss != UINT_E_MAX) return truss - 1;
      return truss; 
    };
    
    //    std::cout << "starting decrements" << std::endl;
    decrement_t.start();
    parallel_for(0, rem_edges.size(), 1, [&](size_t i) {
      edge_t id = rem_edges[i];
      uintE u = trussness_multi.u_for_id(id);
      uintE v = std::get<0>(trussness_multi.big_table[id]);

      auto to_link = [&](edge_t index) {
        if (index != id && still_active[index] != 0) {
          connect_while_peeling.link(id, index, cores_func);
        }
      };

      truss_utils::decrement_trussness<edge_t, trussness_t>(
          GA, id, u, v, decr_tab, get_trussness_and_id, k, use_pnd, inline_hierarchy, to_link);
    });
    decrement_t.stop();
    //    std::cout << "finished decrements" << std::endl;

    auto decr_edges = decr_tab.entries();
    parallel_for(0, decr_edges.size(), [&](size_t i) {
      auto id_and_ct = decr_edges[i];
      edge_t id = std::get<0>(id_and_ct);
      uintE triangles_removed = std::get<1>(id_and_ct);
      uintE current_deg = std::get<1>(trussness_multi.big_table[id]);
      assert(current_deg > k);
      uintE new_deg = std::max(current_deg - triangles_removed, k);
      std::get<1>(trussness_multi.big_table[id]) = new_deg;  // update
      std::get<1>(decr_edges[i]) = b.get_bucket(current_deg, new_deg);
    });

    auto rebucket_edges =
        parlay::filter(decr_edges, [&](const std::tuple<edge_t, uintE>& eb) {
          return std::get<1>(eb) != UINT_E_MAX;
        });
    auto edges_moved_f = [&](size_t i) {
      return std::optional<std::tuple<edge_t, bucket_t>>(rebucket_edges[i]);
    };

    bt.start();
    b.update_buckets(edges_moved_f, rebucket_edges.size());
    bt.stop();

     parallel_for (0, rem_edges.size(), [&] (size_t j) {
        still_active[rem_edges[j]] = 2;
     });

    //    auto apply_f = [&](const std::tuple<edge_t, uintE>& p)
    //        -> const Maybe<std::tuple<edge_t, bucket_t> > {
    //      edge_t id = std::get<0>(p);
    //      uintE triangles_removed = std::get<1>(p);
    //      uintE current_deg = std::get<1>(trussness_multi.big_table[id]);
    //      if (current_deg > k) {
    //        uintE new_deg = std::max(current_deg - triangles_removed, k);
    //        std::get<1>(trussness_multi.big_table[id]) = new_deg;
    //        bucket_t bkt = b.get_bucket(current_deg, new_deg);
    //        return wrap(id, bkt);
    //      }
    //      return Maybe<std::tuple<edge_t, bucket_t>>();
    //    };
    //
    //    em_t.start();
    //    sequence<edge_t> edge_seq = decrement_tab.entries();
    //    auto res = em.template edgeMapCount(edge_seq, apply_f);
    //    em_t.stop();
    //
    //    auto rebucket_f = [&] (size_t i) -> Maybe<std::tuple<edge_t,
    //    bucket_t>> {
    //      auto ret = Maybe<std::tuple<edge_t, bucket_t>>();
    //      ret.t = res.second[i];
    //      ret.exists = true;
    //      return ret;
    //    };
    //    auto edges_moved_map = parlay::delayed_seq<Maybe<std::tuple<edge_t,
    //    bucket_t>>>(res.first, rebucket_f);
    //    auto edges_moved_f = [&] (size_t i) { return edges_moved_map[i]; };
    //    bt.start();
    //    b.update_buckets(edges_moved_f, edges_moved_map.size());
    //    bt.stop();

    // Unmark edges removed in this round, and decrement their trussness.
    parallel_for(0, rem_edges.size(), [&](size_t i) {
      edge_t id = rem_edges[i];
      std::get<1>(trussness_multi.big_table[id]) -= 1;
    });

    // Clear the table storing the edge decrements.
    decr_tab.clear_table();
    iter++;

    if (use_compact) { 
    del_edges.copyIn(rem_edges, rem_edges.size());

    if (del_edges.size > 2 * GA.n) {
      ct.start();
      // compact
      std::cout << "compacting, " << del_edges.size << std::endl;
      // map over both endpoints, update counts using histogram
      // this is really a uintE seq, but edge_t >= uintE, and this way we can
      // re-use the same histogram structure.
      auto decr_seq = sequence<edge_t>(2 * del_edges.size);
      parallel_for(0, del_edges.size, [&](size_t i) {
        size_t fst = 2 * i;
        size_t snd = fst + 1;
        edge_t id = del_edges.A[i];
        uintE u = trussness_multi.u_for_id(id);
        uintE v = std::get<0>(trussness_multi.big_table[id]);
        decr_seq[fst] = u;
        decr_seq[snd] = v;
      });
      std::cout << "compacting 1, " << del_edges.size << std::endl;

      // returns only those vertices that have enough degree lost to warrant
      // packing them out. Again note that edge_t >= uintE
      auto apply_vtx_f = [&](const std::tuple<edge_t, uintE>& p)
          -> const std::optional<std::tuple<edge_t, uintE>> {
            uintE id = std::get<0>(p);
            uintE degree_lost = std::get<1>(p);
            actual_degree[id] -= degree_lost;
            // compare with GA.V[id]. this is the current space used for this
            // vtx.
            return std::nullopt;
          };
      std::cout << "compacting 2, " << del_edges.size << std::endl;

      em_t.start();
      auto vs = vertexSubset(GA.n, std::move(decr_seq));
      std::cout << "compacting 3, " << del_edges.size << std::endl;
      auto cond_f = [&](const uintE& u) { return true; };
      nghCount(GA, vs, cond_f, apply_vtx_f, em, no_output);
      em_t.stop();
      std::cout << "compacting 4, " << del_edges.size << std::endl;

      auto all_vertices =
          parlay::delayed_seq<uintE>(GA.n, [&](size_t i) { return i; });
      auto to_pack_seq = parlay::filter(all_vertices, [&](uintE u) {
        return 4 * actual_degree[u] >= GA.get_vertex(u).out_degree();
      });
      auto to_pack = vertexSubset(GA.n, std::move(to_pack_seq));
      std::cout << "compacting 5, " << del_edges.size << std::endl;

      auto pack_predicate = [&](const uintE& u, const uintE& ngh,
                                const W& wgh) {
        // return true iff edge is still alive
        trussness_t t_u_ngh;
        edge_t edgeid;
        std::tie(t_u_ngh, edgeid) = get_trussness_and_id(u, ngh);
        return t_u_ngh >= k;
      };
      edgeMapFilter(GA, to_pack, pack_predicate, pack_edges | no_output);

      del_edges.size = 0;  // reset dyn_arr
      ct.stop();
      std::cout << "Finished compacting." << std::endl;
    }
    }
    
    prev_bkt = k;
  }

  double tt2 = t2.stop();
  std::cout << "### Peel Running Time: " << tt2 << std::endl;

  peeling_t.stop();
  peeling_t.next("peeling time");
  ct.next("Compaction time");
  bt.next("Bucketing time");
  em_t.next("EdgeMap time");
  decrement_t.next("Decrement trussness time");

  // == Important: The actual trussness is the stored trussness value + 1.
  // Edges with trussness 0 had their values stored as
  // std::numeric_limits<int>::max()
  std::cout << "iters = " << iter << std::endl;

  return trussness_multi;
}

template <class Graph>
void KTruss_connect(Graph& GA, size_t num_buckets, bool inline_hierarchy, bool efficient_inline_hierarchy) {
    if (efficient_inline_hierarchy) inline_hierarchy = true;

  truss_utils::multi_table<uintE, uintE, std::function<size_t(size_t)>> multitable;
  
  EfficientConnectWhilePeeling ecwp;
  ConnectWhilePeeling connect_with_peeling;

  if (!efficient_inline_hierarchy) multitable = KTruss_ht(GA, connect_with_peeling, num_buckets, inline_hierarchy, false, false);
  else multitable = KTruss_ht(GA, ecwp, num_buckets, inline_hierarchy, false, false);

  std::vector<uintE> connect;
  if (!inline_hierarchy) {
    std::cout << "Running Connectivity" << std::endl;
    timer t3; t3.start();
    connect = construct_nd_connectivity(GA, multitable);
    double tt3 = t3.stop();
    std::cout << "### Connectivity Running Time: " << tt3 << std::endl;
  } else {
    std::cout << "Constructing tree" << std::endl;
    timer t3; t3.start();
    if (!efficient_inline_hierarchy) connect = construct_nd_connectivity_from_connect(connect_with_peeling, multitable);
    else connect = construct_nd_connectivity_from_connect(ecwp, multitable);
    double tt3 = t3.stop();
    std::cout << "### Connectivity Tree Running Time: " << tt3 << std::endl;
  }
}

}  // namespace gbbs
