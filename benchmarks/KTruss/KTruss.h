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

namespace ktruss {


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
    if (truss != UINT_E_MAX) return truss + 1;
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

unsigned approx_nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

template <class Graph, class CWP>
truss_utils::multi_table<uintE, uintE, std::function<size_t(size_t)>> KTruss_approx(Graph& GA, CWP& connect_while_peeling, 
  bool inline_hierarchy, std::string& approx_out_str, size_t num_buckets=16, double delta = 0.1, bool use_pow = false) {
  using W = typename Graph::weight_type;
  size_t n_edges = GA.m / 2;

  using edge_t = uintE;
  using bucket_t = uintE;
  using trussness_t = uintE;

      std::cout << "Delta which is really eps: " << delta << std::endl;
    // We're gonna ignore eps
    // delta is really eps

  uintE schooser = 3;

  double one_plus_delta = log(schooser + delta);
  auto get_bucket = [&](size_t deg) -> uintE {
    return ceil(log(1 + deg) / one_plus_delta);
  };

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


  // Initially stores #triangles incident/edge.
  initialize_trussness_values(GA, trussness_multi, false);

  // Initialize the bucket structure. #ids = trussness table size
  std::cout << "multi_size = " << trussness_multi.size() << std::endl;
  auto multi_size = trussness_multi.size();

  sequence<uintE> trussness_multi_capped(multi_size, UINT_E_MAX);


  timer t2; t2.start();

  auto get_trussness_and_id = [&trussness_multi](uintE u, uintE v) {
    // Precondition: uv is an edge in G.
    edge_t id = trussness_multi.idx(u, v);
    trussness_t truss = std::get<1>(trussness_multi.big_table[id]);
    return std::make_tuple(truss, id);
  };

  connect_while_peeling.initialize(multi_size);


  GA.mapEdges([&](const uintE& u, const uintE& v, const W& wgh) {
    if (u < v) {
      auto truss_tuple = get_trussness_and_id(u, v);
      auto bucket = std::get<0>(truss_tuple);
      if (bucket != 0 && bucket != UINT_E_MAX && bucket != std::numeric_limits<int>::max()) bucket = get_bucket(bucket);
      trussness_multi_capped[std::get<1>(truss_tuple)] = bucket;
    }
  });

  auto get_bkt = parlay::delayed_seq<uintE>(multi_size, [&](size_t i) -> uintE {
    return trussness_multi_capped[i];  // the trussness.
  });
  auto b = make_buckets<edge_t, bucket_t>(trussness_multi_capped.size(), get_bkt,
                                          increasing, num_buckets);

  // Stores edges idents that lose a triangle, including duplicates (MultiSet)
  auto hash_edge_id = [&](const edge_t& e) { return parlay::hash32(e); };
  auto decr_source_table = gbbs::make_sparse_table<edge_t, uintE>(
      1 << 20, std::make_tuple(std::numeric_limits<edge_t>::max(), (uintE)0),
      hash_edge_id);

  //auto del_edges = gbbs::dyn_arr<edge_t>(6 * GA.n);
  auto actual_degree = sequence<uintE>::from_function(
      GA.n, [&](size_t i) { return GA.get_vertex(i).out_degree(); });

  size_t data_structure_size = sizeof(decltype(trussness_multi)) + sizeof(std::tuple<uintE, uintE>) * trussness_multi.big_size;
  std::cout << "Data Structure Size: " << data_structure_size << std::endl; fflush(stdout);

  timer em_t, decrement_t, bt, ct, peeling_t;
  peeling_t.start();
  size_t finished = 0, k_max = 0;
  size_t iter = 0;
  uintE prev_bkt = 0;

  char* still_active = (char*) calloc(multi_size, sizeof(char));

  size_t rounds = 0;
  bucket_t cur_bkt = 0;
  bucket_t max_bkt = 0;
  double max_density = 0;
  bool use_max_density = false;
  auto num_entries = n_edges;

  size_t cur_inner_rounds = 0;
  size_t max_inner_rounds = log(num_entries) / log(1.0 + delta / schooser);


  while (finished < num_entries) {
    bt.start();
    auto bkt = b.next_bucket();
    bt.stop();
    auto rem_edges = bkt.identifiers;
    if (rem_edges.size() == 0) {
      continue;
    }

    uintE k = bkt.id;
    finished += rem_edges.size();
    
    if (prev_bkt != k) {
      cur_inner_rounds = 0;
    }
     // Check if we hit the threshold for inner peeling rounds.
    if (cur_inner_rounds == max_inner_rounds) {
      // new re-insertions will go to at least bucket k (one greater than before).
      k++;
      cur_inner_rounds = 0;
    }
     uintE lower_bound = ceil(pow((schooser + delta), k-1));

    if (inline_hierarchy && prev_bkt != k && k != 0) {
      connect_while_peeling.init(k);
    }
    k_max = std::max((uintE)k_max, (uintE)k);

    //std::cout << "k = " << k << " iter = " << iter << " #edges = " << rem_edges.size() << std::endl;

    if (k == 0) { // || finished == n_edges
      // No triangles incident to these edges. We set their trussness to MAX,
      // which is safe since there are no readers until we output.
      parallel_for(0, rem_edges.size(), [&](size_t i) {
        edge_t id = rem_edges[i];
        still_active[rem_edges[i]] = 2;
        std::get<1>(trussness_multi.big_table[id]) =
            std::numeric_limits<int>::max();  // UINT_E_MAX is reserved
        trussness_multi_capped[id] = std::numeric_limits<int>::max();
      });
      continue;
    }
    iter++;

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
      else if (still_active[a] == 1) return k;
      auto truss = trussness_multi_capped[a];
      if (truss == std::numeric_limits<int>::max()) return 0;
      //if (truss != UINT_E_MAX) return truss + 1;
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
          GA, id, u, v, decr_tab, get_trussness_and_id, k, false, inline_hierarchy, to_link, still_active);
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
      uintE new_deg = std::max((bucket_t) current_deg - triangles_removed, (bucket_t) lower_bound);
      std::get<1>(trussness_multi.big_table[id]) = current_deg - triangles_removed; //new_deg;  // update
      uintE old_deg = trussness_multi_capped[id];
      uintE new_bkt = std::max((bucket_t) get_bucket(new_deg),(bucket_t) k);
      trussness_multi_capped[id] = new_bkt;
      std::get<1>(decr_edges[i]) = b.get_bucket(old_deg, new_bkt);
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

    // Unmark edges removed in this round, and decrement their trussness.
    parallel_for(0, rem_edges.size(), [&](size_t i) {
      edge_t id = rem_edges[i];
      still_active[id] = 2;
    });

    // Clear the table storing the edge decrements.
    decr_tab.clear_table();
    prev_bkt = k;

    rounds++;
    cur_inner_rounds++;
  
  }


  double tt2 = t2.stop();
  std::cout << "### Real Peel Running Time: " << tt2 << std::endl;

  peeling_t.stop();
  peeling_t.next("peeling time");
  ct.next("Compaction time");
  bt.next("Bucketing time");
  em_t.next("EdgeMap time");
  decrement_t.next("Decrement trussness time");

  std::cout.precision(17);
  std::cout << "rho: " << iter << std::endl;
  std::cout << "clique core: " << k_max << std::endl;
  if (use_max_density) std::cout << "max density: " << max_density << std::endl;

  if (approx_out_str != "") {
    std::cout << approx_out_str << std::endl;
    std::cout << "Printing" << std::endl;
    std::ofstream file{approx_out_str};

    for (size_t i = 0; i < GA.n; i++) {
      auto f = [&](const uintE& u, const uintE& v, const W& wgh){
        auto truss_tuple = get_trussness_and_id(u, v);
        file << u << ", " << v << ": " << std::get<0>(truss_tuple) << std::endl;
      };
      GA.get_vertex(i).out_neighbors().map(f, false);
    }
    file.close();

    std::cout << "Finished printing" << std::endl;
  }

  return trussness_multi;
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

  //auto del_edges = gbbs::dyn_arr<edge_t>(6 * GA.n);
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
      else if (still_active[a] == 1) return k;
      auto truss = std::get<1>(trussness_multi.big_table[a]);
      if (truss == std::numeric_limits<int>::max()) return 0;
      if (truss != UINT_E_MAX) return truss + 1;
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

    // Unmark edges removed in this round, and decrement their trussness.
    parallel_for(0, rem_edges.size(), [&](size_t i) {
      edge_t id = rem_edges[i];
      still_active[id] = 2;
      std::get<1>(trussness_multi.big_table[id]) -= 1;
    });

    // Clear the table storing the edge decrements.
    decr_tab.clear_table();
    iter++;
    
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
void KTruss_connect(Graph& GA, size_t num_buckets, bool inline_hierarchy, bool efficient_inline_hierarchy, bool use_approx, std::string& approx_out_str, double delta) {
    if (efficient_inline_hierarchy) inline_hierarchy = true;

  truss_utils::multi_table<uintE, uintE, std::function<size_t(size_t)>> multitable;
  
  EfficientConnectWhilePeeling ecwp;
  ConnectWhilePeeling connect_with_peeling;

  if (!use_approx) {
    if (!efficient_inline_hierarchy) multitable = KTruss_ht(GA, connect_with_peeling, num_buckets, inline_hierarchy, false, false);
    else multitable = KTruss_ht(GA, ecwp, num_buckets, inline_hierarchy, false, false);
  } else {
    if (!efficient_inline_hierarchy) multitable = KTruss_approx(GA, connect_with_peeling, inline_hierarchy, approx_out_str, num_buckets, delta, false);
    else multitable = KTruss_approx(GA, ecwp, inline_hierarchy, approx_out_str, num_buckets, delta, false);
  }

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

} // namespace ktruss

}  // namespace gbbs
