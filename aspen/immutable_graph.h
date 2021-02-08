#pragma once

#include "build.h"

namespace aspen {


template <class weight>
struct symmetric_graph {

  // using weight = empty;  // placeholder for now.

  struct edge_entry {
    using key_t = vertex_id;   // a vertex_id
    using val_t = weight;        // placeholder
    static inline bool comp(key_t a, key_t b) {return a < b;}
  };

#ifdef USE_PAM
  using edge_tree = pam_map<edge_entry>;
#else
  // todo: use diffencoding.
  //using edge_tree = cpam::pam_map<edge_entry>;
  using edge_tree = cpam::diff_encoded_map<edge_entry>;
#endif

  struct vertex_entry {
    using key_t = vertex_id;
    using val_t = edge_tree;
    using aug_t = edge_id;
    static inline bool comp(key_t a, key_t b) {return a < b;}
    static aug_t get_empty() {return 0;}
    static aug_t from_entry(const key_t& k, const val_t& v) {
      return v.size();}
    static aug_t combine(aug_t a, aug_t b) {return a + b;}
  };

  using vertex_tree = aug_map<vertex_entry>;


  using ngh_and_weight = std::tuple<vertex_id, weight>;
  using edge = std::pair<vertex_id, ngh_and_weight>;
  using G = symmetric_graph;


  vertex_tree V;

  // Build from a static graph (todo).
  template <class Graph>
  symmetric_graph(Graph& G) {
    using W = typename Graph::weight_type;
    static_assert(std::is_same<W, weight>());
    auto edges = build::graph_to_edges(G);
    std::cout << "Edges.size = " << edges.size() << std::endl;
    V = from_edges(edges);
  }

  // Build from a sequence of edges.
  symmetric_graph(parlay::sequence<edge>& edges) {
    V = from_edges(edges);
  }

  static vertex_tree from_edges(parlay::sequence<edge>& edges) {
    auto reduce = [&] (parlay::slice<ngh_and_weight*,ngh_and_weight*> R) {
      return edge_tree(R);
    };
    vertex_tree vertices;
    return vertex_tree::multi_insert_reduce(vertices, edges, reduce);
  }

  // Reserve space for n vertices and m edges.
  void reserve(size_t n, size_t m) {
    // todo
  }

  auto get_vertex(vertex_id v) {
    // todo
  }

  auto edge_exists(edge e) {
    // todo
  }

  void print_stats() {
//    auto count_edges_fn = [&] (const auto& et) {
////      return 0;
//      //const auto& incident = std::get<1>(std::get<0>(et));
//      const auto& incident = std::get<1>(et);
//      return incident.size();
//    };
//    size_t edge_size = V.size_in_bytes(count_edges_fn);
    size_t sz = 0;
    size_t edges_bytes = 0;
    auto f = [&] (const auto& et) {
      const auto& incident = std::get<1>(et);
      auto noop = [] (const auto& q) { return 0; };
      size_t edges_size = incident.size();
      edges_bytes += incident.size_in_bytes(noop);
      if (edges_size < 2*cpam::utils::compression_block_size) {
	assert(incident.root_is_compressed());
      }
    };
    vertex_tree::foreach_seq(V, f);
    std::cout << "num_edges = " << sz << std::endl;
    std::cout << "edges_size = " << edges_bytes << std::endl;
//    auto count_vertex_fn = [&] (const auto& et) {
//      return 0;
//    };
//    size_t vertex_size = V.size_in_bytes(count_vertex_fn);
//    edge_size -= vertex_size;
//    std::cout << "Vertex tree requires " << vertex_size << " bytes." << std::endl;
//    std::cout << "Edge trees require " << edge_size << " bytes." << std::endl;
  }

};



}  // namespace aspen
