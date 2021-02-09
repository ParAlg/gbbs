#pragma once

#include "build.h"

namespace aspen {

template <class weight>
struct symmetric_graph {

  struct edge_entry {
    using key_t = vertex_id;   // a vertex_id
    using val_t = weight;        // placeholder
    static inline bool comp(key_t a, key_t b) {return a < b;}
    using entry_t = std::tuple<key_t, val_t>;
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
    using entry_t = std::tuple<key_t, val_t>;
  };
  using vertex_tree = aug_map<vertex_entry>;

  struct neighbors {
    edge_tree& tree;
    vertex_id id;
    neighbors(vertex_id id, edge_tree& tree) : id(id), tree(tree) {}

    template <class F, class G>
    void copy(size_t offset, F& f, G& g) {
      auto map_f = [&] (const auto& et, size_t i) {
        auto [ngh, wgh] = et;
        auto val = f(id, ngh, wgh);
        g(ngh, offset + i, val);
      };
      tree.foreach_index(tree, map_f);
    }

    template <class F>
    void map_index(F& f) {
      auto map_f = [&] (const auto& et, size_t i) {
        auto [ngh, wgh] = et;
        f(id, ngh, wgh, i);
      };
      tree.foreach_index(tree, map_f);
    }

    template <class F>
    void map(F& f) {
      auto map_f = [&] (const auto& et, size_t i) {
        auto [ngh, wgh] = et;
        auto val = f(id, ngh, wgh);
      };
      tree.foreach_index(tree, map_f);
    }

    template <class F>
    void foreach_cond(F& f) {
      auto map_f = [&] (const auto& et, size_t i) -> bool {
        auto [ngh, wgh] = et;
        return f(id, ngh, wgh);
      };
      tree.foreach_cond(tree, map_f);
    }

  };

  struct vertex {
    vertex_id id;
    edge_tree tree;
    size_t out_degree() { return tree.size(); }
    size_t in_degree() { return tree.size(); }
    auto out_neighbors() { return neighbors(id, tree); }
    auto in_neighbors() { return neighbors(id, tree); }
    vertex(vertex_id id, edge_tree&& tree) : id(id), tree(tree) {}
    // todo: map, etc, necessary for edgeMap.
  };

  using ngh_and_weight = std::tuple<vertex_id, weight>;
  using edge = std::pair<vertex_id, ngh_and_weight>;
  using maybe_vertex = std::optional<vertex>;
  using weight_type = weight;
  using G = symmetric_graph<weight>;

  vertex_tree V;

  // Build from a static graph.
  template <class Graph>
  symmetric_graph(Graph& GA) {
    using W = typename Graph::weight_type;
    static_assert(std::is_same<W, weight>());
    auto edges = build::graph_to_edges(GA);
    std::cout << "Edges.size = " << edges.size() << std::endl;
    G::reserve(GA.n, GA.m);
    V = from_edges(edges);
  }

  // Build from a sequence of edges.
  symmetric_graph(parlay::sequence<edge>& edges) {
    V = from_edges(edges);
  }

  vertex_tree& get_vertices() {
    return V;
  }

  size_t num_vertices() {
    return V.size();
  }

  size_t num_edges() {
    return V.aug_val();
  }

  vertex get_vertex(vertex_id v) {
    auto opt = V.find(v);
    if (opt.has_value()) {
      return vertex(v, std::move(*opt));
    }
    return vertex(v, edge_tree());
  }

  template <class F>
  void map_vertices(F& f) {
    using entry_t = typename vertex_entry::entry_t;
    auto map_f = [&] (entry_t vtx_entry, size_t i) {
      vertex_id v = std::get<0>(vtx_entry);
      auto vtx = vertex(v, std::move(std::get<1>(vtx_entry)));
      f(vtx);
    };
    V.foreach_index(V, map_f, 0, 1);
  }


  auto edge_exists(edge e) {
    // todo
  }

  static vertex_tree from_edges(parlay::sequence<edge>& edges) {
    auto reduce = [&] (parlay::slice<ngh_and_weight*,ngh_and_weight*> R) {
      return edge_tree(R);
    };
    vertex_tree vertices;
    return vertex_tree::multi_insert_reduce(vertices, edges, reduce);
  }

  // Reserve space for n vertices and m edges.
  static void reserve(size_t n, size_t m) {
    vertex_tree::reserve(n);
    edge_tree::reserve(m);
  }


  void print_stats() {
#ifndef USE_PAM
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
#else
#endif
  }

};



}  // namespace aspen
