#pragma once

#include "gbbs/vertex_subset.h"
#include "gbbs/flags.h"

namespace aspen {

// Defines edge_map, vertex_map, adds high level traversal primitives over the
// underlying graph
template <class graph>
struct traversable_graph : private graph {
  using G = graph;
  using G::G; // imports graph's constructors

  using vertex_tree = typename G::vertex_tree;
  using edge_tree = typename G::edge_tree;
  using vertex = typename G::vertex;
  using weight_type = typename graph::weight_type;
  using Empty = gbbs::empty;

  template <class Data>
  using vertexSubsetData = gbbs::vertexSubsetData<Data>;

  using flags = gbbs::flags;

  template <class Data, class VS, class F>
  auto edgeMapSparse(VS& vs, const F& f, const flags& fl) {
    using S = typename vertexSubsetData<Data>::S;
    size_t n = num_vertices();
    assert(gbbs::should_output(fl));

    auto offsets = parlay::sequence<edge_id>::from_function(
        vs.size(), [&](size_t i) {
      return (fl & gbbs::in_edges) ? get_vertex(vs.vtx(i)).in_degree()
                             : get_vertex(vs.vtx(i)).out_degree();
    });
    size_t outEdgeCount =
      parlay::scan_inplace(parlay::make_slice(offsets));
    auto outEdges = parlay::sequence<S>::uninitialized(outEdgeCount);
    std::cout << "OutEdgeCount = " << outEdgeCount << std::endl;
    parlay::parallel_for(0, vs.size(), [&] (size_t i) {
      vertex_id v = vs.vtx(i);
      edge_id o = offsets[i];
      auto vtx = get_vertex(v);
      auto& neighbors = (fl & gbbs::in_edges) ? vtx.out_neighbors() : vtx.in_neighbors();
      S* out_edges = outEdges.begin();
      auto g = gbbs::get_emsparse_gen_full<Data>(out_edges);
      auto map_f = [&] (const auto& et, size_t i) {
	auto m = f.updateAtomic(v, std::get<0>(et), std::get<1>(et));
	g(std::get<0>(et), o + i, m);
      };
      neighbors.foreach_index(neighbors, map_f);
    });
    auto filtered = parlay::filter(outEdges, [&] (const auto& e) {
      return e != UINT_E_MAX;
    });
    std::cout << "filtered size = " << filtered.size() << std::endl;
    return vertexSubsetData<Data>(n, std::move(filtered));
  }

  template <class Data, class VS, class F>
  auto edgeMapSparse(VS& vs, const F& f, const flags& fl) {
    using D = typename vertexSubsetData<Data>::D;
    size_t n = GA.n;
    assert(gbbs::should_output(fl));

  }

  template <class Data, class VS, class F>
  auto edgeMapData(VS& vs, const F& f, long threshold = -1, const flags& fl = 0) {
    size_t n = num_vertices();
    size_t m = num_edges();
    if (threshold == -1) threshold = n / 20;
    if (vs.size() == 0) return vertexSubsetData<Data>(n);

    size_t out_degrees = 0;
    if (vs.out_degrees_set()) {
      out_degrees = vs.get_out_degrees();
    } else {
      vs.toSparse();
      auto degree_f = [&](size_t i) {
        auto vertex = get_vertex(vs.vtx(i));
        return (fl & gbbs::in_edges) ? vertex.in_degree() : vertex.out_degree();
      };
      auto degree_im = parlay::delayed_seq<size_t>(vs.size(), degree_f);
      out_degrees = pbbslib::reduce_add(degree_im);
      vs.set_out_degrees(out_degrees);
    }
    if (out_degrees == 0) return vertexSubsetData<Data>(n);
    if (vs.size() + out_degrees > threshold && !(fl & gbbs::no_dense)) {
      return edgeMapDense(vs, f, fl);
    }
    return edgeMapSparse<Data, VS, F>(vs, f, fl);
  }

  template <class VS, class F>
  auto edgeMap(VS& vs, const F& f, long threshold = -1, const flags& fl = 0) {
    return edgeMapData<Empty, VS, F>(vs, f, threshold, fl);
  }


  using G::print_stats;
  using G::num_vertices;
  using G::num_edges;
  using G::get_vertex;
};


}  // namespace aspen
