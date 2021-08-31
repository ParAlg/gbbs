#include "edge_array.h"
#include "edge_map_data.h"
#include "edge_map_reduce.h"
#include "flags.h"
#include "graph.h"
#include "graph_mutation.h"
#include "macros.h"
#include "vertex_subset.h"

namespace gbbs {

// =================== VertexSubset Operators: Mapping =====================
// Applies the the supplied (update, cond) operators on edges incident to the
// input vertex_subset, aggregating results at the neighbors of this vset.
// This is the specialized version of the general nghMap function where the
// output is a plain vertex_subset.
// F is a struct providing the following functions:
//   update : (uintE * uintE * W) -> bool
//   updateAtomic : (uintE * uintE * W) -> bool
//   cond : uintE -> bool
template <class Graph, class VS, class F>
inline vertexSubsetData<gbbs::empty> neighbor_map(Graph& G, VS& vs, F f,
                                                  intT threshold = -1,
                                                  flags fl = 0) {
  return edgeMapData<gbbs::empty>(G, vs, f, threshold, fl);
}

// Generalized version of edgeMap. Takes an input vertex_subset and
// aggregates results at the neighbors of the vset. See above for the
// requirements on F.
template <class Graph,  // underlying graph class
          class Data,   // data associated with vertices in the output
          class VS,     // input vertex_subset type
          class F>      // map function type
inline vertexSubsetData<Data>
nghMapData(Graph& G, VS& vs, F f, intT threshold = -1, flags fl = 0) {
  return edgeMapData<Data>(G, vs, f, threshold, fl);
}

// =================== VertexSubset Operators: Counting =====================
template <class Graph,
          class Cond,   // function from uintE -> bool
          class Apply,  // function from std::tuple<uintE, uintE> ->
                        // std::optional<std::tuple<uintE, Data>>
          class VS>
inline vertexSubsetData<uintE> nghCount(Graph& G, VS& vs, Cond cond_f,
                                        Apply apply_f,
                                        hist_table<uintE, uintE>& ht,
                                        flags fl = 0) {
  return edgeMapCount<uintE, Cond, Apply, VS>(G, vs, cond_f, apply_f, ht, fl);
}

template <class Graph, class Cond, class VS>
inline vertexSubsetData<uintE> nghCount(Graph& G, VS& vs, Cond cond_f,
                                        hist_table<uintE, uintE>& ht,
                                        flags fl = 0) {
  auto apply_f = [&](const std::tuple<uintE, uintE>& ct) {
    return std::optional<std::tuple<uintE, uintE>>(ct);
  };
  return edgeMapCount<uintE, Cond, decltype(apply_f), VS>(G, vs, cond_f,
                                                          apply_f, ht, fl);
}

template <class Graph, class P, class VS>
inline vertexSubsetData<uintE> srcCount(Graph& G, VS& vs, P p, flags fl = 0) {
  return edgeMapFilter(G, vs, p, fl);
}

// =================== VertexSubset Operators: Packing =====================
template <class Graph, class P>
vertexSubsetData<uintE> srcPack(Graph& G, vertexSubset& vs, P p, flags fl = 0) {
  return packEdges(G, vs, p, fl);
}

// ======================= Graph Operators: Filtering ========================
// Filters the symmetric graph, G, with a predicate function pred.  Note
// that the predicate does not have to be symmetric, i.e. f(u,v) is
// not necesssarily equal to f(v,u), but we only represent the out-edges of
// this (possibly) directed graph. For convenience in cases where the graph
// needed is symmetric, we coerce this to a symmetric_graph.
template <
    template <class inner_wgh> class vtx_type, class wgh_type, typename P,
    typename std::enable_if<
        std::is_same<vtx_type<wgh_type>, symmetric_vertex<wgh_type>>::value,
        int>::type = 0>
static inline symmetric_graph<symmetric_vertex, wgh_type> filterGraph(
    symmetric_graph<vtx_type, wgh_type>& G, P& pred) {
  auto ret = filter_graph<vtx_type, wgh_type>(G, pred);
  auto newN = std::get<0>(ret);
  auto newM = std::get<1>(ret);
  auto newVData = std::get<2>(ret);
  auto newEdges = std::get<3>(ret);

  assert(newN == G.num_vertices());
  return symmetric_graph<symmetric_vertex, wgh_type>(
      newVData, newN, newM,
      [=]() {
        gbbs::free_array(newVData, newN);
        gbbs::free_array(newEdges, newM);
      },
      newEdges);
}

template <
    template <class inner_wgh> class vtx_type, class wgh_type, typename P,
    typename std::enable_if<
        std::is_same<vtx_type<wgh_type>, csv_bytepd_amortized<wgh_type>>::value,
        int>::type = 0>
static inline symmetric_graph<csv_byte, wgh_type> filterGraph(
    symmetric_graph<vtx_type, wgh_type>& G, P& pred) {
  auto ret = filter_graph<vtx_type, wgh_type>(G, pred);
  auto newN = std::get<0>(ret);
  auto newM = std::get<1>(ret);
  auto newVData = std::get<2>(ret);
  auto newEdges = std::get<3>(ret);

  assert(newN == G.num_vertices());
  return symmetric_graph<csv_byte, wgh_type>(newVData, newN, newM,
                                             [=]() {
                                               gbbs::free_array(newVData, newN);
                                               gbbs::free_array(newEdges, newM);
                                             },
                                             newEdges);
}

// Used by MST and MaximalMatching
// Predicate returns three values:
// 0 : keep in graph, do not return in edge array
// 1 : remove from graph, do not return in edge array
// 2 : remove from graph, return in edge array
// Cost: O(n+m) work
template <class Graph, class P>
edge_array<typename Graph::weight_type> filterEdges(Graph& G, P& pred,
                                                    flags fl = 0) {
  return filter_edges(G, pred, fl);
}

template <class Graph, class P>
edge_array<typename Graph::weight_type> filterAllEdges(Graph& G, P& pred,
                                                       flags fl = 0) {
  return filter_all_edges(G, pred, fl);
}

template <class Graph, class P>
edge_array<typename Graph::weight_type> sampleEdges(Graph& G, P& pred) {
  return sample_edges(G, pred);
}

}  // namespace gbbs
