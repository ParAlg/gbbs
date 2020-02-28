// Helper functions for SCAN logic.
//
// Mainly contains template functions that would otherwise clutter about the
// main SCAN header file.
#pragma once

#include <utility>

#include "ligra/graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/undirected_edge.h"

namespace indexed_scan {

namespace internal {

using StructuralSimilarities =
  sparse_table<UndirectedEdge, float, std::hash<UndirectedEdge>>;
using VertexSet =
  sparse_table<uintE, pbbslib::empty, decltype(&pbbslib::hash64_2)>;
using Weight = pbbslib::empty;

struct NeighborSimilarity {
  // Vertex ID.
  uintE neighbor;
  // Similarity of neighbor vertex to some original reference vertex.
  float similarity;
};
bool operator==(const NeighborSimilarity&, const NeighborSimilarity&);
std::ostream& operator<<(std::ostream& os, const NeighborSimilarity&);

using NeighborOrder = pbbs::sequence<pbbs::sequence<NeighborSimilarity>>;

struct CoreThreshold {
  uintE vertex_id;
  // Maximum value of the SCAN parameter epsilon for which `vertex_id` is a core
  // vertex (given some fixed reference value for SCAN parameter mu).
  float threshold;
};
bool operator==(const CoreThreshold&, const CoreThreshold&);
std::ostream& operator<<(std::ostream& os, const CoreThreshold&);

class CoreOrder {
 public:
  explicit CoreOrder(const NeighborOrder& neighbor_order);

  // Return all vertices that are cores under SCAN parameters `mu` and
  // `epsilon`.
  pbbs::sequence<uintE> GetCores(uint64_t mu, float epsilon) const;

 private:
  const size_t num_vertices_;
  pbbs::sequence<pbbs::sequence<CoreThreshold>> order_{};
};

// Creates a `VertexSet` for holding up to `capacity` elements.
inline VertexSet MakeVertexSet(const size_t capacity) {
  return
    make_sparse_table<uintE, pbbslib::empty, decltype(&pbbslib::hash64_2)>(
      // Adding 1 avoids having small tables completely full.
      capacity + 1, {UINT_E_MAX, pbbslib::empty{}}, pbbslib::hash64_2);
}

// Finds the least `i` such that `predicate(sequence[i])` is false. If
// `predicate(sequence[i])` is true for all `i`, then this returns
// `sequence.size()`.
//
// `sequence` and `predicate` must be partitioned such that there is some `i`
// (which will be the return value) for which `predicate(sequence[i])` is true
// for all j < i and for which `predicate(sequence[i])` is false for all j >= i.
//
// `predicate` should take a `SeqElement` and return a boolean.
//
// Running time is O(log [return value]).
template <class SeqElement, class Func>
size_t BinarySearch(
    const pbbs::sequence<SeqElement>& sequence,
    Func&& predicate) {
  // Start off the binary search with an exponential search so that running time
  // while be O(log [return value]) rather than O(log sequence.size()).
  // In practice this probably doesn't matter much....
  size_t hi{0};
  while (hi < sequence.size() && predicate(sequence[hi])) {
    hi = 2 * hi + 1;
  }
  const size_t lo{hi / 2};
  if (hi > sequence.size()) {
    hi = sequence.size();
  }

  return pbbs::binary_search(sequence.slice(lo, hi), predicate);
}

// Compute structural similarities (as defined by SCAN) between each pair of
// adjacent vertices.
//
// The structural similarity between two vertices u and v is
//   (size of intersection of closed neighborhoods of u and v) /
//   (geometric mean of size of closed neighborhoods of u and of v)
// where the closed neighborhood of a vertex x consists of all neighbors of x
// along with x itself.
template <template <typename WeightType> class VertexType>
StructuralSimilarities ComputeStructuralSimilarities(
    symmetric_graph<VertexType, Weight>* graph) {
  using Vertex = VertexType<Weight>;

  StructuralSimilarities similarities{
    graph->m,
    std::make_pair(UndirectedEdge{UINT_E_MAX, UINT_E_MAX}, 0.0),
    std::hash<UndirectedEdge>{}};

  pbbs::sequence<VertexSet> adjacency_list{
    pbbs::sequence<VertexSet>::no_init(graph->n)};
  parallel_for(0, graph->n, [&graph, &adjacency_list](const size_t vertex_id) {
    Vertex vertex{graph->get_vertex(vertex_id)};
    auto* const neighbors{&adjacency_list[vertex_id]};
    *neighbors = MakeVertexSet(vertex.getOutDegree());

    const auto update_adjacency_list{[&neighbors](
        const uintE source_vertex,
        const uintE neighbor_vertex,
        const Weight weight) {
      neighbors->insert({neighbor_vertex, pbbslib::empty{}});
    }};
    vertex.mapOutNgh(vertex_id, update_adjacency_list);
  });

  graph->map_edges([&graph, &adjacency_list, &similarities](
        const uintE u_id, const uintE v_id, const Weight) {
      // Only perform this computation once for each undirected edge
      if (u_id < v_id) {
        Vertex u{graph->get_vertex(u_id)};
        Vertex v{graph->get_vertex(v_id)};
        const auto& u_neighbors{adjacency_list[u_id]};
        const auto& v_neighbors{adjacency_list[v_id]};

        const bool u_degree_is_smaller{u.getOutDegree() < v.getOutDegree()};
        const uintE smaller_degree_vertex_id{u_degree_is_smaller ? u_id : v_id};
        Vertex* smaller_degree_vertex{u_degree_is_smaller ? &u : &v};
        const auto& larger_degree_vertex_neighbors{
          u_degree_is_smaller ? v_neighbors : u_neighbors
        };

        const auto is_shared_neighbor{
          [&](const uintE, const uintE neighbor, const Weight) {
            return larger_degree_vertex_neighbors.contains(neighbor);
          }};
        const auto add_monoid{pbbslib::addm<size_t>()};
        const size_t num_shared_neighbors{
          smaller_degree_vertex->template reduceOutNgh<size_t>(
              smaller_degree_vertex_id,
              is_shared_neighbor,
              add_monoid)};

        // The neighborhoods we've computed are open neighborhoods -- since
        // structural similarity uses closed neighborhoods, we need to adjust
        // the number and denominator a little.
        similarities.insert(
            {UndirectedEdge{u_id, v_id},
            (num_shared_neighbors + 2) /
                (sqrt(u.getOutDegree() + 1) * sqrt(v.getOutDegree() + 1))});
      }
  });

  return similarities;
}

// Computes an adjacency list for the graph in which each neighbor list is
// sorted by descending structural similarity with the source vertex.
//
// The output adjacency list `NO` is such that `NO[v][i]` is a pair `{u, sigma}`
// where `sigma` is the structural similarity between `v` and `u` and where `u`
// is the neighbor of `v` with the (zero-indexed) `i`-th  highest structural
// similarity with `v`.
//
// Unlike the presentation in "Efficient Structural Graph Clustering:  An
// Index-Based Approach", the neighbor list for a vertex `v` will not contain
// `v` itself, unless `(v, v)` is explicitly given as an edge in `graph`.
template <template <typename WeightType> class VertexType>
NeighborOrder ComputeNeighborOrder(
    symmetric_graph<VertexType, Weight>* graph,
    const StructuralSimilarities& similarities) {
  using Vertex = VertexType<Weight>;

  NeighborOrder neighbor_order{
    graph->n,
    [&graph](const size_t i) {
      return pbbs::sequence<NeighborSimilarity>::no_init(
        graph->get_vertex(i).getOutDegree());
    }
  };
  par_for(0, graph->n, [&](const uintE v) {
    Vertex vertex{graph->get_vertex(v)};
    auto* const v_order{&neighbor_order[v]};

    par_for(0, vertex.getOutDegree(), [&](const size_t i) {
      const uintE neighbor{vertex.getOutNeighbor(i)};
      constexpr float kNotFound{-1.0};
      const float similarity{
        similarities.find(UndirectedEdge{v, neighbor}, kNotFound)};
      (*v_order)[i] = NeighborSimilarity{
          .neighbor = neighbor, .similarity = similarity};
    });

    // Sort by descending structural similarity.
    const auto compare_similarities_descending{
      [](const NeighborSimilarity& a, const NeighborSimilarity& b) {
        return a.similarity > b.similarity;
      }};
    pbbs::sample_sort_inplace(
        v_order->slice(),
        compare_similarities_descending);
  });
  return neighbor_order;
}

// Returns a sequence CO where CO[i] for i >= 2 is a list of vertices that
// can be a core when the SCAN parameter mu is set to i. The vertices in
// CO[i] are sorted by their core threshold values, the maximum value of SCAN
// parameter epsilon such that the vertex is a core when mu == i.
//
// CO[0] and CO[1] are left empty --- when mu is less than 2, all vertices are
// always cores and have a core threshold of 1.
pbbs::sequence<pbbs::sequence<CoreThreshold>> ComputeCoreOrder(
    const NeighborOrder& neighbor_order);

}  // namespace internal

}  // namespace indexed_scan
