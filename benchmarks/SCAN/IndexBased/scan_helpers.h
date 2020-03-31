// Helper functions for SCAN logic.
//
// Mainly contains template functions that would otherwise clutter about the
// main SCAN header file.
#pragma once

#include <array>
#include <atomic>
#include <utility>

#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"
#include "ligra/graph.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/undirected_edge.h"
#include "pbbslib/get_time.h"

namespace indexed_scan {

namespace internal {

using StructuralSimilarities =
  sparse_table<UndirectedEdge, float, std::hash<UndirectedEdge>>;
using VertexSet =
  sparse_table<uintE, pbbslib::empty, decltype(&pbbslib::hash64_2)>;
using NoWeight = pbbslib::empty;

struct NeighborSimilarity {
  // Vertex ID.
  uintE neighbor;
  // Similarity of neighbor vertex to some original reference vertex.
  float similarity;
};
// Beware that this equality operator compares the floating-point field
// approximately. This is convenient for unit tests but might not be appropriate
// for other uses.
bool operator==(const NeighborSimilarity&, const NeighborSimilarity&);
std::ostream& operator<<(std::ostream& os, const NeighborSimilarity&);

using NeighborOrder = pbbs::sequence<pbbs::sequence<NeighborSimilarity>>;

struct CoreThreshold {
  uintE vertex_id;
  // Maximum value of the SCAN parameter epsilon for which `vertex_id` is a core
  // vertex (given some fixed reference value for SCAN parameter mu).
  float threshold;
};
// Beware that this equality operator compares the floating-point field
// approximately. This is convenient for unit tests but might not be appropriate
// for other uses.
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

// Prints the total time captured by `timer` to stderr if macro
// SCAN_DETAILED_TIMES is defined, otherwise does nothing.
void ReportTime(const timer&);

// Creates a `VertexSet` for holding up to `capacity` elements.
VertexSet MakeVertexSet(const size_t capacity);

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

  return lo +
    pbbs::binary_search(sequence.slice(lo, hi), std::forward<Func>(predicate));
}

// Computes an adjacency list for the graph in which each neighbor list is
// sorted by descending structural similarity with the source vertex.
//
// The output adjacency list `NO` is such that `NO[v][i]` is a pair `{u, sigma}`
// where `sigma` is the structural similarity between `v` and `u` and where `u`
// is the neighbor of `v` with the (zero-indexed) `i`-th  highest structural
// similarity with `v`.
//
// The structural similarity between two vertices u and v is
//   (size of intersection of closed neighborhoods of u and v) /
//   (geometric mean of size of closed neighborhoods of u and of v)
// where the closed neighborhood of a vertex x consists of all neighbors of x
// along with x itself.
//
// Unlike the presentation in "Efficient Structural Graph Clustering:  An
// Index-Based Approach", the neighbor list for a vertex `v` will not contain
// `v` itself, unless `(v, v)` is explicitly given as an edge in `graph`.
template <template <typename> class VertexTemplate>
NeighborOrder ComputeNeighborOrder(
    symmetric_graph<VertexTemplate, NoWeight>* graph) {
  using Vertex = VertexTemplate<NoWeight>;

  timer shared_neighbors_timer{"Compute shared neighbors time"};

  pbbs::sequence<std::atomic<uintE>> counters(
      graph->m, [](size_t) { return std::atomic<uintE>{0}; });
  pbbs::sequence<uintT> offsets(
      graph->n,
      [&](const size_t i) { return graph->get_vertex(i).getOutDegree(); });
  pbbslib::scan_add_inplace(offsets);
  // Maps an edge {u, v} to an index into `counters`. The counter will hold the
  // number of shared neighbors between u and v, i.e. the number of triangles
  // involving edge {u, v}.
  sparse_table<UndirectedEdge, uintT, std::hash<UndirectedEdge>>
    shared_neighbor_counters{
      graph->m / 2,
      std::make_pair(UndirectedEdge{UINT_E_MAX, UINT_E_MAX}, UINT_T_MAX),
      std::hash<UndirectedEdge>{}};
  par_for(0, graph->n, [&](const size_t v_id) {
    Vertex v{graph->get_vertex(v_id)};
    const auto initialize_counter{[&](
        const uintE source_vertex_id,
        const uintE neighbor_vertex_id,
        const uintE neighbor_index,
        NoWeight) {
      if (source_vertex_id < neighbor_vertex_id) {
        const uintT counter_index{offsets[source_vertex_id] + neighbor_index};
        shared_neighbor_counters.insert(
          {UndirectedEdge{source_vertex_id, neighbor_vertex_id},
           counter_index});
      }
    }};
    v.mapOutNghWithIndex(v_id, initialize_counter);
  });

  const auto update_shared_neighbor_counts{[&](
      const uintE vertex_1,
      const uintE vertex_2,
      const uintE vertex_3) {
    constexpr uintT kNotFound{UINT_T_MAX};
    const std::array<UndirectedEdge, 3> triangle_edges{{
       {vertex_1, vertex_2}, {vertex_1, vertex_3}, {vertex_2, vertex_3}}};
    for (const auto& edge : triangle_edges) {
      const uintT counter_index{shared_neighbor_counters.find(edge, kNotFound)};
      counters[counter_index]++;
    }
  }};
  Triangle_degree_ordering(*graph, update_shared_neighbor_counts);

  internal::ReportTime(shared_neighbors_timer);
  timer neighbor_order_timer{"Construct neighbor order time"};

  NeighborOrder neighbor_order{
    graph->n,
    [&graph](const size_t i) {
      return pbbs::sequence<NeighborSimilarity>::no_init(
        graph->get_vertex(i).getOutDegree());
    }
  };
  par_for(0, graph->n, [&](const uintE v_id) {
    Vertex v{graph->get_vertex(v_id)};
    const float v_neighborhood_sqrt{sqrtf(v.getOutDegree() + 1)};
    auto* const v_order{&neighbor_order[v_id]};
    const auto update_v_order{[&](
        const uintE source_vertex_id,
        const uintE u_vertex_id,
        const uintE u_index,
        NoWeight) {
      // The shared neighbor counts we've computed use open neighborhoods --
      // since structural similarity uses closed neighborhoods, we need to
      // adjust the number and denominator a little.
      constexpr uintT kNotFound{UINT_T_MAX};
      const uintT counter_index{shared_neighbor_counters.find(
          UndirectedEdge{source_vertex_id, u_vertex_id}, kNotFound)};
      // SCAN structural similarities are defined using _closed_ neighborhoods,
      // hence the need to to adjust these values by `+ 1` and `+ 2`.
      const uintE num_shared_neighbors{counters[counter_index] + 2};
      const float u_neighborhood_sqrt{
        sqrtf(graph->get_vertex(u_vertex_id).getOutDegree() + 1)};
      (*v_order)[u_index] = NeighborSimilarity{
          .neighbor = u_vertex_id,
          .similarity =
            num_shared_neighbors / (u_neighborhood_sqrt * v_neighborhood_sqrt)
      };
    }};
    v.mapOutNghWithIndex(v_id, update_v_order);

    // Sort by descending structural similarity.
    const auto compare_similarities_descending{
      [](const NeighborSimilarity& a, const NeighborSimilarity& b) {
        return a.similarity > b.similarity;
      }};
    pbbs::sample_sort_inplace(
        v_order->slice(),
        compare_similarities_descending);
  });
  internal::ReportTime(neighbor_order_timer);
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
