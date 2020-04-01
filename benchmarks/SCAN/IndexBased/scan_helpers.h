// Helper functions for SCAN logic.
//
// Mainly contains template functions that would otherwise clutter about the
// main SCAN header file.
#pragma once

#include <atomic>
#include <functional>
#include <tuple>
#include <utility>

#include "benchmarks/SCAN/IndexBased/intersect.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"
#include "ligra/graph.h"
#include "ligra/graph_mutation.h"
#include "ligra/pbbslib/sparse_table.h"
#include "ligra/undirected_edge.h"
#include "pbbslib/get_time.h"
#include "pbbslib/sample_sort.h"

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

// For each edge {u, v}, counts the number of shared neighbors between u and
// v (where we don't include u`and v`themselves in the count).
//
// Arguments:
//   graph
//     Graph on which to count shared neighbors. The neighbor lists for each
//     vertex must be sorted by ascending neighbor ID.
//
// Returns:
//   `graph->m`-length sequence where each entry is (u, v, <number of shared
//   neighbors between u and v>) where u and v are adjacent. Ordering of entries
//   is arbitrary.
template <template <typename> class VertexTemplate>
pbbs::sequence<std::tuple<uintE, uintE, uintE>>
CountSharedNeighbors(symmetric_graph<VertexTemplate, NoWeight>* graph) {
  // Counting the neighbors shared between adjacent vertices u and v is the same
  // as counting the number of triangles that edge {u, v} appears in.
  //
  // The implementation is borrowed from `Triangle_degree_ordering()` in
  // `benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h` --- see the
  // associated paper "Multicore triangle computations without tuning."

  // Create a directed version of `graph`, pointing edges from lower degree
  // vertices to higher degree vertices in order to bound the maximum degree of
  // each vertex.
  uintE* vertex_degree_ranking = rankNodes(*graph, graph->n);
  const auto filter_predicate{[&](const uintE u, const uintE v, NoWeight) {
    return vertex_degree_ranking[u] < vertex_degree_ranking[v];
  }};
  auto directed_graph{filter_graph(*graph, filter_predicate)};

  // Estimate the amount of work for each vertex for better load balancing.
  auto parallel_work{sequence<size_t>(directed_graph.n)};
  {
    const auto get_edge_work{[&](const uintE u, const uintE v, NoWeight) {
      return directed_graph.get_vertex(v).getOutDegree();
    }};
    const auto addition{pbbslib::addm<size_t>()};
    par_for(0, directed_graph.n, [&](const size_t i) {
      parallel_work[i] =
        directed_graph.get_vertex(i).template reduceOutNgh<size_t>(
            i, get_edge_work, addition);
    });
  }
  const size_t total_work{pbbslib::scan_add_inplace(parallel_work)};
  constexpr size_t work_block_size{50000};
  const size_t num_blocks{(total_work + work_block_size - 1) / work_block_size};

  // Each counter in `counters` holds the number of shared neighbors in `graph`
  // between u and v for some edge {u, v}.
  pbbs::sequence<std::atomic<uintE>> counters(
      directed_graph.m, [](size_t) { return std::atomic<uintE>{0}; });
  // We use `counter_offsets` to be able to index into `counters` for each edge.
  pbbs::sequence<uintT> counter_offsets(
      directed_graph.n,
      [&](const size_t i) {
        return directed_graph.get_vertex(i).getOutDegree();
      });
  pbbslib::scan_add_inplace(counter_offsets);

  // Find triangles of the following form:
  //        w
  //       ^ ^
  //      /   \.
  //     u --> v
  // There's a bijection between triangles of this form in `directed_graph` and
  // undirected triangles in `graph`.
  const auto run_intersection_on_block{[&]
      (const size_t block_start, const size_t block_end) {
    for (size_t vertex_id = block_start; vertex_id < block_end; vertex_id++) {
      auto vertex{directed_graph.get_vertex(vertex_id)};
      const uintT vertex_counter_offset{counter_offsets[vertex_id]};
      const auto intersect{[&](
          const uintE v_id,
          const uintE neighbor_id,
          const uintE v_to_neighbor_index,
          NoWeight) {
        auto neighbor{directed_graph.get_vertex(neighbor_id)};
        const uintT neighbor_counter_offset{counter_offsets[neighbor_id]};
        const auto update_counters{[&](
            const uintE shared_neighbor,
            const uintE vertex_to_shared_index,
            const uintE neighbor_to_shared_index) {
          counters[vertex_counter_offset + vertex_to_shared_index]++;
          counters[neighbor_counter_offset + neighbor_to_shared_index]++;
        }};
        counters[vertex_counter_offset + v_to_neighbor_index] +=
          intersect_f_with_index_par(
              &vertex, &neighbor, vertex_id, neighbor_id, update_counters);
      }};
      constexpr bool kParallel{false};
      vertex.mapOutNghWithIndex(vertex_id, intersect, kParallel);
    }
  }};
  constexpr size_t kGranularity{1};
  par_for(0, num_blocks, kGranularity, [&](size_t i) {
    const size_t work_start{i * work_block_size};
    const size_t work_end{work_start + work_block_size};
    constexpr auto less_fn{std::less<size_t>()};
    const size_t block_start{
      pbbs::binary_search(parallel_work, work_start, less_fn)};
    const size_t block_end{
      pbbs::binary_search(parallel_work, work_end, less_fn)};
    run_intersection_on_block(block_start, block_end);
  });

  pbbs::sequence<std::tuple<uintE, uintE, uintE>> shared_neighbor_counts(
      graph->m);
  par_for(0, directed_graph.n, [&](const size_t vertex_id) {
    auto vertex{directed_graph.get_vertex(vertex_id)};
    const uintT vertex_counter_offset{counter_offsets[vertex_id]};
    const auto create_edge_counts{[&](
        const uintE v_id,
        const uintE neighbor_id,
        const uintE v_to_neighbor_index,
        NoWeight) {
      const uintT counter_index{vertex_counter_offset + v_to_neighbor_index};
      const uintE shared_neighbor_count{counters[counter_index]};
      shared_neighbor_counts[2 * counter_index] =
        std::make_tuple(v_id, neighbor_id, shared_neighbor_count);
      shared_neighbor_counts[2 * counter_index + 1] =
        std::make_tuple(neighbor_id, v_id, shared_neighbor_count);
    }};
    vertex.mapOutNghWithIndex(vertex_id, create_edge_counts);
  });

  directed_graph.del();
  pbbs::free_array(vertex_degree_ranking);
  return shared_neighbor_counts;
}

// Computes an adjacency list for the graph in which each neighbor list is
// sorted by descending structural similarity with the source vertex.
//
// The neighbor lists for each vertex in the graph must be sorted by ascending
// neighbor ID.
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

  timer shared_neighbors_timer{
    "Compute neighbor order - count shared neighbors time"};
  pbbs::sequence<std::tuple<uintE, uintE, uintE>> shared_neighbor_counts{
    CountSharedNeighbors(graph)};
  internal::ReportTime(shared_neighbors_timer);
  timer neighbor_order_timer{"Compute neighbor order - construct order time"};

  // Assuming that `graph`'s neighbor lists are sorted by ascending neighbor
  // IDs, a sorted `shared_neighbor_counts` will match up nicely with `graph`'s
  // neighbor lists to allow us to retrieve shared neighbor counts for each
  // edge.
  pbbs::sample_sort_inplace(
      shared_neighbor_counts.slice(),
      [](const std::tuple<uintE, uintE, uintE>& a,
         const std::tuple<uintE, uintE, uintE>& b) {
        return a < b;
      });
  pbbs::sequence<uintT> count_offsets(
      graph->n,
      [&](const size_t i) { return graph->get_vertex(i).getOutDegree(); });
  pbbslib::scan_add_inplace(count_offsets);

  NeighborOrder neighbor_order{
    graph->n,
    [&graph](const size_t i) {
      return pbbs::sequence<NeighborSimilarity>::no_init(
        graph->get_vertex(i).getOutDegree());
    }
  };
  par_for(0, graph->n, [&](const uintE v_id) {
    Vertex v{graph->get_vertex(v_id)};
    const uintT v_count_offset{count_offsets[v_id]};
    const float v_neighborhood_sqrt{sqrtf(v.getOutDegree() + 1)};
    auto* const v_order{&neighbor_order[v_id]};
    const auto update_v_order{[&](
        const uintE _v_id,
        const uintE u_vertex_id,
        const uintE v_to_u_index,
        NoWeight) {
      // SCAN structural similarities are defined using _closed_ neighborhoods,
      // hence the need to to adjust these values by `+ 1` and `+ 2`.
      const uintE num_shared_neighbors{
        std::get<2>(shared_neighbor_counts[v_count_offset + v_to_u_index]) + 2};
      const float u_neighborhood_sqrt{
        sqrtf(graph->get_vertex(u_vertex_id).getOutDegree() + 1)};
      (*v_order)[v_to_u_index] = NeighborSimilarity{
          .neighbor = u_vertex_id,
          .similarity =
            num_shared_neighbors / (v_neighborhood_sqrt * u_neighborhood_sqrt)
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
