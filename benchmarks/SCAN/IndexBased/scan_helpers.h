// Helper functions for SCAN logic.
//
// Mainly contains template functions that would otherwise clutter about the
// main SCAN header file.
#pragma once

#include <atomic>
#include <functional>
#include <utility>

#include "benchmarks/SCAN/IndexBased/intersect.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"
#include "cereal/access.hpp"
#include "gbbs/graph.h"
#include "gbbs/graph_mutation.h"
#include "gbbs/undirected_edge.h"
#include "pbbslib/get_time.h"
#include "pbbslib/sample_sort.h"

namespace indexed_scan {

namespace internal {

using NoWeight = pbbslib::empty;

struct EdgeSimilarity {
  // Source vertex ID.
  uintE source;
  // Neighbor vertex ID.
  uintE neighbor;
  // Similarity of source vertex to neighbor vertex.
  float similarity;
};
bool operator==(const EdgeSimilarity&, const EdgeSimilarity&);
std::ostream& operator<<(std::ostream& os, const EdgeSimilarity&);

// An adjacency list for the graph in which each vertex's neighbor list is
// sorted by descending structural similarity.
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
class NeighborOrder {
 public:
  // Constructor.
  //
  // The neighbor lists for each vertex in the graph must be sorted by ascending
  // neighbor ID.
  template <template <typename> class VertexTemplate>
  explicit NeighborOrder(symmetric_graph<VertexTemplate, NoWeight>* graph);

  NeighborOrder();

  // Get all structural similarity scores from vertex `source` to its neighbors,
  // sorted by descending similarity.
  const pbbs::range<EdgeSimilarity*>& operator[](size_t source) const;

  bool empty() const;
  // Returns the number of vertices.
  size_t size() const;

  pbbs::range<EdgeSimilarity*>* begin() const;
  pbbs::range<EdgeSimilarity*>* end() const;

 private:
  friend class cereal::access;
  friend bool operator==(const NeighborOrder&, const NeighborOrder&);

  // For cereal serialization.
  template<class CerealArchive>
  void save(CerealArchive& archive) const;
  template<class CerealArchive>
  void load(CerealArchive& archive);

  // Holds similarity scores for all edges, sorted by source and then by
  // similarity.
  pbbs::sequence<EdgeSimilarity> similarities_;
  pbbs::sequence<pbbs::range<EdgeSimilarity*>> similarities_by_source_;
};

struct CoreThreshold {
  uintE vertex_id;
  // Maximum value of the SCAN parameter epsilon for which `vertex_id` is a core
  // vertex (given some fixed reference value for SCAN parameter mu).
  float threshold;

 private:
  friend class cereal::access;

  // For cereal serialization.
  template<class CerealArchive>
  void serialize(CerealArchive& archive) {
    archive(vertex_id, threshold);
  }
};
bool operator==(const CoreThreshold&, const CoreThreshold&);
std::ostream& operator<<(std::ostream& os, const CoreThreshold&);

class CoreOrder {
 public:
  explicit CoreOrder(const NeighborOrder& neighbor_order);

  CoreOrder();

  // Return all vertices that are cores under SCAN parameters `mu` and
  // `epsilon`.
  pbbs::sequence<uintE> GetCores(uint64_t mu, float epsilon) const;

 private:
  friend class cereal::access;
  friend bool operator==(const CoreOrder&, const CoreOrder&);

  // For cereal serialization.
  template<class CerealArchive>
  void save(CerealArchive& archive) const;
  template<class CerealArchive>
  void load(CerealArchive& archive);

  size_t num_vertices_;
  pbbs::sequence<pbbs::sequence<CoreThreshold>> order_{};
};

// Prints the total time captured by `timer` to stderr if macro
// SCAN_DETAILED_TIMES is defined, otherwise does nothing.
void ReportTime(const timer&);

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
template <class Seq, class Func>
size_t BinarySearch(const Seq& sequence, Func&& predicate) {
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

// Returns a sequence CO where CO[i] for i >= 2 is a list of vertices that
// can be a core when the SCAN parameter mu is set to i. The vertices in
// CO[i] are sorted by their core threshold values, the maximum value of SCAN
// parameter epsilon such that the vertex is a core when mu == i.
//
// CO[0] and CO[1] are left empty --- when mu is less than 2, all vertices are
// always cores and have a core threshold of 1.
pbbs::sequence<pbbs::sequence<CoreThreshold>> ComputeCoreOrder(
    const NeighborOrder& neighbor_order);

template <template <typename> class VertexTemplate>
NeighborOrder::NeighborOrder(
    symmetric_graph<VertexTemplate, NoWeight>* graph) {
  timer function_timer{"Construct neighbor order"};

  // To compute structural similarities, we need to count shared neighbors
  // between two vertices. Counting the neighbors shared between adjacent
  // vertices u and v is the same as counting the number of triangles that the
  // edge {u, v} appears in.
  //
  // The triangle counting logic here is borrowed from
  // `Triangle_degree_ordering()` in
  // `benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h`. We modify it to
  // maintain triangle counts for each edge.

  // Create a directed version of `graph`, pointing edges from lower degree
  // vertices to higher degree vertices in order to bound the maximum degree of
  // each vertex.
  uintE* vertex_degree_ranking{rankNodes(*graph, graph->n)};
  const auto filter_predicate{[&](const uintE u, const uintE v, NoWeight) {
    return vertex_degree_ranking[u] < vertex_degree_ranking[v];
  }};
  auto directed_graph{graph->filterGraph(*graph, filter_predicate)};

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
  const size_t num_work_blocks{
    (total_work + work_block_size - 1) / work_block_size};

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
          NoWeight,
          const uintE v_to_neighbor_index) {
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
  par_for(0, num_work_blocks, kGranularity, [&](size_t i) {
    const size_t work_start{i * work_block_size};
    const size_t work_end{work_start + work_block_size};
    constexpr auto less_fn{std::less<size_t>()};
    const size_t block_start{
      pbbs::binary_search(parallel_work, work_start, less_fn)};
    const size_t block_end{
      pbbs::binary_search(parallel_work, work_end, less_fn)};
    run_intersection_on_block(block_start, block_end);
  });

  // Convert shared neighbor counts into structural similarities for each edge.
  similarities_ = pbbs::sequence<EdgeSimilarity>(graph->m);
  par_for(0, directed_graph.n, [&](const size_t vertex_id) {
    const uintT v_counter_offset{counter_offsets[vertex_id]};
    const float v_neighborhood_sqrt{
      sqrtf(graph->get_vertex(vertex_id).getOutDegree() + 1)};
    const auto compute_similarity{[&](
        const uintE v_id,
        const uintE u_id,
        NoWeight,
        const uintE v_to_u_index) {
      // SCAN structural similarities are defined using _closed_ neighborhoods,
      // hence the need to to adjust these values by `+ 1` and `+ 2`.
      const float u_neighborhood_sqrt{
        sqrtf(graph->get_vertex(u_id).getOutDegree() + 1)};
      const uintT counter_index{v_counter_offset + v_to_u_index};
      const uintE num_shared_neighbors{counters[counter_index] + 2};
      const float structural_similarity{
        num_shared_neighbors / (v_neighborhood_sqrt * u_neighborhood_sqrt)};
      similarities_[2 * counter_index] =
        {.source = v_id, .neighbor = u_id, .similarity = structural_similarity};
      similarities_[2 * counter_index + 1] =
        {.source = u_id, .neighbor = v_id, .similarity = structural_similarity};
    }};
    directed_graph.get_vertex(vertex_id).mapOutNghWithIndex(
        vertex_id, compute_similarity);
  });

  pbbs::sample_sort_inplace(
      similarities_.slice(),
      [](const EdgeSimilarity& left, const EdgeSimilarity& right) {
        // Sort by ascending source, then descending similarity.
        return std::tie(left.source, right.similarity) <
          std::tie(right.source, left.similarity);
      });
  pbbs::sequence<uintT> source_offsets(
      graph->n,
      [&](const size_t i) { return graph->get_vertex(i).getOutDegree(); });
  pbbslib::scan_add_inplace(source_offsets);
  similarities_by_source_ = pbbs::sequence<pbbs::range<EdgeSimilarity*>>(
      graph->n,
      [&](const size_t i) {
        return similarities_.slice(
          source_offsets[i],
          i + 1 == graph->n ? similarities_.size() : source_offsets[i + 1]);
      });

  directed_graph.del();
  pbbs::free_array(vertex_degree_ranking);
  internal::ReportTime(function_timer);
}

template<class CerealArchive>
void NeighborOrder::save(CerealArchive& archive) const {
  // Store a neighbor order compactly by storing the offset at which each source
  // vertex first appears in `similarities_` and storing {neighbor, similarity}
  // for each edge in `similarities_`.

  const size_t num_vertices{similarities_by_source_.size()};
  const size_t num_edges{similarities_.size()};
  archive(num_vertices, num_edges);
  for (const auto& range : similarities_by_source_) {
    archive(std::distance(similarities_.begin(), range.begin()));
  }
  for (const auto& edge : similarities_) {
    archive(edge.neighbor, edge.similarity);
  }
}

template<class CerealArchive>
void NeighborOrder::load(CerealArchive& archive) {
  size_t num_vertices;
  size_t num_edges;
  archive(num_vertices, num_edges);
  similarities_ = pbbs::sequence<EdgeSimilarity>::no_init(num_edges);
  similarities_by_source_ =
    pbbs::sequence<pbbs::range<EdgeSimilarity*>>::no_init(num_vertices);
  if (num_vertices > 0) {
    size_t prev_offset;
    archive(prev_offset);
    for (size_t i{1}; i < num_vertices; i++) {
      size_t offset;
      archive(offset);
      similarities_by_source_[i - 1] = similarities_.slice(prev_offset, offset);
      prev_offset = offset;
    }
    similarities_by_source_[num_vertices - 1] =
      similarities_.slice(prev_offset, num_edges);
  }
  uintE source{0};
  for (size_t i{0}; i < num_edges; i++) {
    while (static_cast<size_t>(std::distance(
          similarities_.begin(),
          similarities_by_source_[source].end())) <= i) {
      ++source;
    }
    uintE neighbor;
    float similarity;
    archive(neighbor, similarity);
    similarities_[i] = EdgeSimilarity{
      .source = source,
      .neighbor = neighbor,
      .similarity = similarity};
  }
}

template<class CerealArchive>
void CoreOrder::save(CerealArchive& archive) const {
  archive(num_vertices_, order_.size());
  for (const auto& cores : order_) {
    archive(cores.size());
    for (const auto& core : cores) {
      archive(core);
    }
  }
}

template<class CerealArchive>
void CoreOrder::load(CerealArchive& archive) {
  size_t max_mu;
  archive(num_vertices_, max_mu);
  order_ = sequence<pbbs::sequence<CoreThreshold>>(max_mu);
  for (auto& cores : order_) {
    size_t num_cores;
    archive(num_cores);
    cores = pbbs::sequence<CoreThreshold>::no_init(num_cores);
    for (auto& core : cores) {
      archive(core);
    }
  }
}

}  // namespace internal

}  // namespace indexed_scan
