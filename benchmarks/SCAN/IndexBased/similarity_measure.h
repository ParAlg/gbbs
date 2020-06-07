// TODO add comment
#pragma once

#include <atomic>
#include <functional>
#include <limits>

#include "benchmarks/SCAN/IndexBased/intersect.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"
#include "gbbs/bridge.h"
#include "gbbs/graph_mutation.h"
#include "pbbslib/monoid.h"
#include "pbbslib/random.h"
#include "pbbslib/seq.h"
#include "pbbslib/sequence_ops.h"
#include "pbbslib/utilities.h"

namespace scan {

struct EdgeSimilarity {
  // Source vertex ID.
  uintE source;
  // Neighbor vertex ID.
  uintE neighbor;
  // Similarity of source vertex to neighbor vertex.
  float similarity;
};
// Beware that this equality operator compares the floating-point field
// approximately. This is convenient for unit tests but might not be appropriate
// for other uses (e.g., this operator is not transitive).
bool operator==(const EdgeSimilarity&, const EdgeSimilarity&);
std::ostream& operator<<(std::ostream& os, const EdgeSimilarity&);

// TODO add comment
class CosineSimilarity {
 public:
  CosineSimilarity() = default;

  template <template <typename> class VertexTemplate>
  pbbs::sequence<EdgeSimilarity>
  AllEdges(symmetric_graph<VertexTemplate, pbbs::empty>* graph) const;
};

// TODO add comment
struct ApproxCosineSimilarity {
 public:
  ApproxCosineSimilarity(uint32_t num_samples, size_t random_seed);

  template <template <typename> class VertexTemplate>
  pbbs::sequence<EdgeSimilarity>
  AllEdges(symmetric_graph<VertexTemplate, pbbs::empty>* graph) const;

 private:
  const uint32_t num_samples_;
  const size_t random_seed_;
};

//////////////
// Internal //
//////////////

namespace internal {

pbbs::sequence<float>
RandomNormalNumbers(const size_t num_numbers, const pbbs::random rng);

}  // namespace internal

template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> CosineSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph) const {
  // To compute cosine similarities, we need to count shared neighbors between
  // two vertices. Counting the neighbors shared between adjacent vertices u and
  // v is the same as counting the number of triangles that the edge {u, v}
  // appears in.
  //
  // The triangle counting logic here is borrowed from
  // `Triangle_degree_ordering()` in
  // `benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h`. We modify it to
  // maintain triangle counts for each edge.

  // Create a directed version of `graph`, pointing edges from lower degree
  // vertices to higher degree vertices in order to bound the maximum degree of
  // each vertex.
  uintE* vertex_degree_ranking{rankNodes(*graph, graph->n)};
  const auto filter_predicate{[&](const uintE u, const uintE v, pbbs::empty) {
    return vertex_degree_ranking[u] < vertex_degree_ranking[v];
  }};
  auto directed_graph{graph->filterGraph(*graph, filter_predicate)};

  // Estimate the amount of work for each vertex for better load balancing.
  auto parallel_work{sequence<size_t>(directed_graph.n)};
  {
    const auto get_edge_work{[&](const uintE u, const uintE v, pbbs::empty) {
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
          pbbs::empty,
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
          internal::intersect_f_with_index_par(
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

  pbbs::sequence<EdgeSimilarity> similarities(graph->m);
  // Convert shared neighbor counts into structural similarities for each edge.
  par_for(0, directed_graph.n, [&](const size_t vertex_id) {
    const uintT v_counter_offset{counter_offsets[vertex_id]};
    const float v_neighborhood_sqrt{
      sqrtf(graph->get_vertex(vertex_id).getOutDegree() + 1)};
    const auto compute_similarity{[&](
        const uintE v_id,
        const uintE u_id,
        pbbs::empty,
        const uintE v_to_u_index) {
      // SCAN structural similarities are defined using _closed_ neighborhoods,
      // hence the need to to adjust these values by `+ 1` and `+ 2`.
      const float u_neighborhood_sqrt{
        sqrtf(graph->get_vertex(u_id).getOutDegree() + 1)};
      const uintT counter_index{v_counter_offset + v_to_u_index};
      const uintE num_shared_neighbors{counters[counter_index] + 2};
      const float structural_similarity{
        num_shared_neighbors / (v_neighborhood_sqrt * u_neighborhood_sqrt)};
      similarities[2 * counter_index] =
        {.source = v_id, .neighbor = u_id, .similarity = structural_similarity};
      similarities[2 * counter_index + 1] =
        {.source = u_id, .neighbor = v_id, .similarity = structural_similarity};
    }};
    directed_graph.get_vertex(vertex_id).mapOutNghWithIndex(
        vertex_id, compute_similarity);
  });

  directed_graph.del();
  pbbs::free_array(vertex_degree_ranking);
  return similarities;
}

template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> ApproxCosineSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph) const {
  using Vertex = VertexTemplate<pbbs::empty>;
  // Approximates cosine similarity using SimHash (c.f. "Similarity Estimation
  // Techniques from Rounding Algorithms" by Moses Charikar).

  // TODO more comments explaining what's happening here

  const size_t num_vertices{graph->n};
  const pbbs::sequence<float> normals{internal::RandomNormalNumbers(
      num_vertices * num_samples_, pbbs::random{random_seed_})};
  const auto addition_monoid{pbbs::addm<float>{}};
  // TODO performance optimization: chunking bools into fixed size bitsets or
  // `uint64_t`s
  const pbbs::sequence<pbbs::sequence<bool>> vertex_fingerprints{
    num_vertices,
    [&](const size_t vertex_id) {
      Vertex vertex{graph->get_vertex(vertex_id)};
      return pbbs::sequence<bool>{
        num_samples_,
        [&](const size_t sample_id) {
          const size_t offset{num_vertices * sample_id};
          const auto neighbor_to_normal{
            [&](uintE, const uintE neighbor, pbbs::empty) {
              return normals[offset + neighbor];
            }};
          return (normals[offset + vertex_id] +
              vertex.template reduceOutNgh<float>(
              vertex_id, neighbor_to_normal, addition_monoid)) >= 0;
        }};
    }};

  // TODO(tomtseng): this work is duplicated in
  // indexed_scan::internal::NeighborOrder::NeighborOrder
  pbbs::sequence<uintT> vertex_offsets{
      graph->n,
      [&](const size_t i) { return graph->get_vertex(i).getOutDegree(); }};
  pbbslib::scan_add_inplace(vertex_offsets);

  pbbs::sequence<EdgeSimilarity> undirected_similarities{
    pbbs::sequence<EdgeSimilarity>::no_init(graph->m)};
  par_for(0, undirected_similarities.size(), [&](const size_t i) {
    undirected_similarities[i].similarity =
      std::numeric_limits<float>::quiet_NaN();
  });
  // Get similarities for all edges (u, v) where u < v.
  const auto compute_similarity{[&](
      const uintE vertex_id,
      const uintE neighbor_id,
      pbbs::empty,
      const uintE neighbor_index) {
    if (vertex_id < neighbor_id) {
      const pbbs::sequence<bool>& vertex_fingerprint{
        vertex_fingerprints[vertex_id]};
      const pbbs::sequence<bool>& neighbor_fingerprint{
        vertex_fingerprints[neighbor_id]};
      const auto fingerprint_xor{
        pbbs::delayed_seq<std::remove_const<decltype(num_samples_)>::type>(
          num_samples_,
          [&](const size_t i) {
            return vertex_fingerprint[i] != neighbor_fingerprint[i];
          })};
      const float angle_estimate{static_cast<float>(
          pbbslib::reduce_add(fingerprint_xor) * M_PI / num_samples_)};
      const float cosine_similarity_estimate{std::cos(angle_estimate)};
      undirected_similarities[vertex_offsets[vertex_id] + neighbor_index] = {
        .source = vertex_id,
        .neighbor = neighbor_id,
        .similarity = cosine_similarity_estimate};
    }
  }};
  par_for(0, num_vertices, [&](const size_t vertex_id) {
    graph->get_vertex(vertex_id).mapOutNghWithIndex(
        vertex_id, compute_similarity);
  });
  // Copy similarities for edges (u, v) where u > v.
  const size_t half_num_edges{graph->m / 2};
  pbbs::sequence<EdgeSimilarity> similarities{
    pbbs::sequence<EdgeSimilarity>::no_init(graph->m)};
  const auto is_valid_similarity{
    [](const EdgeSimilarity& edge) { return !std::isnan(edge.similarity); }};
  pbbs::filter_out(
      undirected_similarities, similarities.slice(), is_valid_similarity);
  par_for(0, half_num_edges, [&](const size_t i) {
      const EdgeSimilarity& edge{similarities[i]};
      similarities[i + half_num_edges] = {
        .source = edge.neighbor,
        .neighbor = edge.source,
        .similarity = edge.similarity};
  });
  return similarities;
}

}  // namespace scan
