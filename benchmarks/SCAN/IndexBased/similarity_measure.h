// Similarity measures for determining the similarity of adjacent pairs of
// vertices.
//
// TODO(tomtseng) needs refactoring -- this code has so much copy-and-paste code
// and is confusing
#pragma once

#include <array>
#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <type_traits>

#include "benchmarks/SCAN/IndexBased/intersect.h"
#include "benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h"
#include "gbbs/bridge.h"
#include "gbbs/graph_mutation.h"
#include "pbbslib/monoid.h"
#include "pbbslib/random.h"
#include "pbbslib/seq.h"
#include "pbbslib/sequence_ops.h"
#include "pbbslib/utilities.h"

namespace gbbs {

namespace scan {

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

/////////////////////////
// Similarity measures //
/////////////////////////
// The similarity measure classes implement the following functions:
//   // Returns a `graph->m`-length sequence containing the similarity score
//   // between every adjacent pair of vertices in the graph. The neighbor lists
//   // for each vertex of the graph must be sorted by ascending neighbor ID.
//   template <class Graph>
//   pbbs::sequence<EdgeSimilarity> AllEdges(Graph* graph) const;


// The cosine similarity between two adjacent vertices u and v is
//   (size of intersection of the closed neighborhoods of u and v) /
//   (geometric mean of size of the closed neighborhoods of u and of v)
// where the closed neighborhood of a vertex x consists of all neighbors of x
// along with x itself.
//
// How is this a cosine similarity? We can treat the neighborhood of a vertex
// v as an n-dimensional vector where the i-th entry of the vector is 1 if
// vertex i is in vertex v's neighborhood and is 0 otherwise. Then the cosine
// similarity between the vectors corresponding to the neighborhoods of two
// vertices u and v is the equation above.
//
// Cosine similarity also works for weighted graphs, in which case the
// definition is
//   (sum(weight({u, z}) * weight({v, z}) for z in the intersection of the
//      closed neighborhoods of u and v) / (norm(u) * norm(v)))
// where for any v,
//   norm(v) = sqrt(sum(weight({v, z})^2 for z in the closed neighborhood of v))
// and
//   weight({v, v}) = 1.
class CosineSimilarity {
 public:
  CosineSimilarity() = default;

  template <template <typename> class VertexTemplate, typename Weight>
  pbbs::sequence<EdgeSimilarity>
  AllEdges(symmetric_graph<VertexTemplate, Weight>* graph) const;
};

// The Jaccard similarity between two adjacent vertices u and v is
//   (size of intersection of the closed neighborhoods of u and v) /
//   (size of union of the closed neighborhoods of u and of v)
// where the closed neighborhood of a vertex x consists of all neighbors of x
// along with x itself.
class JaccardSimilarity {
 public:
  JaccardSimilarity() = default;

  template <template <typename> class VertexTemplate>
  pbbs::sequence<EdgeSimilarity>
  AllEdges(symmetric_graph<VertexTemplate, pbbs::empty>* graph) const;
};

// This is an approximate version of `CosineSimilarity`. Increasing
// `num_samples` increases the approximation accuracy.
//
// Let `m` be the number of undirected edges in the graph, and let `a` and `b`
// be in the range (0, 1). Then, if we replace the random number generator used
// within the code with perfectly random number generator, then picking
//   num_samples = 1.5 * pi ^ 2 * ln(2 * m / a) / b ^ 2
// gives that with probability at least `1 - a`, each edge receives the correct
// cosine similarity with absolute error up to `b`. In practice, setting
// num_samples so high is probably excessive.
//
// This is a biased estimate of the cosine similarity.
//
// This is really only helpful for graphs with lots of high degree vertices.
// Otherwise, the cost to approximate similarities with enough samples to have
// good accuracy outweighs the cost to compute similarities exactly.
struct ApproxCosineSimilarity {
 public:
  ApproxCosineSimilarity(uint32_t num_samples, size_t random_seed);

  // When `random_seed` is fixed, the output of `AllEdges` is deterministic.
  template <template <typename> class VertexTemplate, typename Weight>
  pbbs::sequence<EdgeSimilarity>
  AllEdges(symmetric_graph<VertexTemplate, Weight>* graph) const;

 private:
  const uint32_t num_samples_;
  const size_t random_seed_;
};

// This is an approximate version of `JaccardSimilarity`. Increasing
// `num_samples` increases the approximation accuracy.
//
// This is really only helpful for graphs with lots of high degree vertices.
// Otherwise, the cost to approximate similarities with enough samples to have
// good accuracy outweighs the cost to compute similarities exactly.
struct ApproxJaccardSimilarity {
 public:
  ApproxJaccardSimilarity(uint32_t num_samples, size_t random_seed);

  // When `random_seed` is fixed, the output of `AllEdges` is deterministic.
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

// Type for counting the number or weight of shared neighbors between two
// vertices.
template <typename Weight>
struct EdgeCounter {
  using type = int64_t;  // weighted case
};

template <>
struct EdgeCounter<pbbs::empty> {
  using type = uintE;  // unweighted case
};

// Compute (numerator / denominator), rounding up if there's any remainder.
// Numerator must be positive.
constexpr uint64_t
DivideRoundingUp(const size_t numerator, const size_t denominator) {
  return (numerator - 1) / denominator + 1;
}

// Pseudorandomly generate `num_numbers` random normal numbers, each with zero
// mean and unit variance.
pbbs::sequence<float> RandomNormalNumbers(size_t num_numbers, pbbs::random rng);

// Create a directed version of `graph`, pointing edges from lower degree
// vertices to higher degree vertices. This upper bounds the out-degree of each
// vertex in the directed graph with `sqrt(graph->m)`.
template<template <typename> class VertexTemplate, typename Weight>
auto DirectGraphByDegree(symmetric_graph<VertexTemplate, Weight>* graph) {
  uintE* vertex_degree_ranking{rankNodes(*graph, graph->n)};
  const auto filter_predicate{[&](const uintE u, const uintE v, Weight) {
    return vertex_degree_ranking[u] < vertex_degree_ranking[v];
  }};
  auto directed_graph{filterGraph(*graph, filter_predicate)};
  pbbs::free_array(vertex_degree_ranking);
  return directed_graph;
}

// Returns a sequence `vertex_offsets` such that if there is another sequence
// `edges` consisting of the out-edges of `*graph` sorted by source vertex, then
// `vertex_offsets[i]` is the first appearance of vertex i as a source vertex.
template <template <typename> class VertexTemplate, typename Weight>
pbbs::sequence<uintT>
VertexOutOffsets(symmetric_graph<VertexTemplate, Weight>* graph) {
  pbbs::sequence<uintT> vertex_offsets{
      graph->n,
      [&](const size_t i) { return graph->get_vertex(i).getOutDegree(); }};
  pbbslib::scan_add_inplace(vertex_offsets);
  return vertex_offsets;
}

// Returns a `graph->m`-length sequence containing the similarity score
// between every adjacent pair of vertices u and v. The similarity score is
// provided by `neighborhood_sizes_to_similarity` and must be a function of the
// sizes of the neighborhoods of u and v and the size of the intersection of the
// neighborhoods.
//
// Arguments:
//   graph
//     Graph on which to compute similarities.
//   neighborhood_sizes_to_similarity: (uintE, uintE, uintE) -> float
//     Function that computes the similarity between adjacent vertices u and v,
//     taking (size of u's neighborhood, size of v's neighborhood, size of the
//     the intersection of the two neighborhoods) as arguments. The function
//     should be symmetric, i.e., give the same output when u and v are swapped.
template <template <typename> class VertexTemplate, class F>
pbbs::sequence<EdgeSimilarity> AllEdgeNeighborhoodSimilarities(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph,
    F&& neighborhood_sizes_to_similarity) {
  // Counting the neighbors shared between adjacent vertices u and v is the same
  // as counting the number of triangles that the edge {u, v} appears in.
  //
  // The triangle counting logic here is borrowed from
  // `Triangle_degree_ordering()` in
  // `benchmarks/TriangleCounting/ShunTangwongsan15/Triangle.h`. We modify it to
  // maintain triangle counts for each edge.

  auto directed_graph{DirectGraphByDegree(graph)};
  // Each counter in `counters` holds the number of shared neighbors in `graph`
  // between u and v for some edge {u, v}.
  pbbs::sequence<std::atomic<uintE>> counters(
      directed_graph.m, [](size_t) { return std::atomic<uintE>{0}; });
  // We use `counter_offsets` to index into `counters` for each edge.
  const pbbs::sequence<uintT> counter_offsets{
    internal::VertexOutOffsets(&directed_graph)};
  // Find triangles of the following form:
  //        w
  //       ^ ^
  //      /   \.
  //     u --> v
  // There's a bijection between triangles of this form in `directed_graph` and
  // undirected triangles in `graph`.
  par_for(0, graph->n, [&](const size_t vertex_id) {
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
            &vertex, &neighbor, update_counters);
    }};
    constexpr bool kParallel{false};
    vertex.mapOutNghWithIndex(vertex_id, intersect, kParallel);
  });

  pbbs::sequence<EdgeSimilarity> similarities(graph->m);
  // Convert shared neighbor counts into similarities for each edge.
  par_for(0, directed_graph.n, [&](const size_t vertex_id) {
    const uintT v_counter_offset{counter_offsets[vertex_id]};
    const uintE v_degree{graph->get_vertex(vertex_id).getOutDegree()};
    const auto compute_similarity{[&](
        const uintE v_id,
        const uintE u_id,
        pbbs::empty,
        const uintE v_to_u_index) {
      const uintT counter_index{v_counter_offset + v_to_u_index};
      const uintE num_shared_neighbors{counters[counter_index]};
      const uintE u_degree{graph->get_vertex(u_id).getOutDegree()};
      const float similarity{neighborhood_sizes_to_similarity(
          v_degree, u_degree, num_shared_neighbors)};
      similarities[2 * counter_index] =
        {.source = v_id, .neighbor = u_id, .similarity = similarity};
      similarities[2 * counter_index + 1] =
        {.source = u_id, .neighbor = v_id, .similarity = similarity};
    }};
    directed_graph.get_vertex(vertex_id).mapOutNghWithIndex(
        vertex_id, compute_similarity);
  });

  directed_graph.del();
  return similarities;
}

// Implementation of ApproxCosineSimilarities::AllEdges.
//
// `degree_threshold` is a threshold so that we only approximate the similarity
// score between two vertices if their degrees are high enough. (When the
// degrees are low, it's cheap to compute the similarity exactly.)
template <template <typename> class VertexTemplate, typename Weight>
pbbs::sequence<EdgeSimilarity> ApproxCosineEdgeSimilarities(
    symmetric_graph<VertexTemplate, Weight>* graph,
    const uint32_t num_samples,
    const size_t degree_threshold,
    const size_t random_seed) {
  // Approximates cosine similarity using SimHash (c.f. "Similarity Estimation
  // Techniques from Rounding Algorithms" by Moses Charikar).
  //
  // The idea is that we can estimate the angle between two n-dimensional
  // vectors by drawing a random n-dimensional hyperplane and determining which
  // side of the hyperplane the vectors fall on. The larger the angle between
  // the two vectors, the more likely that the two vectors will fall on
  // opposite sides of the hyperplane. Repeat this for several random
  // hyperplanes.
  // Represent a hyperplane by a vector orthogonal to that hyperplane. Generate
  // that uniformly random orthogonal vector by drawing i.i.d. normal variables
  // for each dimension. Determine which side of the hyperplane vectors fall on
  // by taking the dot product with the orthogonal vector.
  //
  // For edges between high degree vertices, estimate the similarity with
  // SimHash. For edges with a low degree vertex, compute the similarity exactly
  // with triangle counting like in `AllEdgeNeighborhoodSimilarities()`.
  using Vertex = VertexTemplate<Weight>;
  // We compute `num_samples_` hyperplanes and, to sketch a vertex's
  // neighborhood vector, we compute `num_samples_` bits representing the sign
  // of the vector's dot product with each hyperplane. For efficiency, we store
  // the bits in chunks of `kBitArraySize` rather than one-by-one.
  using BitArray = uint64_t;
  constexpr size_t kBitArraySize{sizeof(BitArray) * 8};

  // Computing random normal numbers is expensive, so we precompute which
  // vertices need assignments of normal numbers for Minhash fingerprinting.
  pbbs::sequence<bool> needs_fingerprint_seq(graph->n, false);
  pbbs::sequence<uintE> needs_normals_seq(graph->n, 0U);
  par_for(0, graph->n, [&](const size_t vertex_id) {
    Vertex vertex{graph->get_vertex(vertex_id)};
    if (vertex.getOutDegree() >= degree_threshold) {
      // Vertex should be fingerprinted if both it and one of its neighbors
      // has high degree. If a vertex needs to be fingerprinted, then normal
      // random numbers should be generated for it and its neighbors.
      bool skip_fingerprint{true};
      const auto check_degree_threshold{
        [&](uintE, const uintE neighbor_id, Weight) {
          if (skip_fingerprint &&
              graph->v_data[neighbor_id].degree >= degree_threshold) {
            skip_fingerprint = false;
          }
        }};
      vertex.mapOutNgh(vertex_id, check_degree_threshold);
      if (!skip_fingerprint) {
        needs_fingerprint_seq[vertex_id] = true;
        needs_normals_seq[vertex_id] = true;
        const auto set_needs_normals{
          [&](uintE, const uintE neighbor_id, Weight) {
            needs_normals_seq[neighbor_id] = true;
          }};
        vertex.mapOutNgh(vertex_id, set_needs_normals);
      }
    }
  });
  // repurpose `needs_normals_seq` to serve as the index of a vertex into
  // `normals`
  const uintE num_needs_normals{pbbslib::scan_add_inplace(needs_normals_seq)};
  const pbbs::sequence<float> normals{RandomNormalNumbers(
      num_needs_normals * num_samples, pbbs::random{random_seed})};

  const size_t num_vertices{graph->n};
  const size_t num_bit_arrays{
    internal::DivideRoundingUp(num_samples, kBitArraySize)};
  // Simhash fingerprints.
  const pbbs::sequence<pbbs::sequence<BitArray>> vertex_fingerprints{
    num_vertices,
    [&](const size_t vertex_id) {
      if (!needs_fingerprint_seq[vertex_id]) {
        return pbbs::sequence<BitArray>{};
      }
      Vertex vertex{graph->get_vertex(vertex_id)};
      const uintE vertex_normal_offset{
        num_samples * needs_normals_seq[vertex_id]};
      return pbbs::sequence<BitArray>{
        num_bit_arrays,
        [&](const size_t bit_array_id) {
          BitArray bits{0};
          const size_t max_bit_id{
            bit_array_id + 1 == num_bit_arrays
            ? kBitArraySize - (num_bit_arrays * kBitArraySize - num_samples)
            : kBitArraySize};
          const size_t bits_offset{bit_array_id * kBitArraySize};
          std::array<float, kBitArraySize> hyperplane_dot_products;
          for (size_t bit_id = 0; bit_id < max_bit_id; bit_id++) {
            hyperplane_dot_products[bit_id] =
              normals[vertex_normal_offset + bits_offset + bit_id];
          }
          const auto update_dot_products{[&](
              uintE,
              const uintE neighbor_id,
              [[maybe_unused]] const Weight weight) {
            for (size_t bit_id = 0; bit_id < max_bit_id; bit_id++) {
              if constexpr (std::is_same<Weight, pbbslib::empty>::value) {
                // unweighted case
                hyperplane_dot_products[bit_id] +=
                  normals[num_samples * needs_normals_seq[neighbor_id]
                    + bits_offset + bit_id];
              } else {  // weighted case
                hyperplane_dot_products[bit_id] +=
                  normals[num_samples * needs_normals_seq[neighbor_id]
                    + bits_offset + bit_id] * weight;
              }
            }
          }};
          constexpr bool kParallel{false};
          vertex.mapOutNgh(vertex_id, update_dot_products, kParallel);
          for (size_t bit_id = 0; bit_id < max_bit_id; bit_id++) {
            if (hyperplane_dot_products[bit_id] >= 0) {
              bits |= (1LL << bit_id);
            }
          }
          return bits;
        }};
    }};

  auto directed_graph{DirectGraphByDegree(graph)};
  using Counter = typename EdgeCounter<Weight>::type;
  // For unweighted graphs, each counter in `counters` holds the number of
  // shared neighbors in `graph` between u and v for some edge {u, v} -- the
  // numerator in CosineSimilarity(u, v).
  //
  // For weighted graphs, each entry still corresponds to the numerator in
  // CosineSimilarity(u, v). Even if the graph is weighted with floating-point
  // numbers, we still make the counter an integer type since
  // std::atomic<double> doesn't have an atomic increment operation in C++17,
  // and we instead round the edge weights to integers.
  pbbs::sequence<std::atomic<Counter>> counters(
      directed_graph.m, [](size_t) { return std::atomic<Counter>{0}; });
  // For floating-point-weighted graphs, we multiply all edge weights by this
  // factor and round to the nearest integer. Multiplying all edge weights by a
  // constant doesn't affect cosine similarities. The multiplication preserves
  // more accuracy after rounding.
  constexpr int64_t kWeightFactor =
    std::is_floating_point<Weight>::value ? 1000 : 1;
  // We use `counter_offsets` to index into `counters` for each edge.
  const pbbs::sequence<uintT> counter_offsets{
    internal::VertexOutOffsets(&directed_graph)};
  // Find triangles of the following form:
  //        w
  //       ^ ^
  //      /   \.
  //     u --> v
  // Count each of these triangles to get the number of shared neighbors between
  // vertices. However, we skip pairs of vertices that have high degree in the
  // original, undirected graph.
  par_for(0, graph->n, [&](const size_t vertex_id) {
    auto vertex{directed_graph.get_vertex(vertex_id)};
    const bool vertex_is_high_degree{
      graph->v_data[vertex_id].degree >= degree_threshold};
    if (vertex_is_high_degree) {
      // Since all edges in the directed graph point towards higher degree
      // vertices, if the current vertex is high degree, so are all its directed
      // neighbors. Skip these pairs of edges since we'll approximate their
      // similarities.
      return;
    }

    const uintT vertex_counter_offset{counter_offsets[vertex_id]};
    const auto intersect{[&](
        const uintE v_id,
        const uintE neighbor_id,
        [[maybe_unused]] const Weight weight,
        const uintE v_to_neighbor_index) {
      auto neighbor{directed_graph.get_vertex(neighbor_id)};
      const bool neighbor_is_high_degree{
        graph->v_data[neighbor_id].degree >= degree_threshold};
      const uintT neighbor_counter_offset{counter_offsets[neighbor_id]};

      if constexpr (std::is_same<Weight, pbbslib::empty>::value) {
        // unweighted case
        const auto update_counters{[&](
            const uintE shared_neighbor,
            const uintE vertex_to_shared_index,
            const uintE neighbor_to_shared_index) {
          counters[vertex_counter_offset + vertex_to_shared_index]++;
          if (!(neighbor_is_high_degree &&
              graph->v_data[shared_neighbor].degree >= degree_threshold)) {
            counters[neighbor_counter_offset + neighbor_to_shared_index]++;
          }
        }};
        counters[vertex_counter_offset + v_to_neighbor_index] +=
          internal::intersect_f_with_index_par(
              &vertex, &neighbor, update_counters);
      } else if constexpr (std::is_floating_point<Weight>::value) {
        // weighted case
        const auto update_counters{[&](
          const uintE shared_neighbor,
          const uintE vertex_to_shared_index,
          const uintE neighbor_to_shared_index,
          const Weight weight_1,
          const Weight weight_2) {
        counters[vertex_counter_offset + vertex_to_shared_index] +=
          kWeightFactor * kWeightFactor * weight * weight_2;
        if (!(neighbor_is_high_degree &&
            graph->v_data[shared_neighbor].degree >= degree_threshold)) {
          counters[neighbor_counter_offset + neighbor_to_shared_index] +=
            kWeightFactor * kWeightFactor * weight * weight_1;
        }
      }};
      counters[vertex_counter_offset + v_to_neighbor_index] +=
        kWeightFactor * kWeightFactor *
        internal::intersect_f_with_index_par(
            &vertex, &neighbor, update_counters);
      }
    }};
    constexpr bool kParallel{false};
    vertex.mapOutNghWithIndex(vertex_id, intersect, kParallel);
  });

  const pbbs::sequence<double> norms{[&]() -> pbbs::sequence<double> {
    if constexpr (std::is_same<Weight, pbbslib::empty>::value) {  // unweighted
      return {};
    } else {  // weighted case
      return pbbs::sequence<double>{
        graph->n,
        [&](const size_t vertex_id) {
          double norm = kWeightFactor * kWeightFactor;  // add self-loop weight
          auto vertex = graph->get_vertex(vertex_id);
          constexpr bool kParallel{false};
          const auto update_norm{
            [&](uintE, uintE, const Weight weight) {
              norm += kWeightFactor * kWeightFactor * weight * weight;
            }};
          vertex.mapOutNgh(vertex_id, update_norm, kParallel);
          return std::sqrt(norm);
        }};
    }
  }()};

  pbbs::sequence<EdgeSimilarity> similarities(graph->m);
  // Convert shared neighbor counts into similarities for each edge.
  par_for(0, directed_graph.n, [&](const size_t vertex_id) {
    const uintT v_counter_offset{counter_offsets[vertex_id]};
    const uintE v_degree{graph->get_vertex(vertex_id).getOutDegree()};
    const bool vertex_is_high_degree{v_degree >= degree_threshold};
    const pbbs::sequence<uint64_t>& vertex_fingerprint{
      vertex_fingerprints[vertex_id]};
    const auto compute_similarity{[&](
        const uintE v_id,
        const uintE u_id,
        [[maybe_unused]] const Weight weight,
        const uintE v_to_u_index) {
      const uintT counter_index{v_counter_offset + v_to_u_index};
      float similarity{-1};
      if (vertex_is_high_degree) {  // approximate similarity
        const pbbs::sequence<BitArray>& neighbor_fingerprint{
          vertex_fingerprints[u_id]};
        const auto fingerprint_xor{
          pbbs::delayed_seq<std::remove_const<decltype(num_samples)>::type>(
            vertex_fingerprint.size(),
            [&](const size_t i) {
              return __builtin_popcountll(
                  vertex_fingerprint[i] ^ neighbor_fingerprint[i]);
            })};
        const float angle_estimate{static_cast<float>(
            pbbslib::reduce_add(fingerprint_xor) * M_PI / num_samples)};
        similarity = std::cos(angle_estimate);
      } else {  // exact similarity
        if constexpr (std::is_same<Weight, pbbslib::empty>::value) {
          // unweighted case
          const uintE num_shared_neighbors{counters[counter_index]};
          const uintE u_degree{graph->get_vertex(u_id).getOutDegree()};
          similarity = (num_shared_neighbors + 2) /
            (sqrtf(v_degree + 1) * sqrtf(u_degree + 1));
        } else {  // weighted case
          // additional term is to account for self-loop edge in closed
          // neighborhood
          const double shared_weight{
            counters[counter_index]
              + 2 * kWeightFactor * kWeightFactor * weight};
          similarity = static_cast<float>(
              shared_weight / (norms[v_id] * norms[u_id]));
        }
      }
      similarities[2 * counter_index] =
        {.source = v_id, .neighbor = u_id, .similarity = similarity};
      similarities[2 * counter_index + 1] =
        {.source = u_id, .neighbor = v_id, .similarity = similarity};
    }};
    directed_graph.get_vertex(vertex_id).mapOutNghWithIndex(
        vertex_id, compute_similarity);
  });

  directed_graph.del();
  return similarities;
}

// Implementation of ApproxJaccardSimilarities::AllEdges.
//
// `degree_threshold` is a threshold so that we only approximate the similarity
// score between two vertices if their degrees are high enough. (When the
// degrees are low, it's cheap to compute the similarity exactly.)
template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> ApproxJaccardEdgeSimilarities(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph,
    const uint32_t original_num_samples,
    const size_t degree_threshold,
    const size_t random_seed) {
  using Weight = pbbs::empty;
  using Vertex = VertexTemplate<Weight>;
  // For edges between high degree vertices, estimate the Jaccard similarity
  // with a MinHash variant --- see paper "One Permutation Hashing for Efficient
  // Search and Learning." For edges with a low degree vertex, compute the
  // Jaccard similarity exactly with triangle counting like in
  // `AllEdgeNeighborhoodSimilarities()`.
  const uint32_t log_num_samples{
    std::max<uint32_t>(pbbs::log2_up(original_num_samples), 1)};
  const uintE num_samples{static_cast<uintE>(1ULL << log_num_samples)};
  const uintE bucket_mask{num_samples - 1};
  constexpr uintE kEmptyBucket{UINT_E_MAX};
  const size_t num_vertices{graph->n};
  const pbbs::sequence<uintE> vertex_permutation{
    pbbs::random_permutation<uintE>(num_vertices, pbbs::random{random_seed})};

  // Compute MinHash fingerprints for high degree vertices.
  const pbbs::sequence<pbbs::sequence<uintE>> vertex_fingerprints{
    num_vertices,
    [&](const size_t vertex_id) {
      Vertex vertex{graph->get_vertex(vertex_id)};
      if (vertex.getOutDegree() < degree_threshold) {
        return pbbs::sequence<uintE>{};
      }
      bool skip_fingerprint{true};
      const auto check_degree_threshold{
        [&](uintE, const uintE neighbor_id, Weight) {
          if (skip_fingerprint &&
              graph->v_data[neighbor_id].degree >= degree_threshold) {
            skip_fingerprint = false;
          }
        }};
      vertex.mapOutNgh(vertex_id, check_degree_threshold);
      if (skip_fingerprint) {
        return pbbs::sequence<uintE>{};
      }

      pbbs::sequence<uintE> fingerprint(num_samples, kEmptyBucket);
      const auto update_fingerprint{
        [&](uintE, const uintE neighbor, pbbs::empty) {
          const uintE permuted_neighbor{vertex_permutation[neighbor]};
          const uintE bucket_id{permuted_neighbor & bucket_mask};
          const uintE bucket_value{permuted_neighbor >> log_num_samples};
          pbbs::write_min(
              &(fingerprint[bucket_id]), bucket_value, std::less<uintE>{});
        }};
      update_fingerprint(vertex_id, vertex_id, pbbs::empty{});
      vertex.mapOutNgh(vertex_id, update_fingerprint);
      return fingerprint;
    }};

  auto directed_graph{DirectGraphByDegree(graph)};
  // Each counter in `counters` holds the number of shared neighbors in `graph`
  // between u and v for some edge {u, v}.
  pbbs::sequence<std::atomic<uintE>> counters(
      directed_graph.m, [](size_t) { return std::atomic<uintE>{0}; });
  // We use `counter_offsets` to index into `counters` for each edge.
  const pbbs::sequence<uintT> counter_offsets{
    internal::VertexOutOffsets(&directed_graph)};
  // Find triangles of the following form:
  //        w
  //       ^ ^
  //      /   \.
  //     u --> v
  // Count each of these triangles to get the number of shared neighbors between
  // vertices. However, we skip pairs of vertices that have high degree in the
  // original, undirected graph.
  par_for(0, graph->n, [&](const size_t vertex_id) {
    auto vertex{directed_graph.get_vertex(vertex_id)};
    const bool vertex_is_high_degree{
      graph->v_data[vertex_id].degree >= degree_threshold};
    if (vertex_is_high_degree) {
      // Since all edges in the directed graph point towards higher degree
      // vertices, if the current vertex is high degree, so are all its directed
      // neighbors. Skip these pairs of edges since we'll approximate their
      // similarities.
      return;
    }

    const uintT vertex_counter_offset{counter_offsets[vertex_id]};
    const auto intersect{[&](
        const uintE v_id,
        const uintE neighbor_id,
        pbbs::empty,
        const uintE v_to_neighbor_index) {
      auto neighbor{directed_graph.get_vertex(neighbor_id)};
      const bool neighbor_is_high_degree{
        graph->v_data[neighbor_id].degree >= degree_threshold};
      const uintT neighbor_counter_offset{counter_offsets[neighbor_id]};
      const auto update_counters{[&](
          const uintE shared_neighbor,
          const uintE vertex_to_shared_index,
          const uintE neighbor_to_shared_index) {
        counters[vertex_counter_offset + vertex_to_shared_index]++;
        if (!(neighbor_is_high_degree &&
            graph->v_data[shared_neighbor].degree >= degree_threshold)) {
          counters[neighbor_counter_offset + neighbor_to_shared_index]++;
        }
      }};
      counters[vertex_counter_offset + v_to_neighbor_index] +=
        internal::intersect_f_with_index_par(
            &vertex, &neighbor, update_counters);
    }};
    constexpr bool kParallel{false};
    vertex.mapOutNghWithIndex(vertex_id, intersect, kParallel);
  });

  pbbs::sequence<EdgeSimilarity> similarities(graph->m);
  // Convert shared neighbor counts into similarities for each edge.
  par_for(0, directed_graph.n, [&](const size_t vertex_id) {
    const uintT v_counter_offset{counter_offsets[vertex_id]};
    const uintE v_degree{graph->get_vertex(vertex_id).getOutDegree()};
    const bool vertex_is_high_degree{v_degree >= degree_threshold};
    const pbbs::sequence<uintE>& vertex_fingerprint{
      vertex_fingerprints[vertex_id]};
    const auto compute_similarity{[&](
        const uintE v_id,
        const uintE u_id,
        pbbs::empty,
        const uintE v_to_u_index) {
      const uintT counter_index{v_counter_offset + v_to_u_index};
      float similarity{-1};
      if (vertex_is_high_degree) {  // approximate similarity
        const pbbs::sequence<uintE>& neighbor_fingerprint{
          vertex_fingerprints[u_id]};
        const uintE fingerprint_matches{
          pbbslib::reduce_add(
            pbbs::delayed_seq<uintE>(
              vertex_fingerprint.size(),
              [&](const size_t i) {
                return vertex_fingerprint[i] != kEmptyBucket &&
                  vertex_fingerprint[i] == neighbor_fingerprint[i];
              }))};
        const uintE fingerprint_empty_count{
          pbbslib::reduce_add(
            pbbs::delayed_seq<uintE>(
              vertex_fingerprint.size(),
              [&](const size_t i) {
                return vertex_fingerprint[i] == kEmptyBucket &&
                  neighbor_fingerprint[i] == kEmptyBucket;
              }))};
        similarity = fingerprint_matches /
          (static_cast<float>(num_samples - fingerprint_empty_count));
      } else {  // exact similarity
        const uintE num_shared_neighbors{counters[counter_index]};
        const uintE u_degree{graph->get_vertex(u_id).getOutDegree()};
        const uintE neighborhood_union{
          v_degree + u_degree - num_shared_neighbors};
        similarity =
          static_cast<float>((num_shared_neighbors + 2)) / neighborhood_union;
      }
      similarities[2 * counter_index] =
        {.source = v_id, .neighbor = u_id, .similarity = similarity};
      similarities[2 * counter_index + 1] =
        {.source = u_id, .neighbor = v_id, .similarity = similarity};
    }};
    directed_graph.get_vertex(vertex_id).mapOutNghWithIndex(
        vertex_id, compute_similarity);
  });

  directed_graph.del();
  return similarities;
}

}  // namespace internal

template <template <typename> class VertexTemplate, typename Weight>
pbbs::sequence<EdgeSimilarity> CosineSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, Weight>* graph) const {
  if constexpr (std::is_same<Weight, pbbslib::empty>::value) {  // unweighted
    constexpr auto similarity_func{
      [](const uintE neighborhood_size_1,
         const uintE neighborhood_size_2,
         const uintE num_shared_neighbors) {
      // SCAN structural/cosine similarities are defined using _closed_
      // neighborhoods, hence the need to to adjust these values by `+ 1` and
      // `+ 2`.
      return (num_shared_neighbors + 2) /
        (sqrtf(neighborhood_size_1 + 1) *
         sqrtf(neighborhood_size_2 + 1));
    }};
    return internal::AllEdgeNeighborhoodSimilarities(graph, similarity_func);
  } else {  // weighted case
    auto directed_graph{internal::DirectGraphByDegree(graph)};
    // Each counter in `counters` will hold numerator of CosineSimilarity(u, v)
    // for some edge {u, v}.
    pbbs::sequence<std::atomic<int64_t>> counters(
        directed_graph.m, [](size_t) { return std::atomic<int64_t>{0}; });
    // For floating-point-weighted graphs, we multiply all edge weights by this
    // factor and round to the nearest integer. Multiplying all edge weights by
    // a constant doesn't affect cosine similarities. The multiplication
    // preserves more accuracy after rounding.
    constexpr int64_t kWeightFactor =
      std::is_floating_point<Weight>::value ? 1000 : 1;
    // We use `counter_offsets` to index into `counters` for each edge.
    const pbbs::sequence<uintT> counter_offsets{
      internal::VertexOutOffsets(&directed_graph)};
    // Find triangles of the following form:
    //        w
    //       ^ ^
    //      /   \.
    //     u --> v
    // There's a bijection between triangles of this form in `directed_graph`
    // and undirected triangles in `graph`.
    par_for(0, graph->n, [&](const size_t vertex_id) {
      auto vertex{directed_graph.get_vertex(vertex_id)};
      const uintT vertex_counter_offset{counter_offsets[vertex_id]};
      const auto intersect{[&](
          const uintE v_id,
          const uintE neighbor_id,
          const Weight weight,
          const uintE v_to_neighbor_index) {
        auto neighbor{directed_graph.get_vertex(neighbor_id)};
        const uintT neighbor_counter_offset{counter_offsets[neighbor_id]};
        const auto update_counters{[&](
            const uintE shared_neighbor,
            const uintE vertex_to_shared_index,
            const uintE neighbor_to_shared_index,
            const Weight weight_1,
            const Weight weight_2) {
          counters[vertex_counter_offset + vertex_to_shared_index] +=
            kWeightFactor * kWeightFactor * weight * weight_2;
          counters[neighbor_counter_offset + neighbor_to_shared_index] +=
            kWeightFactor * kWeightFactor * weight * weight_1;
        }};
        counters[vertex_counter_offset + v_to_neighbor_index] +=
          kWeightFactor * kWeightFactor *
            internal::intersect_f_with_index_par(
                &vertex, &neighbor, update_counters);
      }};
      constexpr bool kParallel{false};
      vertex.mapOutNghWithIndex(vertex_id, intersect, kParallel);
    });

    const pbbs::sequence<double> norms{
      graph->n,
      [&](const size_t vertex_id) {
        double norm = kWeightFactor * kWeightFactor;  // add self-loop weight
        auto vertex = graph->get_vertex(vertex_id);
        constexpr bool kParallel{false};
        const auto update_norm{
          [&](uintE, uintE, const Weight weight) {
            norm += kWeightFactor * kWeightFactor * weight * weight;
          }};
        vertex.mapOutNgh(vertex_id, update_norm, kParallel);
        return std::sqrt(norm);
      }};

    pbbs::sequence<EdgeSimilarity> similarities(graph->m);
    // Convert shared neighbor counts into similarities for each edge.
    par_for(0, directed_graph.n, [&](const size_t vertex_id) {
      const uintT v_counter_offset{counter_offsets[vertex_id]};
      const auto compute_similarity{[&](
          const uintE v_id,
          const uintE u_id,
          const Weight weight,
          const uintE v_to_u_index) {
        const uintT counter_index{v_counter_offset + v_to_u_index};
        // additional term is to account for self-loop edge in closed
        // neighborhood
        const double shared_weight{
          counters[counter_index]
            + 2 * kWeightFactor * kWeightFactor * weight};
        const float similarity{static_cast<float>(
            shared_weight / (norms[v_id] * norms[u_id]))};
        similarities[2 * counter_index] =
          {.source = v_id, .neighbor = u_id, .similarity = similarity};
        similarities[2 * counter_index + 1] =
          {.source = u_id, .neighbor = v_id, .similarity = similarity};
      }};
      directed_graph.get_vertex(vertex_id).mapOutNghWithIndex(
          vertex_id, compute_similarity);
    });

    directed_graph.del();
    return similarities;
  }
}

template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> JaccardSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph) const {
  constexpr auto similarity_func{
    [](const uintE neighborhood_size_1,
       const uintE neighborhood_size_2,
       const uintE num_shared_neighbors) {
    const uintE neighborhood_union{neighborhood_size_1 + neighborhood_size_2 -
      num_shared_neighbors};
    // The `+ 2` accounts for the Jaccard similarity being computed with respect
    // to closed neighborhoods.
    return static_cast<float>((num_shared_neighbors + 2)) / neighborhood_union;
  }};
  return internal::AllEdgeNeighborhoodSimilarities(graph, similarity_func);
}

template <template <typename> class VertexTemplate, typename Weight>
pbbs::sequence<EdgeSimilarity> ApproxCosineSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, Weight>* graph) const {
  const size_t degree_threshold{static_cast<size_t>(1.5 * num_samples_)};
  return internal::ApproxCosineEdgeSimilarities(
      graph, num_samples_, degree_threshold, random_seed_);
}

template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> ApproxJaccardSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph) const {
  const size_t degree_threshold{static_cast<size_t>(1.0 * num_samples_)};
  return internal::ApproxJaccardEdgeSimilarities(
      graph, num_samples_, degree_threshold, random_seed_);
}

}  // namespace scan

}  // namespace gbbs
