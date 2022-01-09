// Similarity measures for determining the similarity of adjacent pairs of
// vertices.
//
// TODO(tomtseng) needs refactoring -- this code has so much copy-and-paste code
// and is confusing
#pragma once

#include <algorithm>
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
//   sequence<EdgeSimilarity> AllEdges(Graph* graph) const;

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
  sequence<EdgeSimilarity> AllEdges(
      symmetric_graph<VertexTemplate, Weight>* graph) const;
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
  sequence<EdgeSimilarity> AllEdges(
      symmetric_graph<VertexTemplate, gbbs::empty>* graph) const;
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
  sequence<EdgeSimilarity> AllEdges(
      symmetric_graph<VertexTemplate, Weight>* graph) const;

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
  sequence<EdgeSimilarity> AllEdges(
      symmetric_graph<VertexTemplate, gbbs::empty>* graph) const;

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
struct EdgeCounter<gbbs::empty> {
  using type = uintE;  // unweighted case
};

// Compute (numerator / denominator), rounding up if there's any remainder.
// Numerator must be positive.
constexpr uint64_t DivideRoundingUp(const size_t numerator,
                                    const size_t denominator) {
  return (numerator - 1) / denominator + 1;
}

// Pseudorandomly generate `num_numbers` random normal numbers, each with zero
// mean and unit variance.
sequence<float> RandomNormalNumbers(size_t num_numbers, parlay::random rng);

// Create a directed version of `graph`, pointing edges from lower degree
// vertices to higher degree vertices. This upper bounds the out-degree of each
// vertex in the directed graph with `sqrt(graph->m)`.
template <template <typename> class VertexTemplate, typename Weight>
auto DirectGraphByDegree(symmetric_graph<VertexTemplate, Weight>* graph) {
  uintE* vertex_degree_ranking{rankNodes(*graph, graph->n)};
  const auto filter_predicate{[&](const uintE u, const uintE v, Weight) {
    return vertex_degree_ranking[u] < vertex_degree_ranking[v];
  }};
  auto directed_graph{filterGraph(*graph, filter_predicate)};
  gbbs::free_array(vertex_degree_ranking, graph->n);
  return directed_graph;
}

// Returns a sequence `vertex_offsets` such that if there is another sequence
// `edges` consisting of the out-edges of `*graph` sorted by source vertex, then
// `vertex_offsets[i]` is the first appearance of vertex i as a source vertex.
template <template <typename> class VertexTemplate, typename Weight>
sequence<uintT> VertexOutOffsets(
    symmetric_graph<VertexTemplate, Weight>* graph) {
  auto vertex_offsets = sequence<uintT>::from_function(
      graph->n,
      [&](const size_t i) { return graph->get_vertex(i).out_degree(); });
  parlay::scan_inplace(vertex_offsets);
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
sequence<EdgeSimilarity> AllEdgeNeighborhoodSimilarities(
    symmetric_graph<VertexTemplate, gbbs::empty>* graph,
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
  auto counters = sequence<std::atomic<uintE>>(directed_graph.m);
  parallel_for(0, directed_graph.m, [&](size_t i) { counters[i] = 0; });
  // We use `counter_offsets` to index into `counters` for each edge.
  const sequence<uintT> counter_offsets{
      internal::VertexOutOffsets(&directed_graph)};
  // Find triangles of the following form:
  //        w
  //       ^ ^
  //      /   \.
  //     u --> v
  // There's a bijection between triangles of this form in `directed_graph` and
  // undirected triangles in `graph`.
  parallel_for(0, graph->n, [&](const size_t vertex_id) {
    auto vertex{directed_graph.get_vertex(vertex_id)};
    const uintT vertex_counter_offset{counter_offsets[vertex_id]};
    const auto intersect{[&](const uintE v_id, const uintE neighbor_id,
                             gbbs::empty, const uintE v_to_neighbor_index) {
      auto neighbor{directed_graph.get_vertex(neighbor_id)};
      const uintT neighbor_counter_offset{counter_offsets[neighbor_id]};
      const auto update_counters{[&](const uintE shared_neighbor,
                                     const uintE vertex_to_shared_index,
                                     const uintE neighbor_to_shared_index) {
        counters[vertex_counter_offset + vertex_to_shared_index]++;
        counters[neighbor_counter_offset + neighbor_to_shared_index]++;
      }};
      counters[vertex_counter_offset + v_to_neighbor_index] +=
          internal::intersect_f_with_index_par(&vertex, &neighbor,
                                               update_counters);
    }};
    constexpr bool kParallel{false};
    vertex.out_neighbors().map_with_index(intersect, kParallel);
  });

  sequence<EdgeSimilarity> similarities(graph->m);
  // Convert shared neighbor counts into similarities for each edge.
  parallel_for(0, directed_graph.n, [&](const size_t vertex_id) {
    const uintT v_counter_offset{counter_offsets[vertex_id]};
    const uintE v_degree{graph->get_vertex(vertex_id).out_degree()};
    const auto compute_similarity{[&](const uintE v_id, const uintE u_id,
                                      gbbs::empty, const uintE v_to_u_index) {
      const uintT counter_index{v_counter_offset + v_to_u_index};
      const uintE num_shared_neighbors{counters[counter_index]};
      const uintE u_degree{graph->get_vertex(u_id).out_degree()};
      const float similarity{neighborhood_sizes_to_similarity(
          v_degree, u_degree, num_shared_neighbors)};
      similarities[2 * counter_index] = {
          .source = v_id, .neighbor = u_id, .similarity = similarity};
      similarities[2 * counter_index + 1] = {
          .source = u_id, .neighbor = v_id, .similarity = similarity};
    }};
    directed_graph.get_vertex(vertex_id).out_neighbors().map_with_index(
        compute_similarity);
  });

  return similarities;
}

// Implementation of ApproxCosineSimilarities::AllEdges.
//
// `degree_threshold` is a threshold so that we only approximate the similarity
// score between two vertices if their degrees are high enough. (When the
// degrees are low, it's cheap to compute the similarity exactly.)
template <template <typename> class VertexTemplate, typename Weight>
sequence<EdgeSimilarity> ApproxCosineEdgeSimilarities(
    symmetric_graph<VertexTemplate, Weight>* graph, const uint32_t num_samples,
    const size_t degree_threshold, const size_t random_seed) {
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
  sequence<uintE> needs_fingerprint_seq(graph->n, 0U);
  sequence<uintE> needs_normals_seq(graph->n, 0U);
  parallel_for(0, graph->n, [&](const size_t vertex_id) {
    Vertex vertex{graph->get_vertex(vertex_id)};
    if (vertex.out_degree() >= degree_threshold) {
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
      vertex.out_neighbors().map(check_degree_threshold);
      if (!skip_fingerprint) {
        needs_fingerprint_seq[vertex_id] = true;
        needs_normals_seq[vertex_id] = true;
        const auto set_needs_normals{
            [&](uintE, const uintE neighbor_id, Weight) {
              needs_normals_seq[neighbor_id] = true;
            }};
        vertex.out_neighbors().map(set_needs_normals);
      }
    }
  });
  // repurpose `needs_normals_seq` to serve as the index of a vertex into
  // `normals`, and same with `needs_fingerprint_seq`
  const uintE num_needs_normals{parlay::scan_inplace(needs_normals_seq)};
  const sequence<uintE>& normals_indices{needs_normals_seq};
  const sequence<float> normals{RandomNormalNumbers(
      num_needs_normals * num_samples, parlay::random{random_seed})};
  const uintE num_needs_fingerprint{
      parlay::scan_inplace(needs_fingerprint_seq)};
  const sequence<uintE>& fingerprint_indices{needs_fingerprint_seq};

  const size_t num_vertices{graph->n};
  const size_t num_bit_arrays{
      internal::DivideRoundingUp(num_samples, kBitArraySize)};
  // Simhash fingerprints.
  sequence<BitArray> fingerprints{sequence<BitArray>::uninitialized(
      num_needs_fingerprint * num_bit_arrays)};
  parallel_for(0, num_vertices, [&](const size_t vertex_id) {
    const uintE fingerprint_index{fingerprint_indices[vertex_id]};
    const bool needs_fingerprint{
        vertex_id + 1 == num_vertices
            ? fingerprint_index != num_needs_fingerprint
            : fingerprint_index != fingerprint_indices[vertex_id + 1]};
    if (!needs_fingerprint) {
      return;
    }
    Vertex vertex{graph->get_vertex(vertex_id)};
    const uintE vertex_normal_offset{num_samples * normals_indices[vertex_id]};
    const size_t fingerprint_offset{fingerprint_index * num_bit_arrays};
    parallel_for(0, num_bit_arrays, [&](const size_t bit_array_id) {
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
      const auto update_dot_products{[&](uintE, const uintE neighbor_id,
                                         [[maybe_unused]] const Weight weight) {
        for (size_t bit_id = 0; bit_id < max_bit_id; bit_id++) {
          if
            constexpr(std::is_same<Weight, gbbs::empty>::value) {
              // unweighted case
              hyperplane_dot_products[bit_id] +=
                  normals[num_samples * normals_indices[neighbor_id] +
                          bits_offset + bit_id];
            }
          else {  // weighted case
            hyperplane_dot_products[bit_id] +=
                normals[num_samples * normals_indices[neighbor_id] +
                        bits_offset + bit_id] *
                weight;
          }
        }
      }};
      constexpr bool kParallel{false};
      vertex.out_neighbors().map(update_dot_products, kParallel);
      for (size_t bit_id = 0; bit_id < max_bit_id; bit_id++) {
        if (hyperplane_dot_products[bit_id] >= 0) {
          bits |= (1LL << bit_id);
        }
      }
      fingerprints[fingerprint_offset + bit_array_id] = bits;
    });
  });

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
  auto counters = sequence<std::atomic<Counter>>(directed_graph.m);
  parallel_for(0, directed_graph.m, [&](size_t i) { counters[i] = 0; });
  // For floating-point-weighted graphs, we multiply all edge weights by this
  // factor and round to the nearest integer. Multiplying all edge weights by a
  // constant doesn't affect cosine similarities. The multiplication preserves
  // more accuracy after rounding.
  constexpr int64_t kWeightFactor =
      std::is_floating_point<Weight>::value ? 1000 : 1;
  // We use `counter_offsets` to index into `counters` for each edge.
  const sequence<uintT> counter_offsets{
      internal::VertexOutOffsets(&directed_graph)};
  // Find triangles of the following form:
  //        w
  //       ^ ^
  //      /   \.
  //     u --> v
  // Count each of these triangles to get the number of shared neighbors between
  // vertices. However, we skip pairs of vertices that have high degree in the
  // original, undirected graph.
  parallel_for(0, graph->n, [&](const size_t vertex_id) {
    auto vertex{directed_graph.get_vertex(vertex_id)};
    const bool vertex_is_high_degree{graph->v_data[vertex_id].degree >=
                                     degree_threshold};
    if (vertex_is_high_degree) {
      // Since all edges in the directed graph point towards higher degree
      // vertices, if the current vertex is high degree, so are all its directed
      // neighbors. Skip these pairs of edges since we'll approximate their
      // similarities.
      return;
    }

    const uintT vertex_counter_offset{counter_offsets[vertex_id]};
    const auto intersect{[&](const uintE v_id, const uintE neighbor_id,
                             [[maybe_unused]] const Weight weight,
                             const uintE v_to_neighbor_index) {
      auto neighbor{directed_graph.get_vertex(neighbor_id)};
      const bool neighbor_is_high_degree{graph->v_data[neighbor_id].degree >=
                                         degree_threshold};
      const uintT neighbor_counter_offset{counter_offsets[neighbor_id]};

      if
        constexpr(std::is_same<Weight, gbbs::empty>::value) {
          // unweighted case
          const auto update_counters{[&](const uintE shared_neighbor,
                                         const uintE vertex_to_shared_index,
                                         const uintE neighbor_to_shared_index) {
            counters[vertex_counter_offset + vertex_to_shared_index]++;
            if (!(neighbor_is_high_degree &&
                  graph->v_data[shared_neighbor].degree >= degree_threshold)) {
              counters[neighbor_counter_offset + neighbor_to_shared_index]++;
            }
          }};
          counters[vertex_counter_offset + v_to_neighbor_index] +=
              internal::intersect_f_with_index_par(&vertex, &neighbor,
                                                   update_counters);
        }
      else {
        // weighted case
        const auto update_counters{[&](
            const uintE shared_neighbor, const uintE vertex_to_shared_index,
            const uintE neighbor_to_shared_index, const Weight weight_1,
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
            internal::intersect_f_with_index_par(&vertex, &neighbor,
                                                 update_counters);
      }
    }};
    constexpr bool kParallel{false};
    vertex.out_neighbors().map_with_index(intersect, kParallel);
  });

  const sequence<double> norms{[&]() -> sequence<double> {
    if
      constexpr(std::is_same<Weight, gbbs::empty>::value) {  // unweighted
        return {};
      }
    else {  // weighted case
      return sequence<double>::from_function(
          graph->n, [&](const size_t vertex_id) {
            double norm =
                kWeightFactor * kWeightFactor;  // add self-loop weight
            auto vertex = graph->get_vertex(vertex_id);
            constexpr bool kParallel{false};
            const auto update_norm{[&](uintE, uintE, const Weight weight) {
              norm += kWeightFactor * kWeightFactor * weight * weight;
            }};
            vertex.out_neighbors().map(update_norm, kParallel);
            return std::sqrt(norm);
          });
    }
  }()};

  sequence<EdgeSimilarity> similarities(graph->m);
  // Convert shared neighbor counts into similarities for each edge.
  parallel_for(0, directed_graph.n, [&](const size_t vertex_id) {
    const uintT v_counter_offset{counter_offsets[vertex_id]};
    const uintE v_degree{graph->get_vertex(vertex_id).out_degree()};
    const bool vertex_is_high_degree{v_degree >= degree_threshold};
    const size_t vertex_fingerprint_offset{fingerprint_indices[vertex_id] *
                                           num_bit_arrays};
    const auto compute_similarity{[&](const uintE v_id, const uintE u_id,
                                      [[maybe_unused]] const Weight weight,
                                      const uintE v_to_u_index) {
      const uintT counter_index{v_counter_offset + v_to_u_index};
      float similarity{-1};
      if (vertex_is_high_degree) {  // approximate similarity
        const size_t neighbor_fingerprint_offset{fingerprint_indices[u_id] *
                                                 num_bit_arrays};
        const auto fingerprint_xor{
            parlay::delayed_seq<std::remove_const<decltype(num_samples)>::type>(
                num_bit_arrays, [&](const size_t i) {
                  return __builtin_popcountll(
                      fingerprints[vertex_fingerprint_offset + i] ^
                      fingerprints[neighbor_fingerprint_offset + i]);
                })};
        const float angle_estimate{static_cast<float>(
            parlay::reduce(fingerprint_xor) * M_PI / num_samples)};
        similarity = std::max(std::cos(angle_estimate), 0.0f);
      } else {  // exact similarity
        if
          constexpr(std::is_same<Weight, gbbs::empty>::value) {
            // unweighted case
            const uintE num_shared_neighbors{counters[counter_index]};
            const uintE u_degree{graph->get_vertex(u_id).out_degree()};
            similarity = (num_shared_neighbors + 2) /
                         (sqrtf(v_degree + 1) * sqrtf(u_degree + 1));
          }
        else {  // weighted case
          // additional term is to account for self-loop edge in closed
          // neighborhood
          const double shared_weight{
              static_cast<double>(counters[counter_index] +
                                  2 * kWeightFactor * kWeightFactor * weight)};
          similarity =
              static_cast<float>(shared_weight / (norms[v_id] * norms[u_id]));
        }
      }
      similarities[2 * counter_index] = {
          .source = v_id, .neighbor = u_id, .similarity = similarity};
      similarities[2 * counter_index + 1] = {
          .source = u_id, .neighbor = v_id, .similarity = similarity};
    }};
    directed_graph.get_vertex(vertex_id).out_neighbors().map_with_index(
        compute_similarity);
  });

  return similarities;
}

// Implementation of ApproxJaccardSimilarities::AllEdges.
//
// `degree_threshold` is a threshold so that we only approximate the similarity
// score between two vertices if their degrees are high enough. (When the
// degrees are low, it's cheap to compute the similarity exactly.)
template <template <typename> class VertexTemplate>
sequence<EdgeSimilarity> ApproxJaccardEdgeSimilarities(
    symmetric_graph<VertexTemplate, gbbs::empty>* graph,
    const uint32_t original_num_samples, const size_t degree_threshold,
    const size_t random_seed) {
  using Weight = gbbs::empty;
  using Vertex = VertexTemplate<Weight>;
  // For edges between high degree vertices, estimate the Jaccard similarity
  // with a MinHash variant --- see paper "One Permutation Hashing for Efficient
  // Search and Learning."
  // The MinHash variant works as follows to estimate the Jaccard similarity
  // between two sets using k samples: partition the universe U of elements into
  // k equally sized buckets U_1, U_2, ..., U_k (e.g., permute U to assign each
  // element x some value permute(x) in [0, |U|), then assign x to bucket
  // U_{(permute(x) % k) + 1}). The fingerprint of a set S consists of (x_1,
  // ..., x_k) where x_i is the minimum element of U_i that is also in S, or
  // <empty> if no such element exists. Then the estimate of the Jaccard
  // distance between sets with fingerprints (x_1, ..., x_k) and (y_1, ..., y_k)
  // is <number of indices i where x_i == y_i != <empty>> / (k - <number of
  // indices where x_i == y_i == <empty>>).
  //
  // For edges with a low degree vertex, compute the Jaccard similarity exactly
  // with triangle counting like in `AllEdgeNeighborhoodSimilarities()`.
  const uint32_t log_num_samples{
      std::max<uint32_t>(parlay::log2_up(original_num_samples), 1)};
  const uintE num_samples{static_cast<uintE>(1ULL << log_num_samples)};
  const uintE bucket_mask{num_samples - 1};
  constexpr uintE kEmptyBucket{UINT_E_MAX};
  const size_t num_vertices{graph->n};
  const sequence<uintE> vertex_permutation{parlay::random_permutation<uintE>(
      num_vertices, parlay::random{random_seed})};

  auto needs_fingerprint_seq =
      sequence<uintE>::from_function(graph->n, [&](const size_t vertex_id) {
        Vertex vertex{graph->get_vertex(vertex_id)};
        if (vertex.out_degree() < degree_threshold) {
          return false;
        }
        bool needs_fingerprint{false};
        const auto check_degree_threshold{
            [&](uintE, const uintE neighbor_id, Weight) {
              if (!needs_fingerprint &&
                  graph->v_data[neighbor_id].degree >= degree_threshold) {
                needs_fingerprint = true;
              }
            }};
        vertex.out_neighbors().map(check_degree_threshold);
        return needs_fingerprint;
      });
  const uintE num_needs_fingerprint{
      parlay::scan_inplace(needs_fingerprint_seq)};
  const sequence<uintE>& fingerprint_indices{needs_fingerprint_seq};
  // Compute MinHash fingerprints for high degree vertices.
  sequence<uintE> fingerprints(num_needs_fingerprint * num_samples,
                               kEmptyBucket);
  parallel_for(0, num_vertices, [&](const size_t vertex_id) {
    const uintE fingerprint_index{fingerprint_indices[vertex_id]};
    const bool needs_fingerprint{
        vertex_id + 1 == num_vertices
            ? fingerprint_index != num_needs_fingerprint
            : fingerprint_index != fingerprint_indices[vertex_id + 1]};
    if (!needs_fingerprint) {
      return;
    }

    Vertex vertex{graph->get_vertex(vertex_id)};
    const size_t fingerprint_offset{fingerprint_index * num_samples};
    const auto update_fingerprint{
        [&](uintE, const uintE neighbor, gbbs::empty) {
          const uintE permuted_neighbor{vertex_permutation[neighbor]};
          const uintE bucket_id{permuted_neighbor & bucket_mask};
          const uintE bucket_value{permuted_neighbor >> log_num_samples};
          gbbs::write_min(&(fingerprints[fingerprint_offset + bucket_id]),
                          bucket_value, std::less<uintE>{});
        }};
    update_fingerprint(vertex_id, vertex_id, gbbs::empty{});
    vertex.out_neighbors().map(update_fingerprint);
  });

  auto directed_graph{DirectGraphByDegree(graph)};
  // Each counter in `counters` holds the number of shared neighbors in `graph`
  // between u and v for some edge {u, v}.
  auto counters = sequence<std::atomic<uintE>>(directed_graph.m);
  parallel_for(0, directed_graph.m, [&](size_t i) { counters[i] = 0; });
  // We use `counter_offsets` to index into `counters` for each edge.
  const sequence<uintT> counter_offsets{
      internal::VertexOutOffsets(&directed_graph)};
  // Find triangles of the following form:
  //        w
  //       ^ ^
  //      /   \.
  //     u --> v
  // Count each of these triangles to get the number of shared neighbors between
  // vertices. However, we skip pairs of vertices that have high degree in the
  // original, undirected graph.
  parallel_for(0, graph->n, [&](const size_t vertex_id) {
    auto vertex{directed_graph.get_vertex(vertex_id)};
    const bool vertex_is_high_degree{graph->v_data[vertex_id].degree >=
                                     degree_threshold};
    if (vertex_is_high_degree) {
      // Since all edges in the directed graph point towards higher degree
      // vertices, if the current vertex is high degree, so are all its directed
      // neighbors. Skip these pairs of edges since we'll approximate their
      // similarities.
      return;
    }

    const uintT vertex_counter_offset{counter_offsets[vertex_id]};
    const auto intersect{[&](const uintE v_id, const uintE neighbor_id,
                             gbbs::empty, const uintE v_to_neighbor_index) {
      auto neighbor{directed_graph.get_vertex(neighbor_id)};
      const bool neighbor_is_high_degree{graph->v_data[neighbor_id].degree >=
                                         degree_threshold};
      const uintT neighbor_counter_offset{counter_offsets[neighbor_id]};
      const auto update_counters{[&](const uintE shared_neighbor,
                                     const uintE vertex_to_shared_index,
                                     const uintE neighbor_to_shared_index) {
        counters[vertex_counter_offset + vertex_to_shared_index]++;
        if (!(neighbor_is_high_degree &&
              graph->v_data[shared_neighbor].degree >= degree_threshold)) {
          counters[neighbor_counter_offset + neighbor_to_shared_index]++;
        }
      }};
      counters[vertex_counter_offset + v_to_neighbor_index] +=
          internal::intersect_f_with_index_par(&vertex, &neighbor,
                                               update_counters);
    }};
    constexpr bool kParallel{false};
    vertex.out_neighbors().map_with_index(intersect, kParallel);
  });

  sequence<EdgeSimilarity> similarities(graph->m);
  // Convert shared neighbor counts into similarities for each edge.
  parallel_for(0, directed_graph.n, [&](const size_t vertex_id) {
    const uintT v_counter_offset{counter_offsets[vertex_id]};
    const uintE v_degree{graph->get_vertex(vertex_id).out_degree()};
    const bool vertex_is_high_degree{v_degree >= degree_threshold};
    const size_t vertex_fingerprint_offset{fingerprint_indices[vertex_id] *
                                           num_samples};
    const auto compute_similarity{[&](const uintE v_id, const uintE u_id,
                                      gbbs::empty, const uintE v_to_u_index) {
      const uintT counter_index{v_counter_offset + v_to_u_index};
      float similarity{-1};
      if (vertex_is_high_degree) {  // approximate similarity
        const size_t neighbor_fingerprint_offset{fingerprint_indices[u_id] *
                                                 num_samples};
        const uintE fingerprint_matches{parlay::reduce(
            parlay::delayed_seq<uintE>(num_samples, [&](const size_t i) {
              return fingerprints[vertex_fingerprint_offset + i] !=
                         kEmptyBucket &&
                     fingerprints[vertex_fingerprint_offset + i] ==
                         fingerprints[neighbor_fingerprint_offset + i];
            }))};
        const uintE fingerprint_empty_count{parlay::reduce(
            parlay::delayed_seq<uintE>(num_samples, [&](const size_t i) {
              return fingerprints[vertex_fingerprint_offset + i] ==
                         kEmptyBucket &&
                     fingerprints[neighbor_fingerprint_offset + i] ==
                         kEmptyBucket;
            }))};
        similarity =
            fingerprint_matches /
            (static_cast<float>(num_samples - fingerprint_empty_count));
      } else {  // exact similarity
        const uintE num_shared_neighbors{counters[counter_index]};
        const uintE u_degree{graph->get_vertex(u_id).out_degree()};
        const uintE neighborhood_union{v_degree + u_degree -
                                       num_shared_neighbors};
        similarity =
            static_cast<float>((num_shared_neighbors + 2)) / neighborhood_union;
      }
      similarities[2 * counter_index] = {
          .source = v_id, .neighbor = u_id, .similarity = similarity};
      similarities[2 * counter_index + 1] = {
          .source = u_id, .neighbor = v_id, .similarity = similarity};
    }};
    directed_graph.get_vertex(vertex_id).out_neighbors().map_with_index(
        compute_similarity);
  });

  return similarities;
}

}  // namespace internal

template <template <typename> class VertexTemplate, typename Weight>
sequence<EdgeSimilarity> CosineSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, Weight>* graph) const {
  if
    constexpr(std::is_same<Weight, gbbs::empty>::value) {  // unweighted
      constexpr auto similarity_func{[](const uintE neighborhood_size_1,
                                        const uintE neighborhood_size_2,
                                        const uintE num_shared_neighbors) {
        // SCAN structural/cosine similarities are defined using _closed_
        // neighborhoods, hence the need to to adjust these values by `+ 1` and
        // `+ 2`.
        return (num_shared_neighbors + 2) / (sqrtf(neighborhood_size_1 + 1) *
                                             sqrtf(neighborhood_size_2 + 1));
      }};
      return internal::AllEdgeNeighborhoodSimilarities(graph, similarity_func);
    }
  else {  // weighted case
    auto directed_graph{internal::DirectGraphByDegree(graph)};
    // Each counter in `counters` will hold numerator of CosineSimilarity(u, v)
    // for some edge {u, v}.
    auto counters = sequence<std::atomic<uintE>>(directed_graph.m);
    parallel_for(0, directed_graph.m, [&](size_t i) { counters[i] = 0; });
    // For floating-point-weighted graphs, we multiply all edge weights by this
    // factor and round to the nearest integer. Multiplying all edge weights by
    // a constant doesn't affect cosine similarities. The multiplication
    // preserves more accuracy after rounding.
    constexpr int64_t kWeightFactor =
        std::is_floating_point<Weight>::value ? 1000 : 1;
    // We use `counter_offsets` to index into `counters` for each edge.
    const sequence<uintT> counter_offsets{
        internal::VertexOutOffsets(&directed_graph)};
    // Find triangles of the following form:
    //        w
    //       ^ ^
    //      /   \.
    //     u --> v
    // There's a bijection between triangles of this form in `directed_graph`
    // and undirected triangles in `graph`.
    parallel_for(0, graph->n, [&](const size_t vertex_id) {
      auto vertex{directed_graph.get_vertex(vertex_id)};
      const uintT vertex_counter_offset{counter_offsets[vertex_id]};
      const auto intersect{[&](const uintE v_id, const uintE neighbor_id,
                               const Weight weight,
                               const uintE v_to_neighbor_index) {
        auto neighbor{directed_graph.get_vertex(neighbor_id)};
        const uintT neighbor_counter_offset{counter_offsets[neighbor_id]};
        const auto update_counters{
            [&](const uintE shared_neighbor, const uintE vertex_to_shared_index,
                const uintE neighbor_to_shared_index, const Weight weight_1,
                const Weight weight_2) {
              counters[vertex_counter_offset + vertex_to_shared_index] +=
                  kWeightFactor * kWeightFactor * weight * weight_2;
              counters[neighbor_counter_offset + neighbor_to_shared_index] +=
                  kWeightFactor * kWeightFactor * weight * weight_1;
            }};
        counters[vertex_counter_offset + v_to_neighbor_index] +=
            kWeightFactor * kWeightFactor *
            internal::intersect_f_with_index_par(&vertex, &neighbor,
                                                 update_counters);
      }};
      constexpr bool kParallel{false};
      vertex.out_neighbors().map_with_index(intersect, kParallel);
    });

    const auto norms{
        sequence<double>::from_function(graph->n, [&](const size_t vertex_id) {
          double norm = kWeightFactor * kWeightFactor;  // add self-loop weight
          auto vertex = graph->get_vertex(vertex_id);
          constexpr bool kParallel{false};
          const auto update_norm{[&](uintE, uintE, const Weight weight) {
            norm += kWeightFactor * kWeightFactor * weight * weight;
          }};
          vertex.out_neighbors().map(update_norm, kParallel);
          return std::sqrt(norm);
        })};

    sequence<EdgeSimilarity> similarities(graph->m);
    // Convert shared neighbor counts into similarities for each edge.
    parallel_for(0, directed_graph.n, [&](const size_t vertex_id) {
      const uintT v_counter_offset{counter_offsets[vertex_id]};
      const auto compute_similarity{[&](const uintE v_id, const uintE u_id,
                                        const Weight weight,
                                        const uintE v_to_u_index) {
        const uintT counter_index{v_counter_offset + v_to_u_index};
        // additional term is to account for self-loop edge in closed
        // neighborhood
        const double shared_weight{
            static_cast<double>(counters[counter_index] +
                                2 * kWeightFactor * kWeightFactor * weight)};
        const float similarity{
            static_cast<float>(shared_weight / (norms[v_id] * norms[u_id]))};
        similarities[2 * counter_index] = {
            .source = v_id, .neighbor = u_id, .similarity = similarity};
        similarities[2 * counter_index + 1] = {
            .source = u_id, .neighbor = v_id, .similarity = similarity};
      }};
      directed_graph.get_vertex(vertex_id).out_neighbors().map_with_index(
          compute_similarity);
    });

    return similarities;
  }
}

template <template <typename> class VertexTemplate>
sequence<EdgeSimilarity> JaccardSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, gbbs::empty>* graph) const {
  constexpr auto similarity_func{[](const uintE neighborhood_size_1,
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
sequence<EdgeSimilarity> ApproxCosineSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, Weight>* graph) const {
  const size_t degree_threshold{static_cast<size_t>(1.5 * num_samples_)};
  return internal::ApproxCosineEdgeSimilarities(graph, num_samples_,
                                                degree_threshold, random_seed_);
}

template <template <typename> class VertexTemplate>
sequence<EdgeSimilarity> ApproxJaccardSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, gbbs::empty>* graph) const {
  const size_t degree_threshold{static_cast<size_t>(1.0 * num_samples_)};
  return internal::ApproxJaccardEdgeSimilarities(
      graph, num_samples_, degree_threshold, random_seed_);
}

}  // namespace scan

}  // namespace gbbs
