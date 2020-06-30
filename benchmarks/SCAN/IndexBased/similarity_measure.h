// Similarity measures for determining the similarity of adjacent pairs of
// vertices.
#pragma once

#include <atomic>
#include <cmath>
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
class CosineSimilarity {
 public:
  CosineSimilarity() = default;

  template <template <typename> class VertexTemplate>
  pbbs::sequence<EdgeSimilarity>
  AllEdges(symmetric_graph<VertexTemplate, pbbs::empty>* graph) const;
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
// cosine similarity with absolute error up to `b`.
//
// This is a biased estimate of the cosine similarity.
struct ApproxCosineSimilarity {
 public:
  ApproxCosineSimilarity(uint32_t num_samples, size_t random_seed);

  // When `random_seed` is fixed, the output of `AllEdges` is deterministic.
  template <template <typename> class VertexTemplate>
  pbbs::sequence<EdgeSimilarity>
  AllEdges(symmetric_graph<VertexTemplate, pbbs::empty>* graph) const;

 private:
  const uint32_t num_samples_;
  const size_t random_seed_;
};

// This is an approximate version of `JaccardSimilarity`. Increasing
// `num_samples` increases the approximation accuracy.
//
// Let `m` be the number of undirected edges in the graph, and let `a` and `b`
// be in the range (0, 1). Then, if we replace the pseudorandom number generator
// used within the code with perfectly random number generator and replace the
// hash function with a random hash function with no collisions, then picking
//   num_samples = 3 * ln(2 * m / a) / b ^ 2
// gives that with probability at least `1 - a`, each edge receives the correct
// Jaccard similarity with absolute error up to `b`.
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

pbbs::sequence<float>
RandomNormalNumbers(size_t num_numbers, pbbs::random rng);

// Compute (numerator / denominator), rounding up if there's any remainder.
// Numerator must be positive.
constexpr uint64_t
DivideRoundingUp(const size_t numerator, const size_t denominator) {
  return (numerator - 1) / denominator + 1;
}

pbbs::sequence<EdgeSimilarity> BidirectionalSimilarities(
    const size_t num_directed_edges,
    const pbbs::sequence<EdgeSimilarity>& unidirectional_similarities);

// Returns an sequence `vertex_offsets` such that if there is another sequence
// `edges` consisting of the out-edges of `*graph` sorted by source vertex, then
// `vertex_offsets[i]` is the first appearance of vertex i as a source vertex.
template <class Graph>
pbbs::sequence<uintT> VertexOutOffsets(Graph* graph) {
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
//     the intersection of the two negibhorhoods) as arguments. The function
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
  const pbbs::sequence<uintT> counter_offsets{
    internal::VertexOutOffsets(&directed_graph)};
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
  // Convert shared neighbor counts into similarities for each edge.
  par_for(0, directed_graph.n, [&](const size_t vertex_id) {
    const uintT v_counter_offset{counter_offsets[vertex_id]};
    const uintE v_neighborhood{graph->get_vertex(vertex_id).getOutDegree()};
    const auto compute_similarity{[&](
        const uintE v_id,
        const uintE u_id,
        pbbs::empty,
        const uintE v_to_u_index) {
      const uintT counter_index{v_counter_offset + v_to_u_index};
      const uintE num_shared_neighbors{counters[counter_index]};
      const uintE u_neighborhood{graph->get_vertex(u_id).getOutDegree()};
      const float similarity{neighborhood_sizes_to_similarity(
          v_neighborhood, u_neighborhood, num_shared_neighbors)};
      similarities[2 * counter_index] =
        {.source = v_id, .neighbor = u_id, .similarity = similarity};
      similarities[2 * counter_index + 1] =
        {.source = u_id, .neighbor = v_id, .similarity = similarity};
    }};
    directed_graph.get_vertex(vertex_id).mapOutNghWithIndex(
        vertex_id, compute_similarity);
  });

  directed_graph.del();
  pbbs::free_array(vertex_degree_ranking);
  return similarities;
}

// Implementation of ApproxCosineSimilarities::AllEdges.
//
// `exact_threshold` is a threshold where if the sum of neighborhood sizes of
// adjacent vertices u and v is below `exact_threshold`, then we compute
// similarity score between u and v exactly.
// (The cost to computing the similarity score between vertices u and v exactly
// is O(<size of neighborhood of u> + <size of neighborhood of v>), whereas the
// cost to computing the similarity score approximately is O(num_samples_).)
template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> ApproxCosineEdgeSimilarities(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph,
    const uint32_t num_samples,
    const size_t exact_threshold,
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

  using Weight = pbbs::empty;
  using Vertex = VertexTemplate<Weight>;
  // We compute `num_samples_` hyperplanes and, to sketch a vertex's
  // neighborhood vector, we compute `num_samples_` bits representing the sign
  // of the vector's dot product with each hyperplane. For efficiency, we store
  // the bits in chunks of `kBitArraySize` rather than one-by-one.
  using BitArray = uint64_t;
  constexpr size_t kBitArraySize{sizeof(BitArray) * 8};

  const size_t num_vertices{graph->n};
  // <number of vertices>-by-<number of samples> array of normal numbers
  const pbbs::sequence<float> normals{internal::RandomNormalNumbers(
      num_vertices * num_samples, pbbs::random{random_seed})};
  const auto addition_monoid{pbbs::addm<float>{}};
  const size_t num_bit_arrays{
    internal::DivideRoundingUp(num_samples, kBitArraySize)};
  const pbbs::sequence<pbbs::sequence<BitArray>> vertex_fingerprints{
    num_vertices,
    [&](const size_t vertex_id) {
      Vertex vertex{graph->get_vertex(vertex_id)};
      const uintE degree{vertex.getOutDegree()};
      bool skip_fingerprint{true};
      const auto check_exact_threshold{
        [&](uintE, const uintE neighbor_id, Weight) {
          if (skip_fingerprint &&
              (degree + graph->get_vertex(neighbor_id).getOutDegree() >=
               exact_threshold)) {
            skip_fingerprint = false;
          }
        }};
      vertex.mapOutNgh(vertex_id, check_exact_threshold);
      if (skip_fingerprint) {
        return pbbs::sequence<BitArray>{};
      }
      return pbbs::sequence<BitArray>{
        num_bit_arrays,
        [&](const size_t bit_array_id) {
          BitArray bits{0};
          const size_t max_bit_id{
            bit_array_id == num_bit_arrays - 1
            ? kBitArraySize - (num_bit_arrays * kBitArraySize - num_samples)
            : kBitArraySize};
          const size_t bits_offset{bit_array_id * kBitArraySize};
          for (size_t bit_id = 0; bit_id < max_bit_id; bit_id++) {
            const auto neighbor_to_normal_sample{
              [&](uintE, const uintE neighbor, pbbs::empty) {
                return normals[num_samples * neighbor + bits_offset + bit_id];
            }};
            const float hyperplane_dot_product{
              normals[num_samples * vertex_id + bits_offset + bit_id] +
                vertex.template reduceOutNgh<float>(
                  vertex_id, neighbor_to_normal_sample, addition_monoid)};
            if (hyperplane_dot_product >= 0) {
              bits |= (1LL << bit_id);
            }
          }
          return bits;
        }};
    }};

  pbbs::sequence<EdgeSimilarity> undirected_similarities{
    pbbs::sequence<EdgeSimilarity>::no_init(graph->m)};
  par_for(0, undirected_similarities.size(), [&](const size_t i) {
    undirected_similarities[i].similarity =
      std::numeric_limits<float>::quiet_NaN();
  });
  const pbbs::sequence<uintT> vertex_offsets{internal::VertexOutOffsets(graph)};
  // Get similarities for all edges (u, v) where u < v.
  const auto compute_similarity{[&](
      const uintE vertex_id,
      const uintE neighbor_id,
      pbbs::empty,
      const uintE neighbor_index) {
    if (vertex_id < neighbor_id) {
      float similarity_estimate{0};
      Vertex vertex{graph->get_vertex(vertex_id)};
      Vertex neighbor{graph->get_vertex(neighbor_id)};
      const size_t vertex_degree{vertex.getOutDegree()};
      const size_t neighbor_degree{neighbor.getOutDegree()};
      if (vertex_degree + neighbor_degree < exact_threshold) {
        // compute exact similarity
        const size_t num_shared_neighbors{
          vertex.intersect(&neighbor, vertex_id, neighbor_id)};
        similarity_estimate = (num_shared_neighbors + 2) /
          (sqrtf(vertex_degree + 1) * sqrtf(neighbor_degree + 1));
      } else {  // compute approximate similarity
        const pbbs::sequence<BitArray>& vertex_fingerprint{
          vertex_fingerprints[vertex_id]};
        const pbbs::sequence<BitArray>& neighbor_fingerprint{
          vertex_fingerprints[neighbor_id]};
        const auto fingerprint_xor{
          pbbs::delayed_seq<std::remove_const<decltype(num_samples)>::type>(
            vertex_fingerprint.size(),
            [&](const size_t i) {
              return __builtin_popcountll(
                  vertex_fingerprint[i] ^ neighbor_fingerprint[i]);
            })};
        const float angle_estimate{static_cast<float>(
            pbbslib::reduce_add(fingerprint_xor) * M_PI / num_samples)};
        similarity_estimate = std::cos(angle_estimate);
      }
      undirected_similarities[vertex_offsets[vertex_id] + neighbor_index] = {
        .source = vertex_id,
        .neighbor = neighbor_id,
        .similarity = similarity_estimate};
    }
  }};
  par_for(0, num_vertices, [&](const size_t vertex_id) {
    graph->get_vertex(vertex_id).mapOutNghWithIndex(
        vertex_id, compute_similarity);
  });
  return internal::BidirectionalSimilarities(graph->m, undirected_similarities);
}

// Implementation of ApproxJaccardSimilarities::AllEdges.
//
// `exact_threshold` is a threshold where if the sum of neighborhood sizes of
// adjacent vertices u and v is below `exact_threshold`, then we compute
// similarity score between u and v exactly.
// (The cost to computing the similarity score between vertices u and v exactly
// is O(<size of neighborhood of u> + <size of neighborhood of v>), whereas the
// cost to computing the similarity score approximately is O(num_samples_).)
template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> ApproxJaccardEdgeSimilarities(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph,
    const uint32_t num_samples,
    const size_t exact_threshold,
    const size_t random_seed) {
  using Weight = pbbs::empty;
  using Vertex = VertexTemplate<Weight>;
  // Estimate the Jaccard similarity with MinHash.

  const size_t num_vertices{graph->n};
  const auto min_monoid{pbbs::minm<uint64_t>{}};
  const uint64_t random_offset{pbbs::hash64(random_seed)};
  const pbbs::sequence<pbbs::sequence<uint64_t>> vertex_fingerprints{
    num_vertices,
    [&](const size_t vertex_id) {
      Vertex vertex{graph->get_vertex(vertex_id)};
      if (vertex.getOutDegree() < exact_threshold) {
        return pbbs::sequence<uint64_t>{};
      }
      bool skip_fingerprint{true};
      const auto check_exact_threshold{
        [&](uintE, const uintE neighbor_id, Weight) {
          if (skip_fingerprint &&
              graph->get_vertex(neighbor_id).getOutDegree() >=
                exact_threshold) {
            skip_fingerprint = false;
          }
        }};
      vertex.mapOutNgh(vertex_id, check_exact_threshold);
      if (skip_fingerprint) {
        return pbbs::sequence<uint64_t>{};
      }
      return pbbs::sequence<uint64_t>{
        num_samples,
        [&](const size_t sample_id) {
          const auto hash_neighbor{
            [&](uintE, const uintE neighbor, pbbs::empty) {
              return pbbs::hash64_2(
                  random_offset + num_samples * neighbor + sample_id);
          }};
          return std::min(
              pbbs::hash64_2(
                random_offset + num_samples * vertex_id + sample_id),
              vertex.template reduceOutNgh<uint64_t>(
                vertex_id, hash_neighbor, min_monoid));
        }};
    }};

  pbbs::sequence<EdgeSimilarity> undirected_similarities{
    pbbs::sequence<EdgeSimilarity>::no_init(graph->m)};
  par_for(0, undirected_similarities.size(), [&](const size_t i) {
    undirected_similarities[i].similarity =
      std::numeric_limits<float>::quiet_NaN();
  });
  const pbbs::sequence<uintT> vertex_offsets{internal::VertexOutOffsets(graph)};
  // Get similarities for all edges (u, v) where u < v.
  const auto compute_similarity{[&](
      const uintE vertex_id,
      const uintE neighbor_id,
      pbbs::empty,
      const uintE neighbor_index) {
    if (vertex_id < neighbor_id) {
      float similarity_estimate{0};
      Vertex vertex{graph->get_vertex(vertex_id)};
      Vertex neighbor{graph->get_vertex(neighbor_id)};
      const size_t vertex_degree{vertex.getOutDegree()};
      const size_t neighbor_degree{neighbor.getOutDegree()};
      if (vertex_degree < exact_threshold ||
          neighbor_degree < exact_threshold) {
        // compute exact similarity
        constexpr auto no_op{[](uintE, uintE, uintE) {}};
        const size_t num_shared_neighbors{
          vertex.intersect_f_par(&neighbor, vertex_id, neighbor_id, no_op)};
        const size_t union_of_neighbors{
          vertex_degree + neighbor_degree - num_shared_neighbors};
        similarity_estimate =
          static_cast<float>(num_shared_neighbors + 2) /
            union_of_neighbors;
      } else {
        const pbbs::sequence<uint64_t>& vertex_fingerprint{
          vertex_fingerprints[vertex_id]};
        const pbbs::sequence<uint64_t>& neighbor_fingerprint{
          vertex_fingerprints[neighbor_id]};
        const auto fingerprint_matches{
          pbbs::delayed_seq<std::remove_const<decltype(num_samples)>::type>(
            vertex_fingerprint.size(),
            [&](const size_t i) {
              return vertex_fingerprint[i] == neighbor_fingerprint[i];
            })};
        similarity_estimate =
            pbbslib::reduce_add(fingerprint_matches) /
              static_cast<float>(num_samples);
      }
      undirected_similarities[vertex_offsets[vertex_id] + neighbor_index] = {
        .source = vertex_id,
        .neighbor = neighbor_id,
        .similarity = similarity_estimate};
    }
  }};
  par_for(0, num_vertices, [&](const size_t vertex_id) {
    graph->get_vertex(vertex_id).mapOutNghWithIndex(
        vertex_id, compute_similarity);
  });
  return internal::BidirectionalSimilarities(graph->m, undirected_similarities);
}

}  // namespace internal

template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> CosineSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph) const {
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

template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> ApproxCosineSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph) const {
  const size_t exact_threshold{static_cast<size_t>(1.5 * num_samples_)};
  return internal::ApproxCosineEdgeSimilarities(
      graph, num_samples_, exact_threshold, random_seed_);
}

template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> ApproxJaccardSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph) const {
  const size_t exact_threshold{static_cast<size_t>(2.0 * num_samples_)};
  return internal::ApproxJaccardEdgeSimilarities(
      graph, num_samples_, exact_threshold, random_seed_);
}

}  // namespace scan

} // namespace gbbs
