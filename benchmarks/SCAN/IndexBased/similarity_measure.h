// Similarity measures for determining the similarity of adjacent pairs of
// vertices.
#pragma once

#include <array>
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

  template <class Graph>
  pbbs::sequence<EdgeSimilarity>
  AllEdges(Graph* graph) const;
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
// Jaccard similarity with absolute error up to `b`. In practice, setting
// num_samples so high is probably excessive.
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
template <class Graph>
auto DirectGraphByDegree(Graph* graph) {
  uintE* vertex_degree_ranking{rankNodes(*graph, graph->n)};
  const auto filter_predicate{[&](const uintE u, const uintE v, typename Graph::weight_type) {
    return vertex_degree_ranking[u] < vertex_degree_ranking[v];
  }};
  auto directed_graph{filterGraph(*graph, filter_predicate)};
  pbbs::free_array(vertex_degree_ranking);
  return directed_graph;
}

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
template <class Graph, class F>
pbbs::sequence<EdgeSimilarity> AllEdgeNeighborhoodSimilarities(
    Graph* graph,
    F&& neighborhood_sizes_to_similarity) {
  exit(-1);
}

// Implementation of ApproxCosineSimilarities::AllEdges.
//
// `degree_threshold` is a threshold so that we only approximate the similarity
// score between two vertices if their degrees are high enough. (When the
// degrees are low, it's cheap to compute the similarity exactly.)
template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> ApproxCosineEdgeSimilarities(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph,
    const uint32_t num_samples,
    const size_t degree_threshold,
    const size_t random_seed) {
  exit(-1);
}

// Implementation of ApproxJaccardSimilarities::AllEdges.
//
// `degree_threshold` is a threshold so that we only approximate the similarity
// score between two vertices if their degrees are high enough. (When the
// degrees are low, it's cheap to compute the similarity exactly.)
template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> ApproxJaccardEdgeSimilarities(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph,
    const uint32_t num_samples,
    const size_t degree_threshold,
    const size_t random_seed) {
  exit(-1);
}

}  // namespace internal

template <class Graph>
pbbs::sequence<EdgeSimilarity> CosineSimilarity::AllEdges(
    Graph* graph) const {
  using W = typename Graph::weight_type;
  auto directed_graph{internal::DirectGraphByDegree(graph)};
  // Each counter in `counters` holds the number of shared neighbors in `graph`
  // between u and v for some edge {u, v}.
  pbbs::sequence<std::atomic<long long>> counters(
      directed_graph.m, [](size_t) { return std::atomic<long long>{0}; });
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
  par_for(0, graph->n, [&](const size_t vertex_id) {
    auto vertex{directed_graph.get_vertex(vertex_id)};
    const uintT vertex_counter_offset{counter_offsets[vertex_id]};
    const auto intersect{[&](
        const uintE v_id,
        const uintE neighbor_id,
        W w,
        const uintE v_to_neighbor_index) {
      auto neighbor{directed_graph.get_vertex(neighbor_id)};
      const uintT neighbor_counter_offset{counter_offsets[neighbor_id]};
      const auto update_counters{[&](
          const uintE shared_neighbor,
          const uintE vertex_to_shared_index,
          const uintE neighbor_to_shared_index,
          W w1,
          W w2) {
        counters[vertex_counter_offset + vertex_to_shared_index] += static_cast<long long>(1000 * 1000 * w * w2);
        counters[neighbor_counter_offset + neighbor_to_shared_index] += static_cast<long long>(1000 * 1000 * w * w1);
      }};
      counters[vertex_counter_offset + v_to_neighbor_index] +=
        static_cast<long long>(1000 * 1000 *
          internal::intersect_f_with_index_par(
              &vertex, &neighbor, update_counters));
    }};
    constexpr bool kParallel{false};
    vertex.mapOutNghWithIndex(vertex_id, intersect, kParallel);
  });

  pbbs::sequence<double> weighted_deg{
    graph->n,
    [&](const size_t i) {
      double wt = 1000 * 1000;
      auto v = graph->get_vertex(i);
      size_t d = v.getOutDegree();
      auto* nghs = v.getOutNeighbors();
      for (size_t j = 0; j < d; j++) {
        double t = std::get<1>(nghs[j]);
        wt += 1000 * 1000 * t * t;
      }
      return std::sqrt(wt);
    }
  };
  pbbs::sequence<EdgeSimilarity> similarities(graph->m);
  // Convert shared neighbor counts into similarities for each edge.
  par_for(0, directed_graph.n, [&](const size_t vertex_id) {
    const uintT v_counter_offset{counter_offsets[vertex_id]};
    const auto compute_similarity{[&](
        const uintE v_id,
        const uintE u_id,
        W w,
        const uintE v_to_u_index) {
      const uintT counter_index{v_counter_offset + v_to_u_index};
      const double shared_weight{counters[counter_index] + 2 * 1000 * 1000 * w};  // add self edge
      const float similarity{static_cast<float>(shared_weight / (weighted_deg[v_id] * weighted_deg[u_id]))};
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

template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> JaccardSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph) const {
  exit(-1);
}

template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> ApproxCosineSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph) const {
  exit(-1);
}

template <template <typename> class VertexTemplate>
pbbs::sequence<EdgeSimilarity> ApproxJaccardSimilarity::AllEdges(
    symmetric_graph<VertexTemplate, pbbs::empty>* graph) const {
  exit(-1);
}

}  // namespace scan

}  // namespace gbbs
