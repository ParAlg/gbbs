#define NOTMAIN

#include "benchmarks/SCAN/scan.h"

#include <cmath>
#include <atomic>
#include <unordered_set>
#include <utility>
#include <vector>

#include "ligra/bridge.h"
#include "ligra/macros.h"
#include "ligra/vertex.h"
#include "pbbslib/parallel.h"
#include "pbbslib/sample_sort.h"

namespace scan {

namespace {

// Holds a vertex and its degree.
struct VertexDegree {
  uintE vertex_id;
  uintE degree;
};

}  // namespace

namespace internal {

bool operator==(const NeighborSimilarity& a, const NeighborSimilarity& b) {
  return
    std::tie(a.neighbor, a.similarity) == std::tie(b.neighbor, b.similarity);
}

std::ostream&
operator<<(std::ostream& os, const NeighborSimilarity& neighbor_similarity) {
  os << "{neighbor=" << neighbor_similarity.neighbor
     << ", similarity=" << neighbor_similarity.similarity << '}';
  return os;
}

bool operator==(const CoreThreshold& a, const CoreThreshold& b) {
  return
    std::tie(a.vertex_id, a.threshold) == std::tie(b.vertex_id, b.threshold);
}

std::ostream&
operator<<(std::ostream& os, const CoreThreshold& core_threshold) {
  os << "{vertex=" << core_threshold.vertex_id
     << ", threshold=" << core_threshold.threshold << '}';
  return os;
}

// Compute structural similarities (as defined by SCAN) between each pair of
// adjacent vertices.
//
// The structural similarity between two vertices u and v is
//   (size of intersection of closed neighborhoods of u and v) /
//   (geometric mean of size of closed neighborhoods of u and of v)
// where the closed neighborhood of a vertex x consists of all neighbors of x
// along with x itself.
template <class Graph>
StructuralSimilarities ComputeStructuralSimilarities(Graph* graph) {
  using Vertex = typename Graph::vertex;
  using Weight = typename Graph::weight_type;

  StructuralSimilarities similarities{
    graph->m,
    std::make_pair(UndirectedEdge{UINT_E_MAX, UINT_E_MAX}, 0.0),
    std::hash<UndirectedEdge>{}};

  std::vector<sparse_table<
    uintE, pbbslib::empty, std::function<decltype(pbbslib::hash64_2)>>>
    adjacency_list{graph->n};
  parallel_for(0, graph->n, [&graph, &adjacency_list](const size_t vertex_id) {
    Vertex vertex{graph->get_vertex(vertex_id)};
    auto* neighbors{&adjacency_list[vertex_id]};
    *neighbors = make_sparse_table<
      uintE, pbbslib::empty, std::function<decltype(pbbslib::hash64_2)>>(
        // Adding 1 avoids having small tables completely full
        vertex.getOutDegree() + 1,
        {UINT_E_MAX, pbbslib::empty{}},
        pbbslib::hash64_2);

    const auto update_adjacency_list{[&neighbors](
        const uintE source_vertex,
        const uintE neighbor_vertex,
        const Weight weight) {
      neighbors->insert(std::make_pair(neighbor_vertex, pbbslib::empty{}));
    }};
    vertex.mapOutNgh(vertex_id, update_adjacency_list);
  });

  graph->map_edges([&graph, &adjacency_list, &similarities](
        const uintE u_id,
        const uintE v_id,
        const Weight) {
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

        std::atomic<uintE> num_shared_neighbors{0};
        const auto count_shared_neighbors{
          [&larger_degree_vertex_neighbors, &num_shared_neighbors]
          (const uintE, const uintE neighbor, const Weight) {
            if (larger_degree_vertex_neighbors.contains(neighbor)) {
                num_shared_neighbors++;
            }
          }};
        smaller_degree_vertex->mapOutNgh(
            smaller_degree_vertex_id,
            count_shared_neighbors);

        // The neighborhoods we've computed are open neighborhoods -- since
        // structural similarity uses closed neighborhoods, we need to adjust
        // the number and denominator a little.
        similarities.insert({UndirectedEdge{u_id, v_id},
            (num_shared_neighbors + 2) /
                (sqrt(graph->get_vertex(u_id).getOutDegree() + 1) *
                 sqrt(graph->get_vertex(v_id).getOutDegree() + 1))});
      }
  });

  return similarities;
}

template
StructuralSimilarities ComputeStructuralSimilarities(
    symmetric_graph<symmetric_vertex, pbbslib::empty>*);

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
template <class Graph>
NeighborOrder
ComputeNeighborOrder(Graph* graph, const StructuralSimilarities& similarities) {
  using Vertex = typename Graph::vertex;

  NeighborOrder neighbor_order{
    graph->n,
    [&graph](size_t i) {
      return pbbs::sequence<NeighborSimilarity>{
        graph->get_vertex(i).getOutDegree()};
    }
  };

  par_for(0, graph->n, [&](const uintE v) {
    Vertex vertex{graph->get_vertex(v)};
    auto& v_order{neighbor_order[v]};

    par_for(0, vertex.getOutDegree(), [&](const size_t i) {
      const uintE neighbor{vertex.getOutNeighbor(i)};
      const float kNotFound{-1.0};
      const float similarity{
        similarities.find(UndirectedEdge{v, neighbor}, kNotFound)};
      v_order[i] = NeighborSimilarity{
          .neighbor = neighbor, .similarity = similarity};
    });

    // Sort by descending structural similarity
    const auto compare_similarities_descending{[](
        const NeighborSimilarity& a, const NeighborSimilarity& b) {
      return a.similarity > b.similarity;
    }};
    pbbs::sample_sort_inplace(v_order.slice(), compare_similarities_descending);
  });

  return neighbor_order;
}

template
NeighborOrder ComputeNeighborOrder(
    symmetric_graph<symmetric_vertex, pbbslib::empty>*,
    const StructuralSimilarities&);

// Returns a sequence CO where CO[i] for i >= 2 is a list of vertices that
// can be a core when the SCAN parameter mu is set to i. The vertices in
// CO[i], are sorted by their core threshold values, the maximum value of SCAN
// parameter epsilon such that the vertex is a core when mu == i.
//
// CO[0] and CO[1] are left empty --- when mu is less than 2, all vertices are
// always cores and have a core threshold of 1.
CoreOrder ComputeCoreOrder(const NeighborOrder& neighbor_order) {
  if (neighbor_order.empty()) {
    return CoreOrder{};
  }

  pbbs::sequence<VertexDegree> vertex_degrees{
    pbbs::map_with_index<VertexDegree>(
        neighbor_order,
        [](size_t i, const pbbs::sequence<NeighborSimilarity>& similarity) {
          return VertexDegree{
            .vertex_id = static_cast<uintE>(i),
            .degree = static_cast<uintE>(similarity.size())};
        })
  };
  // Sort `vertex_degrees` by ascending degree.
  integer_sort_inplace(
      vertex_degrees.slice(),
      [](const VertexDegree& vertex_degree) { return vertex_degree.degree; });
  const size_t max_degree{vertex_degrees[vertex_degrees.size() - 1].degree};

  // `bucket_offsets[i]` is the first index `j` at which
  // `vertex_degrees[j].degree >= i`.
  pbbs::sequence<uintE> degree_offsets{max_degree + 1};
  const size_t min_degree{vertex_degrees[0].degree};
  par_for(0, min_degree + 1, [&](const size_t j) {
    degree_offsets[j] = 0;
  });
  par_for(1, vertex_degrees.size(), [&](const size_t i) {
    const size_t degree{vertex_degrees[i].degree};
    const size_t prev_degree{vertex_degrees[i-1].degree};
    if (degree != prev_degree) {
      par_for(prev_degree + 1, degree + 1, [&](const size_t j) {
        degree_offsets[j] = i;
      });
    }
  });

  const auto get_core_order{[&](const size_t mu) {
    if (mu <= 1) {
      return pbbs::sequence<CoreThreshold>{};
    }
    // Only vertices with high enough degree can be cores.
    const pbbs::sequence<VertexDegree>& core_vertices{
      vertex_degrees.slice(degree_offsets[mu - 1], vertex_degrees.size())};

    pbbs::sequence<CoreThreshold> core_thresholds{
      pbbs::map<CoreThreshold>(
        core_vertices,
        [&](const VertexDegree& vertex_degree) {
          return CoreThreshold{
            .vertex_id = vertex_degree.vertex_id,
            .threshold =
                neighbor_order[vertex_degree.vertex_id][mu - 2].similarity};
        })};
    // Sort by descending threshold
    const auto compare_threshold_descending{[](
        const CoreThreshold& a, const CoreThreshold& b) {
      return a.threshold > b.threshold;
    }};
    pbbs::sample_sort_inplace(
        core_thresholds.slice(), compare_threshold_descending);
    return core_thresholds;
  }};

  return CoreOrder{max_degree + 2, get_core_order};
}

}  // namespace internal

template <class Graph>
ScanIndex::ScanIndex(Graph* graph) {}

template
ScanIndex::ScanIndex(symmetric_graph<symmetric_vertex, pbbslib::empty>*);

}  // namespace scan
