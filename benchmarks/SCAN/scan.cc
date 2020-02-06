#define NOTMAIN

#include "benchmarks/SCAN/scan.h"

#include <cmath>
#include <algorithm>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"
#include "ligra/bridge.h"
#include "ligra/macros.h"
#include "pbbslib/binary_search.h"
#include "pbbslib/parallel.h"
#include "pbbslib/sample_sort.h"
#include "utils/assert.h"

namespace scan {

namespace {

using Weight = pbbslib::empty;
using Vertex = symmetric_vertex<Weight>;

using DirectedEdge = std::pair<uintE, uintE>;
using VertexSet =
  sparse_table<uintE, pbbslib::empty, decltype(&pbbslib::hash64_2)>;

// Holds a vertex and its degree.
struct VertexDegree {
  uintE vertex_id;
  uintE degree;
};

// Creates a `VertexSet` for holding up to `capacity` elements.
VertexSet MakeVertexSet(const size_t capacity) {
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
  if (sequence.size() < hi) {
    hi = sequence.size();
  }

  return pbbs::binary_search(sequence.slice(lo, hi), predicate);
}

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
StructuralSimilarities ComputeStructuralSimilarities(
    symmetric_graph<symmetric_vertex, pbbslib::empty>* graph) {
  StructuralSimilarities similarities{
    graph->m,
    std::make_pair(UndirectedEdge{UINT_E_MAX, UINT_E_MAX}, 0.0),
    std::hash<UndirectedEdge>{}};

  std::vector<VertexSet> adjacency_list{graph->n};
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

  // TODO(tom.tseng): This all might be overkill --- look at
  // ligra/vertex.h intersection::intersect
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

        const auto is_shared_neighbor{
          [&](const uintE, const uintE neighbor, const Weight) {
            return larger_degree_vertex_neighbors.contains(neighbor);
          }};
        const auto add_monoid{pbbslib::addm<size_t>()};
        const size_t num_shared_neighbors{
          smaller_degree_vertex->reduceOutNgh<size_t>(
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
NeighborOrder ComputeNeighborOrder(
    symmetric_graph<symmetric_vertex, pbbslib::empty>* graph,
    const StructuralSimilarities& similarities) {
  NeighborOrder neighbor_order{
    graph->n,
    [&graph](const size_t i) {
      return pbbs::sequence<NeighborSimilarity>{
        graph->get_vertex(i).getOutDegree()};
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
        [](const size_t i,
           const pbbs::sequence<NeighborSimilarity>& similarity) {
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
  pbbs::sequence<uintE> degree_offsets{
    pbbs::sequence<uintE>::no_init(max_degree + 1)};
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
    // Sort by descending threshold.
    const auto compare_threshold_descending{
      [](const CoreThreshold& a, const CoreThreshold& b) {
        return a.threshold > b.threshold;
      }};
    pbbs::sample_sort_inplace(
        core_thresholds.slice(), compare_threshold_descending);
    return core_thresholds;
  }};

  return CoreOrder{max_degree + 2, get_core_order};
}

}  // namespace internal

ScanIndex::ScanIndex(symmetric_graph<symmetric_vertex, pbbslib::empty>* graph)
  : num_vertices_{graph->n}
  , neighbor_order_{
      internal::ComputeNeighborOrder(
          graph,
          internal::ComputeStructuralSimilarities(graph))}
  , core_order_{internal::ComputeCoreOrder(neighbor_order_)} {}

Clustering ScanIndex::Cluster(const float epsilon, const uint64_t mu) const {
  if (mu <= 1) {
    // Every vertex is a core. Return connected components on edges with
    // similarity > epsilon.
    ABORT("SCAN for `mu <= 1` is unimplemented");
  }
  if (mu >= core_order_.size()) {
    // Nothing is a core. There are no clusters, and every vertex is an outlier.
    return Clustering {
        .num_clusters = 0,
        .clusters_by_vertex =
          pbbs::sequence<VertexType>{num_vertices_, VertexType{Outlier{}}}
    };
  }

  const size_t cores_end{
    BinarySearch<internal::CoreThreshold>(
        core_order_[mu],
        [epsilon](const internal::CoreThreshold& core_threshold) {
          return core_threshold.threshold >= epsilon;
        })};
  const pbbs::range cores{core_order_[mu].slice(0, cores_end)};

  pbbs::sequence<size_t> epsilon_neighborhood_offsets{
      pbbs::map<size_t>(
          cores,
          [&](const internal::CoreThreshold& core_threshold) {
            // Get the number of neighbors of core vertex
            // `core_threshold.vertex_id` that have at least `epsilon`
            // structural similarity with the core.
            return BinarySearch<internal::NeighborSimilarity>(
                neighbor_order_[core_threshold.vertex_id],
                [epsilon](const internal::NeighborSimilarity& ns) {
                  return ns.similarity >= epsilon;
                });
          })};
  const size_t num_core_incident_edges{
    pbbslib::scan_add_inplace(epsilon_neighborhood_offsets)};
  // List of edges with structural similarity at least `epsilon` that are
  // incident on a core.
  // Note: edges of the form (core vertex, non-core) will appear only once in
  // this list. On the other hand, edges of the form (core vertex 1, core
  // vertex 2) will appear twice in this list, once in each direction.
  const pbbs::sequence<DirectedEdge> core_incident_edges{
    pbbs::sequence<DirectedEdge>::no_init(num_core_incident_edges)};
  par_for(0, cores.size(), [&](const size_t i) {
    const size_t offset{epsilon_neighborhood_offsets[i]};
    const size_t size{
      (i == cores.size()
       ? num_core_incident_edges
       : epsilon_neighborhood_offsets[i + 1]) - offset};
    const uintE core_id{cores[i].vertex_id};
    const auto& core_neighbors{neighbor_order_[core_id]};
    par_for(0, size, [&](const size_t j) {
      core_incident_edges[offset + j] =
        std::make_pair(core_id, core_neighbors[j].neighbor);
    });
  });

  VertexSet cores_set{MakeVertexSet(cores.size())};
  par_for(0, cores.size(), [&](const size_t i) {
    cores_set.insert(std::make_pair(cores[i].vertex_id, pbbslib::empty{}));
  });

  // `partitioned_core_edges` is `core_incident_edges` partitioned into edges
  // whose endpoints are both cores and edges that have a non-core endpoint.
  pbbs::sequence<DirectedEdge> partitioned_core_edges{};
  size_t num_core_to_core_edges{0};
  std::tie(partitioned_core_edges, num_core_to_core_edges) =
    pbbs::split_two_with_predicate(
        core_incident_edges,
        [&cores_set](const DirectedEdge& edge) {
          // Only need to check second endpoint. First endpoint is a core.
          return !cores_set.contains(edge.second);
        });
  auto core_to_core_edges{
    partitioned_core_edges.slice(0, num_core_to_core_edges)};

  // Create graph consisting of edges of sufficient similarity between cores.
  // The vertex ids of the cores are kept the same for simplicity, so actually
  // all the non-core vertices are also in the graph as singletons.
  symmetric_graph<symmetric_vertex, pbbslib::empty> core_graph{
    sym_graph_from_edges<pbbslib::empty>(
        core_to_core_edges,
        num_vertices_,
        [](const DirectedEdge& edge) { return edge.first; },
        [](const DirectedEdge& edge) { return edge.second; },
        [](const DirectedEdge& edge) { return pbbslib::empty{}; })};

  Clustering clustering{
    .num_clusters = 0,
    // Mark everything as an outlier for now.
    .clusters_by_vertex =
      pbbs::sequence<VertexType>{num_vertices_, VertexType{Outlier{}}}
  };

  // Get connected components in the resulting core graph. Relabel the resulting
  // component IDs to be contiguous and to ignore the non-core vertices. This
  // identifies the clusters for all the core vertices.
  pbbs::sequence<parent> core_connected_components{
    workefficient_cc::CC(core_graph)};
  pbbs::sequence<uintE> component_relabel_map{num_vertices_, 0U};
  par_for(0, cores.size(), [&](const size_t i) {
    const uintE cluster{core_connected_components[cores[i].vertex_id]};
    if (component_relabel_map[cluster] == 0) {
      component_relabel_map[cluster] = 1;
    }
  });
  clustering.num_clusters = pbbslib::scan_add_inplace(component_relabel_map);
  par_for(0, cores.size(), [&](const size_t i) {
    const uintE core{cores[i].vertex_id};
    clustering.clusters_by_vertex[core] =
      ClusterMember{
        .clusters =
          pbbs::sequence{
            1,
            component_relabel_map[core_connected_components[core]]}
      };
  });

  // Non-core vertices that are epsilon-similar to cores are assigned to the
  // core's cluster. These vertices may belong to more than one cluster.
  auto core_to_noncore_edges{
    partitioned_core_edges.slice(num_core_to_core_edges,
        partitioned_core_edges.size())};
  par_for(0, core_to_noncore_edges.size(), [&](const size_t i) {
    // Replace core vertex ID with its cluster ID.
    core_to_noncore_edges[i].first =
      std::get<ClusterMember>(
          clustering.clusters_by_vertex[core_to_noncore_edges[i].first])
      .clusters[0];
  });
  pbbs::sample_sort_inplace(
      core_to_noncore_edges,
      // Sort first by non-core endpoint.
      [](const DirectedEdge& a, const DirectedEdge& b) {
        return std::tie(a.second, a.first) < std::tie(b.second, b.first);
      });
  // Table storing the last index at which a non-core vertex in
  // core_to_noncore_edges[*].second appears in core_to_noncore_edges.
  auto noncore_ends{
    make_sparse_table<uintE, uintT, decltype(&pbbslib::hash64_2)>(
      std::min<size_t>(num_vertices_, core_to_noncore_edges.size()) + 1,
      {UINT_E_MAX, UINT_T_MAX},
      pbbslib::hash64_2)};
  par_for(0, core_to_noncore_edges.size(), [&](const size_t i) {
    const uintE noncore{core_to_noncore_edges[i].second};
    if (i == core_to_noncore_edges.size() - 1 ||
        noncore != core_to_noncore_edges[i + 1].second) {
      noncore_ends.insert({noncore, i + 1});
    }
  });
  par_for(0, core_to_noncore_edges.size(), [&](const size_t i) {
    const uintE noncore{core_to_noncore_edges[i].second};
    if (i == 0 || noncore != core_to_noncore_edges[i - 1].second) {
      constexpr uintT kDefaultEnd{UINT_T_MAX};
      const uintT noncore_end{noncore_ends.find(noncore, kDefaultEnd)};
      VertexSet clusters{
        MakeVertexSet(
            std::min<size_t>(noncore_end - i, clustering.num_clusters))};
      par_for(i, noncore_end, [&](const size_t j) {
        const uintE cluster{core_to_noncore_edges[i].first};
        if (j == i || cluster != core_to_noncore_edges[i - 1].first) {
          clusters.insert({cluster, pbbslib::empty{}});
        }
      });
      clustering.clusters_by_vertex[noncore] =
        ClusterMember{
          .clusters =
            pbbs::map<uintE>(
                clusters.entries(),
                [](const std::tuple<uintE, pbbslib::empty> kv) {
                  return std::get<0>(kv);
                })
        };
    }
  });

  par_for(0, num_vertices_, [&](const size_t i) {
    auto* const vertex_type{&clustering.clusters_by_vertex[i]};
    if (!std::holds_alternative<ClusterMember>(*vertex_type)) {
      // Determine whether remaining vertex i is a hub or outlier.

      const auto& neighbors{neighbor_order_[i]};
      bool is_hub{false};
      // `candidate_cluster` holds a cluster ID that vertex i is adjacent to, or
      // UINT_E_MAX before such a cluster is found.
      uintE candidate_cluster{UINT_E_MAX};
      par_for(0, neighbors.size(), [&](const size_t j) {
        const ClusterMember* const neighbor_clusters{
          std::get_if<ClusterMember>(
              &clustering.clusters_by_vertex[neighbors[j].neighbor])};
        if (neighbor_clusters != nullptr) {
          if (neighbor_clusters->clusters.size() > 1) {
            if (!is_hub) {
              is_hub = true;
            }
          } else {
            const uintE neighbor_cluster{neighbor_clusters->clusters[0]};
            // If `candidate_cluster` is at its default value of UINT_E_MAX,
            // assign `neighbor_cluster` to it. Otherwise, if it has a value
            // differing from `neighbor_cluster`, then we know vertex i is
            // adjacent to multiple clusters and is a hub.
            if (!(candidate_cluster == UINT_E_MAX
                  && pbbs::atomic_compare_and_swap(
                        &candidate_cluster, UINT_E_MAX, neighbor_cluster))
                && candidate_cluster != neighbor_cluster
                && !is_hub) {
              is_hub = true;
            }
          }
        }

        *vertex_type = is_hub? VertexType{Hub{}} : VertexType{Outlier{}};
      });
    }
  });

  return clustering;
}

}  // namespace scan
