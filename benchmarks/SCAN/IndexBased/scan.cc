#include "benchmarks/SCAN/IndexBased/scan.h"

#include <algorithm>
#include <tuple>
#include <utility>

#include "benchmarks/Connectivity/UnionFind/union_find_rules.h"
#include "gbbs/bridge.h"
#include "gbbs/helpers/sparse_table.h"

namespace gbbs {
namespace indexed_scan {

namespace {

using DirectedEdge = std::pair<uintE, uintE>;
using VertexSet =
    gbbs::sparse_table<uintE, gbbs::empty, decltype(&parlay::hash64_2)>;

// Creates a `VertexSet` for holding up to `capacity` elements.
VertexSet MakeVertexSet(const size_t capacity) {
  return gbbs::make_sparse_table<uintE, gbbs::empty,
                                 decltype(&parlay::hash64_2)>(
      capacity, {UINT_E_MAX, gbbs::empty{}}, parlay::hash64_2);
}

// Identifies the clusters for the core vertices and populates `clustering`
// accordingly. Non-core vertices in `clustering` are left unchanged.
//
// Arguments:
//   neighbor_order
//     Neighbor order for the graph.
//   cores
//     List of core vertices in the graph.
//   core_similar_edge_counts:
//     `core_similar_edge_counts[i]` is the number of epsilon-similar edges
//     incident on vertex `core[i]`.
//   clustering: Output. Will be updated with the cluster IDs only for the core
//     vertices.
void ClusterCores(const internal::NeighborOrder& neighbor_order,
                  const sequence<uintE>& cores,
                  const sequence<uintE>& core_similar_edge_counts,
                  Clustering* clustering) {
  timer function_timer{"Cluster cores time"};

  VertexSet cores_set{MakeVertexSet(cores.size())};
  parallel_for(0, cores.size(), [&](const size_t i) {
    cores_set.insert(std::make_pair(cores[i], gbbs::empty{}));
  });

  // Get connected components induced by sufficiently similar core-to-core
  // edges.
  parallel_for(0, cores.size(), [&](const size_t i) {
    const uintE core{cores[i]};
    (*clustering)[core] = core;
  });
  constexpr auto find{find_variants::find_compress};
  auto unite{unite_variants::Unite<decltype(find)>{find}};
  parallel_for(0, cores.size(), [&](const size_t i) {
    const uintE core{cores[i]};
    const auto& neighbors{neighbor_order[core]};
    constexpr bool kParallelizeInnerLoop{false};
    parallel_for(0, core_similar_edge_counts[i],
                 [&](const size_t j) {
                   const uintE neighbor{neighbors[j].neighbor};
                   if (core > neighbor && cores_set.contains(neighbor)) {
                     unite(core, neighbor, *clustering);
                   }
                 },
                 kParallelizeInnerLoop);
  });
  parallel_for(0, cores.size(), [&](const size_t i) {
    const uintE core{cores[i]};
    (*clustering)[core] = find(core, *clustering);
  });

  internal::ReportTime(function_timer);
}

// Given `clustering` that's populated for all core vertices, populate
// `clustering` for all non-core vertices.
//
// Arguments:
//   neighbor_order
//     Neighbor order for the graph.
//   cores
//     List of core vertices in the graph.
//   core_similar_edge_counts:
//     `core_similar_edge_counts[i]` is the number of epsilon-similar edges
//     incident on vertex `core[i]`.
//   clustering: Output. When function is called, core vertices must be already
//     be assigned a cluster ID. After function is called, non-core vertices
//     that belong to a cluster will be marked accordingly.
void AttachNoncoresToClusters(const internal::NeighborOrder& neighbor_order,
                              const sequence<uintE>& cores,
                              const sequence<uintE>& core_similar_edge_counts,
                              Clustering* clustering) {
  timer function_timer{"Attach non-cores to clusters time"};
  // Attach each non-core to the same cluster as an arbitrary adjacent core.
  //
  // We could do something smarter like attach each non-core to its most
  // structurally similar adjacent core, but let's keep things simple for now.
  parallel_for(0, cores.size(), [&](const size_t i) {
    const uintE core{cores[i]};
    const uintE core_cluster{(*clustering)[core]};
    const auto& neighbors{neighbor_order[core]};
    constexpr bool kParallelizeInnerLoop{false};
    parallel_for(0, core_similar_edge_counts[i],
                 [&](const size_t j) {
                   const uintE neighbor{neighbors[j].neighbor};
                   auto* neighbor_cluster_address{&(*clustering)[neighbor]};
                   if (*neighbor_cluster_address == kUnclustered) {
                     gbbs::atomic_compare_and_swap(neighbor_cluster_address,
                                                   kUnclustered, core_cluster);
                   }
                 },
                 kParallelizeInnerLoop);
  });
  internal::ReportTime(function_timer);
}

// Same as AttachNoncoresToClusters, but with a consistent, deterministic
// result. Does not attempt to be particularly efficient.
void AttachNoncoresToClustersDeterministic(
    const internal::NeighborOrder& neighbor_order, const sequence<uintE>& cores,
    const sequence<uintE>& core_similar_edge_counts, Clustering* clustering) {
  timer function_timer{"Attach non-cores to clusters time"};
  const size_t num_vertices{clustering->size()};
  // Attach each non-core to its most structurally similar adjacent core,
  // breaking ties in favor of the lowest vertex ID. (We break ties by vertex ID
  // instead of cluster ID so that we don't have to assume that cluster IDs are
  // deterministic.)

  // for a non-core vertex v, `tentative_attachments[v]` represents (similarity
  // score, vertex ID) of the most similar adjacent core to v.
  auto tentative_attachments = sequence<std::pair<float, uintE>>::from_function(
      num_vertices,
      [](const size_t i) { return std::make_pair(-1, UINT_E_MAX); });

  parallel_for(0, cores.size(), [&](const uintE i) {
    const uintE core{cores[i]};
    const auto& neighbors{neighbor_order[core]};
    constexpr bool kParallelizeInnerLoop{false};
    parallel_for(
        0, core_similar_edge_counts[i],
        [&](const size_t j) {
          const uintE neighbor{neighbors[j].neighbor};
          if ((*clustering)[neighbor] == kUnclustered) {
            const float similarity{neighbors[j].similarity};
            while (true) {
              const std::pair<float, uintE> current_attachment{
                  tentative_attachments[neighbor]};
              if (similarity > current_attachment.first ||
                  (similarity == current_attachment.first &&
                   core < current_attachment.second)) {
                if (gbbs::atomic_compare_and_swap(
                        &(tentative_attachments[neighbor]), current_attachment,
                        std::make_pair(similarity, core))) {
                  break;  // successful attachment
                }
              } else {
                break;  // existing attachment is better
              }
            }
          }
        },
        kParallelizeInnerLoop);
  });
  parallel_for(0, num_vertices, [&](const size_t vertex_id) {
    if (tentative_attachments[vertex_id].first > 0.0) {
      (*clustering)[vertex_id] =
          (*clustering)[tentative_attachments[vertex_id].second];
    }
  });
  internal::ReportTime(function_timer);
}

}  // namespace

Index::Index()
    : num_vertices_{0U}, neighbor_order_{}, core_order_{neighbor_order_} {}

Clustering Index::Cluster(const uint64_t mu, const float epsilon,
                          const bool get_deterministic_result) const {
  timer preprocessing_timer{"Cluster - additional preprocessing time"};
  const sequence<uintE> cores{core_order_.GetCores(mu, epsilon)};
  if (cores.empty()) {
    // Nothing is a core. There are no clusters, and every vertex is an outlier.
    return Clustering(num_vertices_, kUnclustered);
  }

  sequence<uintE> core_similar_edge_counts{
      parlay::map<uintE>(cores, [&](const uintE vertex) {
        // Get the number of neighbors of `vertex` that have at least
        // `epsilon` structural similarity with the vertex.
        return parlay::binary_search(
            neighbor_order_[vertex],
            [epsilon](const internal::EdgeSimilarity& es) {
              return es.similarity >= epsilon;
            });
      })};
  Clustering clustering(num_vertices_, kUnclustered);
  internal::ReportTime(preprocessing_timer);

  ClusterCores(neighbor_order_, cores, core_similar_edge_counts, &clustering);
  if (get_deterministic_result) {
    AttachNoncoresToClustersDeterministic(
        neighbor_order_, cores, core_similar_edge_counts, &clustering);
  } else {
    AttachNoncoresToClusters(neighbor_order_, cores, core_similar_edge_counts,
                             &clustering);
  }
  return clustering;
}

void Index::Cluster(const uint64_t mu, const sequence<float>& epsilons,
                    const std::function<void(Clustering&&, size_t)> f,
                    const bool get_deterministic_result) const {
  // TODO(tomtseng): please refactor this. this is messy, copy-and-pasted code
  // written in a rush

  sequence<size_t> sorted_epsilon_indices = sequence<size_t>::from_function(
      epsilons.size(), [](const size_t i) { return i; });
  // Sort epsilons in decreasing order --- as epsilon decreases, more
  // core-to-core edges appear.
  parlay::sample_sort_inplace(make_slice(sorted_epsilon_indices),
                              [&](const size_t i, const size_t j) {
                                return epsilons[i] > epsilons[j];
                              });

  sequence<uintE> previous_cores{};
  sequence<uintE> previous_core_similar_edge_counts{};
  VertexSet cores_set{MakeVertexSet(num_vertices_)};
  Clustering previous_core_clustering = sequence<uintE>::from_function(
      num_vertices_, [&](const size_t i) { return kUnclustered; });
  for (size_t i{0}; i < epsilons.size(); i++) {
    const float epsilon{epsilons[sorted_epsilon_indices[i]]};
    Clustering clustering{std::move(previous_core_clustering)};

    sequence<uintE> cores{core_order_.GetCores(mu, epsilon)};
    sequence<uintE> core_similar_edge_counts{
        parlay::map<uintE>(cores, [&](const uintE vertex) {
          // Get the number of neighbors of `vertex` that have at least
          // `epsilon` structural similarity with the vertex.
          return parlay::binary_search(
              neighbor_order_[vertex],
              [epsilon](const internal::EdgeSimilarity& es) {
                return es.similarity >= epsilon;
              });
        })};

    // ClusterCores() logic -- get connected components induced by sufficiently
    // similar core-to-core edges
    parallel_for(previous_cores.size(), cores.size(), [&](const size_t j) {
      const uintE core{cores[j]};
      cores_set.insert(std::make_pair(core, gbbs::empty{}));
      clustering[core] = core;
    });
    constexpr auto find{find_variants::find_compress};
    auto unite{unite_variants::Unite<decltype(find)>{find}};
    parallel_for(0, cores.size(), [&](const size_t j) {
      const uintE core{cores[j]};
      const auto& neighbors{neighbor_order_[core]};
      constexpr bool kParallelizeInnerLoop{false};
      // only operate on new edges for this iteration, hence the ternary
      // statement
      parallel_for(
          j < previous_cores.size() ? previous_core_similar_edge_counts[j] : 0,
          core_similar_edge_counts[j],
          [&](const size_t k) {
            const uintE neighbor{neighbors[k].neighbor};
            // the `core > neighbor` check tries to avoid the extra work of
            // adding
            // both directions of an edge. but sometimes when a new core is
            // added we
            // do need to check the other direction of an edge, hence the `j >=
            // previous_cores.size()` check
            if ((j >= previous_cores.size() || core > neighbor) &&
                cores_set.contains(neighbor)) {
              unite(core, neighbor, clustering);
            }
          },
          kParallelizeInnerLoop);
    });
    parallel_for(0, cores.size(), [&](const size_t j) {
      const uintE core{cores[j]};
      clustering[core] = find(core, clustering);
    });

    previous_core_clustering = clustering;

    if (get_deterministic_result) {
      AttachNoncoresToClustersDeterministic(
          neighbor_order_, cores, core_similar_edge_counts, &clustering);
    } else {
      AttachNoncoresToClusters(neighbor_order_, cores, core_similar_edge_counts,
                               &clustering);
    }

    f(std::move(clustering), sorted_epsilon_indices[i]);

    previous_cores = std::move(cores);
    previous_core_similar_edge_counts = std::move(core_similar_edge_counts);
  }
}

}  // namespace indexed_scan
}  // namespace gbbs
