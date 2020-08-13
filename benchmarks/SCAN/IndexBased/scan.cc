#include "benchmarks/SCAN/IndexBased/scan.h"

#include <algorithm>
#include <tuple>
#include <utility>

#include "benchmarks/Connectivity/UnionFind/union_find_rules.h"
#include "gbbs/bridge.h"
#include "gbbs/pbbslib/sparse_table.h"

namespace gbbs {
namespace indexed_scan {

namespace {

using DirectedEdge = std::pair<uintE, uintE>;
using VertexSet =
  pbbslib::sparse_table<uintE, pbbslib::empty, decltype(&pbbslib::hash64_2)>;

// Creates a `VertexSet` for holding up to `capacity` elements.
VertexSet MakeVertexSet(const size_t capacity) {
  return pbbslib::make_sparse_table<
    uintE, pbbslib::empty, decltype(&pbbslib::hash64_2)>(
      capacity, {UINT_E_MAX, pbbslib::empty{}}, pbbslib::hash64_2);
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
void ClusterCores(
    const internal::NeighborOrder& neighbor_order,
    const pbbs::sequence<uintE>& cores,
    const pbbs::sequence<uintE>& core_similar_edge_counts,
    Clustering* clustering) {
  timer function_timer{"Cluster cores time"};

  VertexSet cores_set{MakeVertexSet(cores.size())};
  par_for(0, cores.size(), [&](const size_t i) {
    cores_set.insert(std::make_pair(cores[i], pbbslib::empty{}));
  });

  // Get connected components induced by sufficiently similar core-to-core
  // edges.
  par_for(0, cores.size(), [&](const size_t i) {
      const uintE core{cores[i]};
      (*clustering)[core] = core;
  });
  constexpr auto find{find_variants::find_compress};
  auto unite{unite_variants::Unite<decltype(find)>{find}};
  par_for(0, cores.size(), [&](const size_t i) {
    const uintE core{cores[i]};
    const auto& neighbors{neighbor_order[core]};
    constexpr bool kParallelizeInnerLoop{false};
    par_for(0, core_similar_edge_counts[i], [&](const size_t j) {
      const uintE neighbor{neighbors[j].neighbor};
      if (core > neighbor && cores_set.contains(neighbor)) {
        unite(core, neighbor, *clustering);
      }
    }, kParallelizeInnerLoop);
  });
  par_for(0, cores.size(), [&](const size_t i) {
      const uintE core{cores[i]};
      (*clustering)[core] = find(core, *clustering);
  });

  cores_set.del();
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
void AttachNoncoresToClusters(
    const internal::NeighborOrder& neighbor_order,
    const pbbs::sequence<uintE>& cores,
    const pbbs::sequence<uintE>& core_similar_edge_counts,
    Clustering* clustering) {
  timer function_timer{"Attach non-cores to clusters time"};
  // Attach each non-core to the same cluster as an arbitrary adjacent core.
  //
  // We could do something smarter like attach each non-core to its most
  // structurally similar adjacent core, but let's keep things simple for now.
  par_for(0, cores.size(), [&](const size_t i) {
    const uintE core{cores[i]};
    const uintE core_cluster{(*clustering)[core]};
    const auto& neighbors{neighbor_order[core]};
    constexpr bool kParallelizeInnerLoop{false};
    par_for(0, core_similar_edge_counts[i], [&](const size_t j) {
      const uintE neighbor{neighbors[j].neighbor};
      auto* neighbor_cluster_address{&(*clustering)[neighbor]};
      if (*neighbor_cluster_address == kUnclustered) {
        pbbs::atomic_compare_and_swap(
            neighbor_cluster_address, kUnclustered, core_cluster);
      }
    }, kParallelizeInnerLoop);
  });
  internal::ReportTime(function_timer);
}

// Same as AttachNoncoresToClusters, but with a consistent, deterministic
// result. Does not attempt to be particularly efficient.
void AttachNoncoresToClustersDeterministic(
    const internal::NeighborOrder& neighbor_order,
    const pbbs::sequence<uintE>& cores,
    const pbbs::sequence<uintE>& core_similar_edge_counts,
    Clustering* clustering) {
  timer function_timer{"Attach non-cores to clusters time"};
  const size_t num_vertices{clustering->size()};
  // Attach each non-core to its most structurally similar adjacent core,
  // breaking ties in favor of the lowest vertex ID. (We break ties by vertex ID
  // instead of cluster ID so that we don't have to assume that cluster IDs are
  // deterministic.)

  // for a non-core vertex v, `tentative_attachments[v]` represents (similarity
  // score, vertex ID) of the most similar adjacent core to v.
  pbbs::sequence<std::pair<float, uintE>> tentative_attachments{
    num_vertices,
    [](const size_t i) { return std::make_pair(-1, UINT_E_MAX); }};

  par_for(0, cores.size(), [&](const uintE i) {
    const uintE core{cores[i]};
    const auto& neighbors{neighbor_order[core]};
    constexpr bool kParallelizeInnerLoop{false};
    par_for(0, core_similar_edge_counts[i], [&](const size_t j) {
      const uintE neighbor{neighbors[j].neighbor};
      if ((*clustering)[neighbor] == kUnclustered) {
        const float similarity{neighbors[j].similarity};
        while (true) {
          const std::pair<float, uintE> current_attachment{
            tentative_attachments[neighbor]};
          if (similarity > current_attachment.first ||
              (similarity == current_attachment.first &&
               core < current_attachment.second)) {
            if (pbbslib::atomic_compare_and_swap(
                &(tentative_attachments[neighbor]),
                current_attachment,
                std::make_pair(similarity, core))) {
              break;  // successful attachment
            }
          } else {
            break;  // existing attachment is better
          }
        }
      }
    }, kParallelizeInnerLoop);
  });
  par_for(0, num_vertices, [&](const size_t vertex_id) {
    if (tentative_attachments[vertex_id].first > 0.0) {
      (*clustering)[vertex_id] =
        (*clustering)[tentative_attachments[vertex_id].second];
    }
  });
  internal::ReportTime(function_timer);
}

}  // namespace

Clustering Index::Cluster(
    const uint64_t mu,
    const float epsilon,
    const bool get_deterministic_result) const {
  const pbbs::sequence<uintE> cores{core_order_.GetCores(mu, epsilon)};
  if (cores.empty()) {
    // Nothing is a core. There are no clusters, and every vertex is an outlier.
    return Clustering(num_vertices_, kUnclustered);
  }

  timer preprocessing_timer{"Cluster - additional preprocessing time"};
  pbbs::sequence<uintE> core_similar_edge_counts{
      pbbs::map<uintE>(
          cores,
          [&](const uintE vertex) {
            // Get the number of neighbors of `vertex` that have at least
            // `epsilon` structural similarity with the vertex.
            return internal::BinarySearch(
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
    AttachNoncoresToClusters(
        neighbor_order_, cores, core_similar_edge_counts, &clustering);
  }
  return clustering;
}

pbbs::sequence<Clustering> Index::Cluster(
    const uint64_t mu,
    const pbbs::sequence<float>& epsilons,
    const bool get_deterministic_result) const {
  // TODO(tomtseng): please refactor this. this is messy, copy-and-pasted code
  // written in a rush

  pbbs::sequence<Clustering> clusterings{
    epsilons.size(), [&](const size_t i) { return Clustering{}; }};
  pbbs::sequence<size_t> sorted_epsilon_indices{
    epsilons.size(), [](const size_t i) { return i; }};
  // Sort epsilons in decreasing order --- as epsilon decreases, more
  // core-to-core edges appear.
  pbbs::sample_sort_inplace(
      sorted_epsilon_indices.slice(),
      [&](const size_t i, const size_t j) {
        return epsilons[i] > epsilons[j];
      });

  pbbs::sequence<uintE> previous_cores{};
  pbbs::sequence<uintE> previous_core_similar_edge_counts{};
  VertexSet cores_set{MakeVertexSet(num_vertices_)};
  Clustering previous_core_clustering{
    num_vertices_,
    [&](const size_t i) { return kUnclustered; }};
  for (size_t i{0}; i < epsilons.size(); i++) {
    const float epsilon{epsilons[sorted_epsilon_indices[i]]};
    Clustering* clustering = &(clusterings[sorted_epsilon_indices[i]]);
    *clustering = previous_core_clustering;

    pbbs::sequence<uintE> cores{core_order_.GetCores(mu, epsilon)};
    pbbs::sequence<uintE> core_similar_edge_counts{
        pbbs::map<uintE>(
            cores,
            [&](const uintE vertex) {
              // Get the number of neighbors of `vertex` that have at least
              // `epsilon` structural similarity with the vertex.
              return internal::BinarySearch(
                  neighbor_order_[vertex],
                  [epsilon](const internal::EdgeSimilarity& es) {
                    return es.similarity >= epsilon;
                  });
            })};

    // ClusterCores() logic -- get connected components induced by sufficiently
    // similar core-to-core edges
    par_for(previous_cores.size(), cores.size(), [&](const size_t j) {
      const uintE core{cores[j]};
      cores_set.insert(std::make_pair(core, pbbslib::empty{}));
      (*clustering)[core] = core;
    });
    constexpr auto find{find_variants::find_compress};
    auto unite{unite_variants::Unite<decltype(find)>{find}};
    par_for(0, cores.size(), [&](const size_t j) {
      const uintE core{cores[j]};
      const auto& neighbors{neighbor_order_[core]};
      constexpr bool kParallelizeInnerLoop{false};
      // only operate on new edges for this iteration, hence the ternary
      // statement
      par_for(
          j < previous_cores.size()
            ? previous_core_similar_edge_counts[j]
            : 0,
          core_similar_edge_counts[j],
          [&](const size_t k) {
        const uintE neighbor{neighbors[k].neighbor};
        // the `core > neighbor` check tries to avoid the extra work of adding
        // both directions of an edge. but sometimes when a new core is added we
        // do need to check the other direction of an edge, hence the `j >=
        // previous_cores.size()` check
        if ((j >= previous_cores.size() || core > neighbor) &&
            cores_set.contains(neighbor)) {
          unite(core, neighbor, *clustering);
        }
      }, kParallelizeInnerLoop);
    });
    par_for(0, cores.size(), [&](const size_t j) {
        const uintE core{cores[j]};
        (*clustering)[core] = find(core, *clustering);
    });

    previous_core_clustering = *clustering;

    if (get_deterministic_result) {
      AttachNoncoresToClustersDeterministic(
          neighbor_order_, cores, core_similar_edge_counts, clustering);
    } else {
      AttachNoncoresToClusters(
          neighbor_order_, cores, core_similar_edge_counts, clustering);
    }

    previous_cores = std::move(cores);
    previous_core_similar_edge_counts = std::move(core_similar_edge_counts);
  }

  cores_set.del();
  return clusterings;
}

}  // namespace indexed_scan
}  // namespace gbbs
