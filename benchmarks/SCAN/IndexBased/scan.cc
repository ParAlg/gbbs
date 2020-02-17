#define NOTMAIN

#include "benchmarks/SCAN/scan.h"

#include <algorithm>
#include <tuple>
#include <utility>

#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"
#include "ligra/bridge.h"
#include "ligra/macros.h"
#include "pbbslib/binary_search.h"
#include "pbbslib/parallel.h"
#include "pbbslib/sample_sort.h"
#include "utils/assert.h"

namespace scan {

namespace {

using DirectedEdge = std::pair<uintE, uintE>;

// Get edges with structural similarity at least `epsilon` that are incident on
// a vertex in `vertices`.
//
// Edges of the form {v, u} where v is in `vertices` but `u` is not will appear
// once in the output in the form (v, u). Edges where both `v` and `u`
// are in `vertices`, on the other hand will appear twice, as (v, u) and as (u,
// v).
pbbs::sequence<DirectedEdge> GetSimilarIncidentEdges(
    const internal::NeighborOrder& neighbor_order,
    const pbbs::sequence<uintE>& vertices,
    const float epsilon) {
  pbbs::sequence<size_t> epsilon_neighborhood_offsets{
      pbbs::map<size_t>(
          vertices,
          [&](const uintE vertex) {
            // Get the number of neighbors of `vertex` that have at least
            // `epsilon` structural similarity with the core.
            return internal::BinarySearch<internal::NeighborSimilarity>(
                neighbor_order[vertex],
                [epsilon](const internal::NeighborSimilarity& ns) {
                  return ns.similarity >= epsilon;
                });
          })};
  const size_t num_incident_edges{
    pbbslib::scan_add_inplace(epsilon_neighborhood_offsets)};
  const pbbs::sequence<DirectedEdge> incident_edges{
    pbbs::sequence<DirectedEdge>::no_init(num_incident_edges)};
  par_for(0, vertices.size(), [&](const size_t i) {
    const uintE vertex{vertices[i]};
    const size_t offset{epsilon_neighborhood_offsets[i]};
    const size_t num_incident_on_vertex{
      (i + 1 == vertices.size()
       ? num_incident_edges
       : epsilon_neighborhood_offsets[i + 1]) - offset};
    const auto& neighbors{neighbor_order[vertex]};
    par_for(0, num_incident_on_vertex, [&](const size_t j) {
      incident_edges[offset + j] =
        std::make_pair(vertex, neighbors[j].neighbor);
    });
  });
  return incident_edges;
}

// Given a list of edges of sufficient similarity between core vertices,
// identifies the clusters for the core vertices and populates
// `clustering->num_clusters` and `clustering->clusters_by_vertex` accordingly.
// Non-core vertices in `clustering.clusters_by_vertex` are left unchanged.
//
// Arguments:
//   num_vertices: Number of vertices in the whole graph.
//   cores: List of IDs of vertices that are SCAN cores.
//   core_to_core_edges: Directed edges between cores of at least epsilon
//     similarity. This list may be modified by this function.
//   clustering: Output clustering, to be updated with the number of clusters as
//     well as the clustering for only the core vertices.
void ClusterCores(
    const uintE num_vertices,
    const pbbs::sequence<uintE>& cores,
    pbbs::range<DirectedEdge*>* core_to_core_edges,
    Clustering* clustering) {
  // Create graph consisting of edges between cores. The vertex ids of the
  // cores are kept the same for simplicity, so actually all the non-core
  // vertices are also in the graph as singletons.
  symmetric_graph<symmetric_vertex, pbbslib::empty> core_graph{
    sym_graph_from_edges<pbbslib::empty>(
        *core_to_core_edges,
        num_vertices,
        [](const DirectedEdge& edge) { return edge.first; },
        [](const DirectedEdge& edge) { return edge.second; },
        [](const DirectedEdge& edge) { return pbbslib::empty{}; })};

  // Get connected components in the resulting core graph. Relabel the resulting
  // component IDs to be contiguous and to ignore the non-core vertices. This
  // identifies the clusters for all the core vertices.
  pbbs::sequence<parent> connected_components{workefficient_cc::CC(core_graph)};
  pbbs::sequence<uintE> component_relabel_map{num_vertices, 0U};
  par_for(0, cores.size(), [&](const size_t i) {
    const uintE cluster{connected_components[cores[i]]};
    if (component_relabel_map[cluster] == 0) {
      component_relabel_map[cluster] = 1;
    }
  });
  clustering->num_clusters = pbbslib::scan_add_inplace(component_relabel_map);
  par_for(0, cores.size(), [&](const size_t i) {
    const uintE core{cores[i]};
    clustering->clusters_by_vertex[core] =
      ClusterMember{
        .clusters =
          pbbs::sequence{1, component_relabel_map[connected_components[core]]}
      };
  });
}

// Given `clustering` where `clustering.clusters_by_vertex` is correct for all
// core vertices and a list of [core -> non-core] edges of sufficient
// similarity, populate `clustering.clusters_by_vertex` for all non-core
// vertices listed in `core_to_noncore_edges`.
//
// Arguments:
// - num_vertices: Number of vertices in the whole graph.
// - core_to_noncore_edges: Directed edges of at least epsilon similarityfrom
//   core to noncore vertices of at least epsilon similarity. This list may be
//   modified by this function.
// - clustering: Output. When function is called, core vertices must be already
//   be marked as ClusterMembers. After function is called, non-core vertices
//   that are ClusterMembers will be identified and marked accordingly.
void AttachNoncoresToClusters(
    const uintE num_vertices,
    pbbs::range<DirectedEdge*>* core_to_noncore_edges,
    Clustering* clustering) {
  par_for(0, core_to_noncore_edges->size(), [&](const size_t i) {
    // Replace core vertex ID with its cluster ID.
    (*core_to_noncore_edges)[i].first =
      std::get<ClusterMember>(
          clustering->clusters_by_vertex[(*core_to_noncore_edges)[i].first])
      .clusters[0];
  });
  pbbs::sample_sort_inplace(
      *core_to_noncore_edges,
      // Sort first by non-core endpoint.
      [](const DirectedEdge& a, const DirectedEdge& b) {
        return std::tie(a.second, a.first) < std::tie(b.second, b.first);
      });
  // Table storing the last index at which a non-core vertex in
  // (*core_to_noncore_edges)[...].second appears in core_to_noncore_edges.
  auto noncore_ends{
    make_sparse_table<uintE, uintT, decltype(&pbbslib::hash64_2)>(
      std::min<size_t>(num_vertices, core_to_noncore_edges->size()) + 1,
      {UINT_E_MAX, UINT_T_MAX},
      pbbslib::hash64_2)};
  par_for(0, core_to_noncore_edges->size(), [&](const size_t i) {
    const uintE noncore{(*core_to_noncore_edges)[i].second};
    if (i == core_to_noncore_edges->size() - 1 ||
        noncore != (*core_to_noncore_edges)[i + 1].second) {
      noncore_ends.insert({noncore, i + 1});
    }
  });
  par_for(0, core_to_noncore_edges->size(), [&](const size_t i) {
    const uintE noncore{(*core_to_noncore_edges)[i].second};
    if (i == 0 || noncore != (*core_to_noncore_edges)[i - 1].second) {
      constexpr uintT kDefaultEnd{UINT_T_MAX};
      const uintT noncore_end{noncore_ends.find(noncore, kDefaultEnd)};
      internal::VertexSet clusters{
      internal::MakeVertexSet(
            std::min<size_t>(noncore_end - i, clustering->num_clusters))};
      par_for(i, noncore_end, [&](const size_t j) {
        const uintE cluster{(*core_to_noncore_edges)[i].first};
        if (j == i || cluster != (*core_to_noncore_edges)[i - 1].first) {
          clusters.insert({cluster, pbbslib::empty{}});
        }
      });
      clustering->clusters_by_vertex[noncore] =
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
}

// Gets the SCAN clusters, populating `clustering`, but does not determine hubs
// or outliers. Hubs and outliers in `clustering->clusters_by_vertex` are left
// unchanged.
//
// Arguments:
// - num_vertices: Number of vertices in the graph.
// - cores: List of IDs for SCAN core vertices.
// - core_similar_incident_edges: List of directed edges incident on cores with
//   similarity at least epsilon, where the first endpoint of each edge is a
//   core.
// - clustering: Output. `clustering->clusters_by_vertex` for vertices that are
//   not members of clusters will be unchanged.
void GetClustersFromCores(
    const uintE num_vertices,
    const pbbs::sequence<uintE>& cores,
    const pbbs::sequence<DirectedEdge>& core_similar_incident_edges,
    Clustering* clustering) {
  internal::VertexSet cores_set{internal::MakeVertexSet(cores.size())};
  par_for(0, cores.size(), [&](const size_t i) {
    cores_set.insert(std::make_pair(cores[i], pbbslib::empty{}));
  });

  // `partitioned_edges` is `core_similar_incident_edges` partitioned into
  // edges whose endpoints are both cores and edges that have a non-core
  // endpoint.
  pbbs::sequence<DirectedEdge> partitioned_edges{};
  size_t num_core_to_core_edges{0};
  std::tie(partitioned_edges, num_core_to_core_edges) =
    pbbs::split_two_with_predicate(
        core_similar_incident_edges,
        [&cores_set](const DirectedEdge& edge) {
          // Only need to check second endpoint. First endpoint is a core.
          return !cores_set.contains(edge.second);
        });

  pbbs::range<DirectedEdge*> core_to_core_edges{
    partitioned_edges.slice(0, num_core_to_core_edges)};
  ClusterCores(num_vertices, cores, &core_to_core_edges, clustering);

  pbbs::range<DirectedEdge*> core_to_noncore_edges{
    partitioned_edges.slice(num_core_to_core_edges, partitioned_edges.size())};
  AttachNoncoresToClusters(num_vertices, &core_to_noncore_edges, clustering);
}

// For every vertex in `clusters_by_vertex` that is not marked as a
// `ClusterMember`, mark it as either a hub or an outlier.
//
// All cluster members in `clusters_by_vertex` must be marked correctly as a
// `ClusterMember` prior to calling this function.
void DetermineOutliersAndHubs(
    const internal::NeighborOrder& adjacency_list,
    pbbs::sequence<VertexType>* clusters_by_vertex) {
  par_for(0, clusters_by_vertex->size(), [&](const size_t i) {
    VertexType* const vertex_type{&(*clusters_by_vertex)[i]};
    if (!std::holds_alternative<ClusterMember>(*vertex_type)) {
      // Determine whether remaining vertex i is a hub or outlier.

      const auto& neighbors{adjacency_list[i]};
      bool is_hub{false};
      // `candidate_cluster` holds a cluster ID that vertex i is adjacent to, or
      // UINT_E_MAX before such a cluster is found.
      uintE candidate_cluster{UINT_E_MAX};
      par_for(0, neighbors.size(), [&](const size_t j) {
        const ClusterMember* const neighbor_clusters{
          std::get_if<ClusterMember>(
              &(*clusters_by_vertex)[neighbors[j].neighbor])};
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
}

}  // namespace

bool operator==(const ClusterMember& a, const ClusterMember& b) {
  return a.clusters == b.clusters;
}

bool operator==(const Hub&, const Hub&) {
  return true;
}

bool operator==(const Outlier&, const Outlier&) {
  return true;
}

bool operator==(const Clustering& a, const Clustering& b) {
  return std::tie(a.num_clusters, a.clusters_by_vertex)
    == std::tie(b.num_clusters, b.clusters_by_vertex);
}

Clustering ScanIndex::Cluster(const uint64_t mu, const float epsilon) const {
  const pbbs::sequence<uintE> cores{core_order_.GetCores(mu, epsilon)};
  if (cores.empty()) {
    // Nothing is a core. There are no clusters, and every vertex is an outlier.
    return Clustering {
        .num_clusters = 0,
        .clusters_by_vertex =
          pbbs::sequence<VertexType>{num_vertices_, VertexType{Outlier{}}}
    };
  }
  const pbbs::sequence<DirectedEdge> core_similar_incident_edges{
    GetSimilarIncidentEdges(neighbor_order_, cores, epsilon)};
  Clustering clustering{
    .num_clusters = 0,
    // Mark everything as an outlier for now.
    .clusters_by_vertex =
      pbbs::sequence<VertexType>{num_vertices_, VertexType{Outlier{}}}
  };
  GetClustersFromCores(
      num_vertices_, cores, core_similar_incident_edges, &clustering);
  DetermineOutliersAndHubs(neighbor_order_, &clustering.clusters_by_vertex);
  return clustering;
}

}  // namespace scan
