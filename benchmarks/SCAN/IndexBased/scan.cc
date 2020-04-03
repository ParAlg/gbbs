#include "benchmarks/SCAN/IndexBased/scan.h"

#include <algorithm>
#include <sstream>
#include <tuple>
#include <utility>

#include "benchmarks/Connectivity/UnionFind/union_find_rules.h"
#include "ligra/bridge.h"
#include "pbbslib/parallel.h"

namespace indexed_scan {

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
  timer function_timer{"Get similar incident edges time"};

  pbbs::sequence<size_t> epsilon_neighborhood_offsets{
      pbbs::map<size_t>(
          vertices,
          [&](const uintE vertex) {
            // Get the number of neighbors of `vertex` that have at least
            // `epsilon` structural similarity with the vertex.
            return internal::BinarySearch(
                neighbor_order[vertex],
                [epsilon](const internal::EdgeSimilarity& es) {
                  return es.similarity >= epsilon;
                });
          })};
  const size_t num_incident_edges{
    pbbslib::scan_add_inplace(epsilon_neighborhood_offsets)};
  pbbs::sequence<DirectedEdge> incident_edges{
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

  internal::ReportTime(function_timer);
  return incident_edges;
}

// Given a list of edges of sufficient similarity between core vertices,
// identifies the clusters for the core vertices and populates `clustering`
// accordingly. Non-core vertices in `clustering` are left unchanged.
//
// Arguments:
//   cores: List of IDs of vertices that are SCAN cores.
//   core_to_core_edges: Directed edges between cores of at least epsilon
//     similarity.
//   clustering: Output. Will be updated with the cluster IDs only for the core
//   vertices.
void ClusterCores(
    const pbbs::sequence<uintE>& cores,
    const pbbs::range<DirectedEdge*>& core_to_core_edges,
    Clustering* clustering) {
  timer function_timer{"Cluster cores time"};
  par_for(0, cores.size(), [&](const size_t i) {
      const uintE core{cores[i]};
      (*clustering)[core] = core;
  });
  constexpr auto find{find_variants::find_compress};
  auto unite{unite_variants::Unite<decltype(find)>{find}};
  par_for(0, core_to_core_edges.size(), [&](const size_t i) {
      const auto [u, v] = core_to_core_edges[i];
      if (u > v) {
        unite(u, v, *clustering);
      }
  });
  par_for(0, cores.size(), [&](const size_t i) {
      const uintE core{cores[i]};
      (*clustering)[core] = find(core, *clustering);
  });
  internal::ReportTime(function_timer);
}

// Given `clustering` that's populated for all core vertices and a list of
// [core -> non-core] edges of sufficient similarity, populate `clustering` for
// all non-core vertices.
//
// Arguments:
//   num_vertices: Number of vertices in the whole graph.
//   core_to_noncore_edges: Directed edges of at least epsilon similarity from
//     core to non-core vertices.
//   clustering: Output. When function is called, core vertices must be already
//     be assigned a cluster ID. After function is called, non-core vertices
//     that belong to a cluster will be marked accordingly.
void AttachNoncoresToClusters(
    const uintE num_vertices,
    const pbbs::range<DirectedEdge*>& core_to_noncore_edges,
    Clustering* clustering) {
  timer function_timer{"Attach non-cores to clusters time"};
  // Attach each non-core to the same cluster as an arbitrary adjacent core.
  //
  // We could do something smarter like attach each non-core to its most
  // structurally similar adjacent core, but let's keep things simple for now.
  par_for(0, core_to_noncore_edges.size(), [&](const size_t i) {
    uintE core, noncore;
    std::tie(core, noncore) = core_to_noncore_edges[i];
    auto* noncore_cluster_address{&(*clustering)[noncore]};
    if (*noncore_cluster_address == kUnclustered) {
      const uintE cluster{(*clustering)[core]};
      pbbs::atomic_compare_and_swap(
          noncore_cluster_address, kUnclustered, cluster);
    }
  });
  internal::ReportTime(function_timer);
}

// Given core vertices, computes the SCAN clusters.
//
// Arguments:
//   num_vertices: Number of vertices in the graph.
//   cores: List of IDs for SCAN core vertices.
//   core_similar_incident_edges: List of directed edges incident on cores with
//     similarity at least epsilon, where the first endpoint of each edge is a
//     core.
Clustering GetClustersFromCores(
    const uintE num_vertices,
    const pbbs::sequence<uintE>& cores,
    const pbbs::sequence<DirectedEdge>& core_similar_incident_edges) {
  timer preprocessing_timer{"Get clusters from cores - preprocessing time"};

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
    pbbs::split_two(
        core_similar_incident_edges,
        pbbs::delayed_seq<bool>(
          core_similar_incident_edges.size(),
          [&](const size_t i) {
            // Only need to check second endpoint. First endpoint is a core.
            return !cores_set.contains(core_similar_incident_edges[i].second);
          }));

  internal::ReportTime(preprocessing_timer);

  Clustering clustering(num_vertices, kUnclustered);

  pbbs::range<DirectedEdge*> core_to_core_edges{
    partitioned_edges.slice(0, num_core_to_core_edges)};
  ClusterCores(cores, core_to_core_edges, &clustering);

  pbbs::range<DirectedEdge*> core_to_noncore_edges{
    partitioned_edges.slice(num_core_to_core_edges, partitioned_edges.size())};
  AttachNoncoresToClusters(num_vertices, core_to_noncore_edges, &clustering);

  return clustering;
}

}  // namespace

std::ostream& operator<<(std::ostream& os, UnclusteredType unclustered_type) {
  switch (unclustered_type) {
    case UnclusteredType::kHub:
      os << "hub";
      break;
    case UnclusteredType::kOutlier:
      os << "outlier";
      break;
  }
  return os;
}

Clustering Index::Cluster(const uint64_t mu, const float epsilon) const {
  const pbbs::sequence<uintE> cores{core_order_.GetCores(mu, epsilon)};
  if (cores.empty()) {
    // Nothing is a core. There are no clusters, and every vertex is an outlier.
    return Clustering(num_vertices_, kUnclustered);
  }
  const pbbs::sequence<DirectedEdge> core_similar_incident_edges{
    GetSimilarIncidentEdges(neighbor_order_, cores, epsilon)};
  return GetClustersFromCores(
      num_vertices_, cores, core_similar_incident_edges);
}

size_t CleanClustering(Clustering* clustering) {
  const size_t num_vertices{clustering->size()};
  pbbs::sequence<uintE> cluster_relabel_map(num_vertices, 0U);
  par_for(0, num_vertices, [&](const size_t i) {
    const uintE cluster_id{(*clustering)[i]};
    if (cluster_id != kUnclustered && cluster_relabel_map[cluster_id] == 0) {
      cluster_relabel_map[cluster_id] = 1;
    }
  });
  const size_t num_clusters{pbbslib::scan_add_inplace(cluster_relabel_map)};
  par_for(0, num_vertices, [&](const size_t i) {
    const uintE cluster_id{(*clustering)[i]};
    if (cluster_id != kUnclustered) {
      (*clustering)[i] = cluster_relabel_map[cluster_id];
    }
  });
  return num_clusters;
}

std::string ClusteringToString(const Clustering& clustering) {
  std::ostringstream str;
  str << "{";
  for (size_t i = 0; i < clustering.size(); i++) {
    str << ' ' << i << ':';
    str << (clustering[i] == kUnclustered
            ? "n/a" : std::to_string(clustering[i]));
  }
  str << " }";
  return str.str();
}

}  // namespace indexed_scan
