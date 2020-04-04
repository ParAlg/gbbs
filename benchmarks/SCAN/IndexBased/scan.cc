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
using VertexSet =
  sparse_table<uintE, pbbslib::empty, decltype(&pbbslib::hash64_2)>;

// Creates a `VertexSet` for holding up to `capacity` elements.
VertexSet MakeVertexSet(const size_t capacity) {
  return make_sparse_table<uintE, pbbslib::empty, decltype(&pbbslib::hash64_2)>(
      // Adding 1 avoids having small tables completely full.
      capacity + 1, {UINT_E_MAX, pbbslib::empty{}}, pbbslib::hash64_2);
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
    const pbbs::sequence<size_t>& core_similar_edge_counts,
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
    const pbbs::sequence<size_t>& core_similar_edge_counts,
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

  timer preprocessing_timer{"Cluster - additional preprocessing time"};
  pbbs::sequence<size_t> core_similar_edge_counts{
      pbbs::map<size_t>(
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
  AttachNoncoresToClusters(
      neighbor_order_, cores, core_similar_edge_counts, &clustering);
  return clustering;
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
