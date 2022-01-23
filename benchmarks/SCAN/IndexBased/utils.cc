#include "benchmarks/SCAN/IndexBased/utils.h"

#include <sstream>

namespace gbbs {
namespace scan {

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

std::string ClusteringToString(const Clustering& clustering) {
  std::ostringstream str;
  str << "{";
  for (size_t i = 0; i < clustering.size(); i++) {
    str << ' ' << i << ':';
    str << (clustering[i] == kUnclustered ? "n/a"
                                          : std::to_string(clustering[i]));
  }
  str << " }";
  return str.str();
}

size_t CompactClustering(Clustering* clustering) {
  const size_t num_vertices{clustering->size()};
  sequence<uintE> cluster_relabel_map(num_vertices, 0U);
  parallel_for(0, num_vertices, [&](const size_t i) {
    const uintE cluster_id{(*clustering)[i]};
    if (cluster_id != kUnclustered && cluster_relabel_map[cluster_id] == 0) {
      cluster_relabel_map[cluster_id] = 1;
    }
  });
  const size_t num_clusters{parlay::scan_inplace(cluster_relabel_map)};
  parallel_for(0, num_vertices, [&](const size_t i) {
    const uintE cluster_id{(*clustering)[i]};
    if (cluster_id != kUnclustered) {
      (*clustering)[i] = cluster_relabel_map[cluster_id];
    }
  });
  return num_clusters;
}

}  // namespace scan
}  // namespace gbbs
