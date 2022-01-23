#include "benchmarks/SCAN/Naive/scan.h"

#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>

#include "gbbs/helpers/undirected_edge.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace gbbs {
namespace gt = graph_test;
namespace n = naive_scan;

// Like naive_scan::Clustering, but uses `std::vector` instead of
// `sequence` for convenience.
using ClusteringArray = std::vector<std::vector<uintE>>;

namespace {

// Converts a `naive_scan::Clustering` or `ClusteringArray` to a string. Useful
// for debugging.
template <class Clustering>
std::string ClusteringToString(const Clustering& clustering) {
  constexpr size_t kWidth{2};
  std::ostringstream str;
  str << std::setw(kWidth);
  str << "{\n";
  for (size_t i = 0; i < clustering.size(); i++) {
    str << "  " << i << ":";
    if (clustering[i].empty()) {
      str << " n/a";
    } else {
      for (uintE cluster : clustering[i]) {
        str << " " << cluster;
      }
    }
    str << "\n";
  }
  str << "}";
  return str.str();
}

// Checks that `actual_clustering` has the expected clusters.
void CheckClustering(const n::Clustering& actual_clustering,
                     const ClusteringArray& expected_clustering) {
  ASSERT_EQ(actual_clustering.size(), expected_clustering.size());
  // <cluster IDs in actual_clustering -> cluster IDs in expected_clustering>
  // map
  std::unordered_map<uintE, uintE> cluster_map;

  for (size_t i = 0; i < actual_clustering.size(); i++) {
    // Every used cluster ID will appear by itself on a core vertex, so finding
    // each actual_clustering[*] of size 1 will find all cluster IDs.
    if (actual_clustering[i].size() == 1) {
      ASSERT_EQ(expected_clustering[i].size(), 1);
      const uintE actual_cluster_id{actual_clustering[i][0]};
      const uintE expected_cluster_id{expected_clustering[i][0]};
      const auto it{cluster_map.find(actual_cluster_id)};
      if (it == cluster_map.end()) {
        cluster_map.emplace(actual_cluster_id, expected_cluster_id);
      } else {
        // If a mapping for this cluster ID already exists, it must be
        // consistent with `expected_clustering`.
        EXPECT_EQ(it->second, expected_cluster_id);
      }
    }
  }

  ClusteringArray actual_clustering_remapped{actual_clustering.size()};
  for (size_t i = 0; i < actual_clustering.size(); i++) {
    for (uintE actual_cluster : actual_clustering[i]) {
      actual_clustering_remapped[i].emplace_back(
          cluster_map.at(actual_cluster));
    }
  }

  EXPECT_EQ(actual_clustering_remapped, expected_clustering);
}

}  // namespace

TEST(Cluster, NullGraph) {
  const size_t kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  constexpr float kMu{2};
  constexpr float kEpsilon{0.5};
  const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

  const ClusteringArray kExpectedClustering{};
  CheckClustering(clustering, kExpectedClustering);
}

TEST(Cluster, EmptyGraph) {
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  constexpr float kMu{2};
  constexpr float kEpsilon{0.5};
  const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

  const ClusteringArray kExpectedClustering(kNumVertices, std::vector<uintE>{});
  CheckClustering(clustering, kExpectedClustering);
}

TEST(Cluster, BasicUsage) {
  // Graph diagram with structural similarity scores labeled:
  //       .71    .67  .63
  //     0 --- 1 ---- 2 -- 5
  //           |    / |
  //       .75 |  .89 | .77
  //           | /    |
  //           3 ---- 4
  //             .87
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {1, 2}, {1, 3}, {2, 3}, {2, 4}, {2, 5}, {3, 4},
  };
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.5};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering(kNumVertices,
                                              std::vector<uintE>{0});
    CheckClustering(clustering, kExpectedClustering);
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.7};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering{{0}, {0}, {0}, {0}, {0}, {}};
    CheckClustering(clustering, kExpectedClustering);
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.73};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering{{}, {0}, {0}, {0}, {0}, {}};
    CheckClustering(clustering, kExpectedClustering);
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.88};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering{{}, {}, {0}, {0}, {}, {}};
    CheckClustering(clustering, kExpectedClustering);
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.95};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering(kNumVertices,
                                              std::vector<uintE>{});
    CheckClustering(clustering, kExpectedClustering);
  }
  {
    constexpr uint64_t kMu{4};
    constexpr float kEpsilon{0.7};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering{{}, {0}, {0}, {0}, {0}, {}};
    CheckClustering(clustering, kExpectedClustering);
  }
}

TEST(Cluster, DisconnectedGraph) {
  // Graph diagram with structural similarity scores labeled:
  //       1.0            .82  .82
  //     0 -- 1    2    3 -- 4 -- 5
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {3, 4}, {4, 5},
  };
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.95};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering{{0}, {0}, {}, {}, {}, {}};
    CheckClustering(clustering, kExpectedClustering);
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.8};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering{{0}, {0}, {}, {1}, {1}, {1}};
    CheckClustering(clustering, kExpectedClustering);
  }
  {
    constexpr uint64_t kMu{3};
    constexpr float kEpsilon{0.8};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering{{}, {}, {}, {0}, {0}, {0}};
    CheckClustering(clustering, kExpectedClustering);
  }
  {
    constexpr uint64_t kMu{4};
    constexpr float kEpsilon{0.8};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering(kNumVertices,
                                              std::vector<uintE>{});
    CheckClustering(clustering, kExpectedClustering);
  }
}

TEST(Cluster, TwoClusterGraph) {
  // Graph diagram with structural similarity scores labeled:
  //           2                   6
  //     .87 /   \ .87       .87 /   \ .87
  //        /     \             /     \.
  // 0 --- 1 ----- 3 --- 4 --- 5 ----- 7 --- 8
  //   .71    .75    .58   .58    .75   .71
  const size_t kNumVertices{9};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {1, 2}, {1, 3}, {2, 3}, {3, 4},
      {4, 5}, {5, 6}, {5, 7}, {6, 7}, {7, 8},
  };
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.73};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering{{},  {0}, {0}, {0}, {},
                                              {1}, {1}, {1}, {}};
    CheckClustering(clustering, kExpectedClustering);
  }
  {
    constexpr uint64_t kMu{4};
    constexpr float kEpsilon{0.5};
    const n::Clustering clustering{n::Cluster(&graph, kMu, kEpsilon)};

    const ClusteringArray kExpectedClustering{{0}, {0}, {0}, {0}, {0, 1},
                                              {1}, {1}, {1}, {1}};
    CheckClustering(clustering, kExpectedClustering);
  }
}
}  // namespace gbbs
