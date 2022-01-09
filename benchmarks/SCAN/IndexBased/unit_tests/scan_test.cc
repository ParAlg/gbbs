#include "benchmarks/SCAN/IndexBased/scan.h"
#include "benchmarks/SCAN/IndexBased/scan_helpers.h"

#include <math.h>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "benchmarks/SCAN/IndexBased/similarity_measure.h"
#include "benchmarks/SCAN/IndexBased/unit_tests/similarity_measure_test_utils.h"
#include "benchmarks/SCAN/IndexBased/utils.h"
#include "gbbs/bridge.h"
#include "gbbs/graph.h"
#include "gbbs/helpers/undirected_edge.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gbbs/vertex.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace gbbs {

using ::testing::AllOf;
using ::testing::ElementsAre;
using ::testing::Eq;
using ::testing::Field;
using ::testing::FloatEq;
using ::testing::IsEmpty;
using ::testing::UnorderedElementsAre;

using ClusterList = std::set<std::set<uintE>>;
using VertexList = std::vector<uintE>;

namespace gt = graph_test;
namespace i = indexed_scan;
namespace ii = indexed_scan::internal;
namespace s = scan;

namespace {

// Googletest-like equality matcher for `indexed_scan::internal::CoreThreshold`.
auto CoreThresholdEq(const uintE expected_vertex,
                     const float expected_threshold) {
  return AllOf(
      Field(&ii::CoreThreshold::vertex_id, Eq(expected_vertex)),
      Field(&ii::CoreThreshold::threshold, FloatEq(expected_threshold)));
}

// Checks that `clustering` has the expected clusters and returns true if the
// check passes.
//
// The implementation may be inefficient, so don't expect this function to work
// well on large graphs.
//
// Arguments:
//   clustering
//     Clustering to check.
//   expected_clusters
//     List of distinct clusters that we expect to see, where each cluster is
//     given as a list of vertex IDs.
bool CheckClustering(const i::Clustering& clustering,
                     const ClusterList& expected_clusters) {
  std::unordered_map<uintE, std::set<uintE>> actual_clusters_map;
  actual_clusters_map.reserve(expected_clusters.size());
  for (size_t v = 0; v < clustering.size(); v++) {
    const uintE cluster_id{clustering[v]};
    if (cluster_id != i::kUnclustered) {
      actual_clusters_map[cluster_id].emplace(v);
    }
  }
  ClusterList actual_clusters;
  for (auto& cluster_kv : actual_clusters_map) {
    actual_clusters.emplace(std::move(cluster_kv.second));
  }

  if (actual_clusters != expected_clusters) {
    std::cerr << "Clusters don't match. Actual clustering:\n"
              << scan::ClusteringToString(clustering) << '\n';
    return false;
  } else {
    return true;
  }
}

// Checks that `DetermineUnclusteredType` identifies hubs and outliers as
// expected.
//
// TODO(tomtseng): It would be cleaner to have a unit test that tests
// `scan::DetermineUnclusteredType` separately instead of having this function.
//
// Arguments:
//   graph
//     Input graph.
//   clustering
//     Valid clustering on the graph.
//   expected_hubs
//     List of vertices that we expect to be hubs.
//   expected_outliers
//     List of vertices that we expect to be outliers.
void CheckUnclusteredVertices(
    symmetric_graph<symmetric_vertex, gbbs::empty>* graph,
    const i::Clustering& clustering, const VertexList& expected_hubs,
    const VertexList& expected_outliers) {
  for (const auto v : expected_hubs) {
    EXPECT_EQ(clustering[v], scan::kUnclustered);
    EXPECT_EQ(
        scan::DetermineUnclusteredType(clustering, graph->get_vertex(v), v),
        scan::UnclusteredType::kHub);
  }
  for (const auto v : expected_outliers) {
    EXPECT_EQ(clustering[v], scan::kUnclustered);
    EXPECT_EQ(
        scan::DetermineUnclusteredType(clustering, graph->get_vertex(v), v),
        scan::UnclusteredType::kOutlier);
  }
}

}  // namespace

TEST(ScanSubroutines, NullGraph) {
  const size_t kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const ii::NeighborOrder neighbor_order{&graph, scan::CosineSimilarity{}};
  EXPECT_THAT(neighbor_order, IsEmpty());

  const auto core_order{ii::ComputeCoreOrder(neighbor_order)};
  EXPECT_THAT(core_order, IsEmpty());
}

TEST(ScanSubroutines, EmptyGraph) {
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const ii::NeighborOrder neighbor_order{&graph, scan::CosineSimilarity{}};
  EXPECT_EQ(neighbor_order.size(), kNumVertices);

  EXPECT_THAT(neighbor_order[0], IsEmpty());
  for (const auto& vertex_order : neighbor_order) {
    EXPECT_THAT(vertex_order, IsEmpty());
  }

  auto core_order{ii::ComputeCoreOrder(neighbor_order)};
  ASSERT_EQ(core_order.size(), 2);
  EXPECT_THAT(core_order[0], IsEmpty());
  EXPECT_THAT(core_order[1], IsEmpty());
}

TEST(ScanSubroutines, BasicUsage) {
  // Graph diagram:
  //     0 --- 1 -- 2 -- 5
  //           |   /|
  //           | /  |
  //           3 -- 4
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {1, 2}, {1, 3}, {2, 3}, {2, 4}, {2, 5}, {3, 4},
  };
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const ii::NeighborOrder neighbor_order{&graph, scan::CosineSimilarity{}};
  ASSERT_EQ(neighbor_order.size(), kNumVertices);
  EXPECT_THAT(neighbor_order[0],
              ElementsAre(EdgeSimilarityEq(0, 1, 2.0 / sqrt(8))));
  EXPECT_THAT(neighbor_order[1],
              ElementsAre(EdgeSimilarityEq(1, 3, 3.0 / sqrt(16)),
                          EdgeSimilarityEq(1, 0, 2.0 / sqrt(8)),
                          EdgeSimilarityEq(1, 2, 3.0 / sqrt(20))));
  EXPECT_THAT(neighbor_order[2],
              ElementsAre(EdgeSimilarityEq(2, 3, 4.0 / sqrt(20)),
                          EdgeSimilarityEq(2, 4, 3.0 / sqrt(15)),
                          EdgeSimilarityEq(2, 1, 3.0 / sqrt(20)),
                          EdgeSimilarityEq(2, 5, 2.0 / sqrt(10))));
  EXPECT_THAT(neighbor_order[3],
              ElementsAre(EdgeSimilarityEq(3, 2, 4.0 / sqrt(20)),
                          EdgeSimilarityEq(3, 4, 3.0 / sqrt(12)),
                          EdgeSimilarityEq(3, 1, 3.0 / sqrt(16))));
  EXPECT_THAT(neighbor_order[4],
              ElementsAre(EdgeSimilarityEq(4, 3, 3.0 / sqrt(12)),
                          EdgeSimilarityEq(4, 2, 3.0 / sqrt(15))));
  EXPECT_THAT(neighbor_order[5],
              ElementsAre(EdgeSimilarityEq(5, 2, 2.0 / sqrt(10))));

  {
    const auto core_order{ii::ComputeCoreOrder(neighbor_order)};
    EXPECT_EQ(core_order.size(), 6);
    EXPECT_THAT(core_order[0], IsEmpty());
    EXPECT_THAT(core_order[1], IsEmpty());
    ASSERT_EQ(core_order[2].size(), 6);
    EXPECT_THAT(core_order[2].cut(0, 2),
                UnorderedElementsAre(CoreThresholdEq(2, 4.0 / sqrt(20)),
                                     CoreThresholdEq(3, 4.0 / sqrt(20))));
    EXPECT_THAT(core_order[2].cut(2, core_order[2].size()),
                ElementsAre(CoreThresholdEq(4, 3.0 / sqrt(12)),
                            CoreThresholdEq(1, 3.0 / sqrt(16)),
                            CoreThresholdEq(0, 2.0 / sqrt(8)),
                            CoreThresholdEq(5, 2.0 / sqrt(10))));
    ASSERT_EQ(core_order[3].size(), 4);
    EXPECT_THAT(core_order[3][0], CoreThresholdEq(3, 3.0 / sqrt(12)));
    EXPECT_THAT(core_order[3].cut(1, 3),
                UnorderedElementsAre(CoreThresholdEq(2, 3.0 / sqrt(15)),
                                     CoreThresholdEq(4, 3.0 / sqrt(15))));
    EXPECT_THAT(core_order[3].cut(3, core_order[3].size()),
                ElementsAre(CoreThresholdEq(1, 2.0 / sqrt(8))));
    ASSERT_EQ(core_order[4].size(), 3);
    EXPECT_THAT(core_order[4][0], CoreThresholdEq(3, 3.0 / sqrt(16)));
    EXPECT_THAT(core_order[4].cut(1, core_order[4].size()),
                UnorderedElementsAre(CoreThresholdEq(1, 3.0 / sqrt(20)),
                                     CoreThresholdEq(2, 3.0 / sqrt(20))));
    EXPECT_THAT(core_order[5], ElementsAre(CoreThresholdEq(2, 2.0 / sqrt(10))));
  }

  {
    const ii::CoreOrder core_order{neighbor_order};
    {
      const uint64_t kMu{2};
      const float kEpsilon{0.5};
      const sequence<uintE> cores{core_order.GetCores(kMu, kEpsilon)};
      EXPECT_THAT(cores, UnorderedElementsAre(0, 1, 2, 3, 4, 5));
    }
    {
      const uint64_t kMu{2};
      const float kEpsilon{0.8};
      const sequence<uintE> cores{core_order.GetCores(kMu, kEpsilon)};
      EXPECT_THAT(cores, UnorderedElementsAre(2, 3, 4));
    }
    {
      const uint64_t kMu{2};
      const float kEpsilon{0.88};
      const sequence<uintE> cores{core_order.GetCores(kMu, kEpsilon)};
      EXPECT_THAT(cores, UnorderedElementsAre(2, 3));
    }
    {
      const uint64_t kMu{2};
      const float kEpsilon{0.9};
      const sequence<uintE> cores{core_order.GetCores(kMu, kEpsilon)};
      EXPECT_THAT(cores, IsEmpty());
    }
  }
}

TEST(ScanSubroutines, DisconnectedGraph) {
  // Graph diagram:
  //     0 -- 1    2    3 -- 4 -- 5
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {3, 4}, {4, 5},
  };
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const ii::NeighborOrder neighbor_order{&graph, scan::CosineSimilarity{}};
  ASSERT_EQ(neighbor_order.size(), kNumVertices);
  EXPECT_THAT(neighbor_order[0], ElementsAre(EdgeSimilarityEq(0, 1, 1.0)));
  EXPECT_THAT(neighbor_order[1], ElementsAre(EdgeSimilarityEq(1, 0, 1.0)));
  EXPECT_THAT(neighbor_order[2], IsEmpty());
  EXPECT_THAT(neighbor_order[3],
              ElementsAre(EdgeSimilarityEq(3, 4, 2.0 / sqrt(6))));
  EXPECT_THAT(neighbor_order[4],
              UnorderedElementsAre(EdgeSimilarityEq(4, 3, 2.0 / sqrt(6)),
                                   EdgeSimilarityEq(4, 5, 2.0 / sqrt(6))));
  EXPECT_THAT(neighbor_order[5],
              ElementsAre(EdgeSimilarityEq(5, 4, 2.0 / sqrt(6))));

  const auto core_order{ii::ComputeCoreOrder(neighbor_order)};
  EXPECT_EQ(core_order.size(), 4);
  EXPECT_THAT(core_order[0], IsEmpty());
  EXPECT_THAT(core_order[1], IsEmpty());
  ASSERT_EQ(core_order[2].size(), 5);
  EXPECT_THAT(
      core_order[2].cut(0, 2),
      UnorderedElementsAre(CoreThresholdEq(0, 1.0), CoreThresholdEq(1, 1.0)));
  EXPECT_THAT(core_order[2].cut(2, core_order[2].size()),
              UnorderedElementsAre(CoreThresholdEq(3, 2.0 / sqrt(6)),
                                   CoreThresholdEq(4, 2.0 / sqrt(6)),
                                   CoreThresholdEq(5, 2.0 / sqrt(6))));
  EXPECT_THAT(core_order[3], ElementsAre(CoreThresholdEq(4, 2.0 / sqrt(6))));
}

TEST(Cluster, NullGraph) {
  const size_t kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const i::Index index{&graph};
  constexpr float kMu{2};
  constexpr float kEpsilon{0.5};
  i::Clustering clustering{index.Cluster(kMu, kEpsilon)};
  EXPECT_THAT(clustering, IsEmpty());
}

TEST(Cluster, EmptyGraph) {
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const i::Index index{&graph};
  constexpr float kMu{2};
  constexpr float kEpsilon{0.5};
  i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

  const ClusterList kExpectedClusters{};
  const VertexList kExpectedHubs{};
  const VertexList kExpectedOutliers{0, 1, 2, 3, 4, 5};
  EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
  CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                           kExpectedOutliers);
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
  const i::Index index{&graph};

  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.5};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{0, 1, 2, 3, 4, 5}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{};
    EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
    CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                             kExpectedOutliers);
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.7};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{0, 1, 2, 3, 4}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{5};
    EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
    CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                             kExpectedOutliers);
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.73};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{1, 2, 3, 4}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 5};
    EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
    CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                             kExpectedOutliers);
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.88};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{2, 3}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 1, 4, 5};
    EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
    CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                             kExpectedOutliers);
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.95};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 1, 2, 3, 4, 5};
    EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
    CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                             kExpectedOutliers);
  }
  {
    constexpr uint64_t kMu{4};
    constexpr float kEpsilon{0.7};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{1, 2, 3, 4}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 5};
    EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
    CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                             kExpectedOutliers);
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
  const i::Index index{&graph};

  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.95};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{0, 1}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{2, 3, 4, 5};
    EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
    CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                             kExpectedOutliers);
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.8};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{0, 1}, {3, 4, 5}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{2};
    EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
    CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                             kExpectedOutliers);
  }
  {
    constexpr uint64_t kMu{3};
    constexpr float kEpsilon{0.8};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{3, 4, 5}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 1, 2};
    EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
    CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                             kExpectedOutliers);
  }
  {
    constexpr uint64_t kMu{4};
    constexpr float kEpsilon{0.8};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 1, 2, 3, 4, 5};
    EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
    CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                             kExpectedOutliers);
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
  const i::Index index{&graph};
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.73};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{1, 2, 3}, {5, 6, 7}};
    const VertexList kExpectedHubs{4};
    const VertexList kExpectedOutliers{0, 8};
    EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
    CheckUnclusteredVertices(&graph, clustering, kExpectedHubs,
                             kExpectedOutliers);
  }
}
}  // namespace gbbs
