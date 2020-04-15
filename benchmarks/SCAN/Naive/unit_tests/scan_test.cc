#include "benchmarks/SCAN/Naive/scan.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ligra/graph_test_utils.h"
#include "ligra/undirected_edge.h"
#include "pbbslib/seq.h"

namespace gt = graph_test;
namespace n = naive_scan;

// Like naive_scan::Clustering, but using `std::vector` instead of
// `pbbs::sequence` for convenience.
using ExpectedClusters = std::vector<std::vector<uintE>>;

namespace {

// Checks that `clustering` has the expected clusters and returns true if the
// check passes.
bool CheckClustering(
    const n::Clustering& actual_clusters,
    const ExpectedClusters& expected_clusters) {
  return false;
}

}  // namespace

TEST(Cluster, NullGraph) {
  const size_t kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  constexpr float kMu{2};
  constexpr float kEpsilon{0.5};
  const n::Clustering clusters{n::Cluster(&graph, kMu, kEpsilon)};
  const ExpectedClusters expected_clusters{};
  EXPECT_TRUE(CheckClustering(clusters, expected_clusters));
}

// TEST(Cluster, EmptyGraph) {
//   const size_t kNumVertices{6};
//   const std::unordered_set<UndirectedEdge> kEdges{};
//   auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

//   const i::Index index{&graph};
//   constexpr float kMu{2};
//   constexpr float kEpsilon{0.5};
//   i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//   const ClusterList kExpectedClusters{};
//   const VertexList kExpectedHubs{};
//   const VertexList kExpectedOutliers{0, 1, 2, 3, 4, 5};
//   EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//   CheckUnclusteredVertices(
//       &graph, clustering, kExpectedHubs, kExpectedOutliers);
// }

// TEST(Cluster, BasicUsage) {
//   // Graph diagram with structural similarity scores labeled:
//   //       .71    .67  .63
//   //     0 --- 1 ---- 2 -- 5
//   //           |    / |
//   //       .75 |  .89 | .77
//   //           | /    |
//   //           3 ---- 4
//   //             .87
//   const size_t kNumVertices{6};
//   const std::unordered_set<UndirectedEdge> kEdges{
//     {0, 1},
//     {1, 2},
//     {1, 3},
//     {2, 3},
//     {2, 4},
//     {2, 5},
//     {3, 4},
//   };
//   auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
//   const i::Index index{&graph};

//   {
//     constexpr uint64_t kMu{2};
//     constexpr float kEpsilon{0.5};
//     const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//     const ClusterList kExpectedClusters{{0, 1, 2, 3, 4, 5}};
//     const VertexList kExpectedHubs{};
//     const VertexList kExpectedOutliers{};
//     EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//     CheckUnclusteredVertices(
//         &graph, clustering, kExpectedHubs, kExpectedOutliers);
//   }
//   {
//     constexpr uint64_t kMu{2};
//     constexpr float kEpsilon{0.7};
//     const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//     const ClusterList kExpectedClusters{{0, 1, 2, 3, 4}};
//     const VertexList kExpectedHubs{};
//     const VertexList kExpectedOutliers{5};
//     EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//     CheckUnclusteredVertices(
//         &graph, clustering, kExpectedHubs, kExpectedOutliers);
//   }
//   {
//     constexpr uint64_t kMu{2};
//     constexpr float kEpsilon{0.73};
//     const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//     const ClusterList kExpectedClusters{{1, 2, 3, 4}};
//     const VertexList kExpectedHubs{};
//     const VertexList kExpectedOutliers{0, 5};
//     EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//     CheckUnclusteredVertices(
//         &graph, clustering, kExpectedHubs, kExpectedOutliers);
//   }
//   {
//     constexpr uint64_t kMu{2};
//     constexpr float kEpsilon{0.88};
//     const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//     const ClusterList kExpectedClusters{{2, 3}};
//     const VertexList kExpectedHubs{};
//     const VertexList kExpectedOutliers{0, 1, 4, 5};
//     EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//     CheckUnclusteredVertices(
//         &graph, clustering, kExpectedHubs, kExpectedOutliers);
//   }
//   {
//     constexpr uint64_t kMu{2};
//     constexpr float kEpsilon{0.95};
//     const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//     const ClusterList kExpectedClusters{};
//     const VertexList kExpectedHubs{};
//     const VertexList kExpectedOutliers{0, 1, 2, 3, 4, 5};
//     EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//     CheckUnclusteredVertices(
//         &graph, clustering, kExpectedHubs, kExpectedOutliers);
//   }
//   {
//     constexpr uint64_t kMu{4};
//     constexpr float kEpsilon{0.7};
//     const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//     const ClusterList kExpectedClusters{{1, 2, 3, 4}};
//     const VertexList kExpectedHubs{};
//     const VertexList kExpectedOutliers{0, 5};
//     EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//     CheckUnclusteredVertices(
//         &graph, clustering, kExpectedHubs, kExpectedOutliers);
//   }
// }

// TEST(Cluster, DisconnectedGraph) {
//   // Graph diagram with structural similarity scores labeled:
//   //       1.0            .82  .82
//   //     0 -- 1    2    3 -- 4 -- 5
//   const size_t kNumVertices{6};
//   const std::unordered_set<UndirectedEdge> kEdges{
//     {0, 1},
//     {3, 4},
//     {4, 5},
//   };
//   auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
//   const i::Index index{&graph};

//   {
//     constexpr uint64_t kMu{2};
//     constexpr float kEpsilon{0.95};
//     const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//     const ClusterList kExpectedClusters{{0, 1}};
//     const VertexList kExpectedHubs{};
//     const VertexList kExpectedOutliers{2, 3, 4, 5};
//     EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//     CheckUnclusteredVertices(
//         &graph, clustering, kExpectedHubs, kExpectedOutliers);
//   }
//   {
//     constexpr uint64_t kMu{2};
//     constexpr float kEpsilon{0.8};
//     const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//     const ClusterList kExpectedClusters{{0, 1}, {3, 4, 5}};
//     const VertexList kExpectedHubs{};
//     const VertexList kExpectedOutliers{2};
//     EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//     CheckUnclusteredVertices(
//         &graph, clustering, kExpectedHubs, kExpectedOutliers);
//   }
//   {
//     constexpr uint64_t kMu{3};
//     constexpr float kEpsilon{0.8};
//     const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//     const ClusterList kExpectedClusters{{3, 4, 5}};
//     const VertexList kExpectedHubs{};
//     const VertexList kExpectedOutliers{0, 1, 2};
//     EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//     CheckUnclusteredVertices(
//         &graph, clustering, kExpectedHubs, kExpectedOutliers);
//   }
//   {
//     constexpr uint64_t kMu{4};
//     constexpr float kEpsilon{0.8};
//     const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//     const ClusterList kExpectedClusters{};
//     const VertexList kExpectedHubs{};
//     const VertexList kExpectedOutliers{0, 1, 2, 3, 4, 5};
//     EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//     CheckUnclusteredVertices(
//         &graph, clustering, kExpectedHubs, kExpectedOutliers);
//   }
// }

// TEST(Cluster, TwoClusterGraph) {
//   // Graph diagram with structural similarity scores labeled:
//   //           2                   6
//   //     .87 /   \ .87       .87 /   \ .87
//   //        /     \             /     \.
//   // 0 --- 1 ----- 3 --- 4 --- 5 ----- 7 --- 8
//   //   .71    .75    .58   .58    .75   .71
//   const size_t kNumVertices{9};
//   const std::unordered_set<UndirectedEdge> kEdges{
//     {0, 1},
//     {1, 2},
//     {1, 3},
//     {2, 3},
//     {3, 4},
//     {4, 5},
//     {5, 6},
//     {5, 7},
//     {6, 7},
//     {7, 8},
//   };
//   auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
//   const i::Index index{&graph};
//   {
//     constexpr uint64_t kMu{2};
//     constexpr float kEpsilon{0.73};
//     const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

//     const ClusterList kExpectedClusters{{1, 2, 3}, {5, 6, 7}};
//     const VertexList kExpectedHubs{4};
//     const VertexList kExpectedOutliers{0, 8};
//     EXPECT_TRUE(CheckClustering(clustering, kExpectedClusters));
//     CheckUnclusteredVertices(
//         &graph, clustering, kExpectedHubs, kExpectedOutliers);
//   }
// }
