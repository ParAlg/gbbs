#include "benchmarks/SCAN/IndexBased/scan.h"
#include "benchmarks/SCAN/IndexBased/scan_helpers.h"

#include <math.h>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <optional>
#include <utility>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ligra/graph.h"
#include "ligra/graph_test_utils.h"
#include "ligra/undirected_edge.h"
#include "ligra/vertex.h"
#include "pbbslib/seq.h"

using ::testing::ElementsAre;
using ::testing::IsEmpty;
using ::testing::UnorderedElementsAre;

using ClusterList = std::set<std::set<uintE>>;
using VertexList = std::vector<uintE>;

namespace gt = graph_test;
namespace i = indexed_scan;
namespace ii = indexed_scan::internal;

namespace {

ii::NeighborSimilarity
MakeNeighborSimilarity(const uintE neighbor, const double similarity) {
  return {.neighbor = neighbor, .similarity = static_cast<float>(similarity)};
}

ii::CoreThreshold
MakeCoreThreshold(const uintE vertex_id, const double threshold) {
  return {.vertex_id = vertex_id, .threshold = static_cast<float>(threshold)};
}

// Checks that the `clustering` has the expected clusters. If an error is found,
// the function returns an error message as a string. If the clustering is
// matches what is expected, then `std::nullopt` is returned.
//
// The implementation may be inefficient, so don't expect this function to work
// well on large graphs.
//
// This function exists because checking the equivalence of two clusterings
// takes some work --- clusterings are equivalent under permutation of cluster
// IDs.
// An alternative interface would be to have this function take
// `(const Clustering& actual_clustering,
//   const Clustering& expected_clustering)`
// as input, but the existing interface is more concise invoke since
// constructing a `Clustering` is verbose.
//
// Arguments
// ---------
// clustering
//   Clustering to check.
// expected_clusters
//   List of distinct clusters that we expect to see, where each cluster is
//   given as a list of vertex IDs.
// expected_hubs
//   List of vertices that we expect to be hubs.
// expected_outliers
//   List of vertices that we expect to be outliers.
std::optional<std::string> CheckClusteringForErrors(
    const i::Clustering& clustering,
    const ClusterList& expected_clusters,
    const VertexList& expected_hubs,
    const VertexList& expected_outliers) {
  std::ostringstream error_msg;

  const size_t num_clusters{expected_clusters.size()};
  if (clustering.num_clusters != num_clusters) {
    error_msg << "Expected " << num_clusters << " clusters "
      << "but found " << clustering.num_clusters << " instead";
    return error_msg.str();
  }

  for (const uintE v : expected_hubs) {
    if (!std::holds_alternative<i::Hub>(clustering.clusters_by_vertex[v])) {
      error_msg << "Expected vertex " << v << " to be a hub, "
        << "but is " << clustering.clusters_by_vertex[v] << " instead";
      return error_msg.str();
    }
  }
  for (const uintE v : expected_outliers) {
    if (!std::holds_alternative<i::Outlier>(clustering.clusters_by_vertex[v])) {
      error_msg << "Expected vertex " << v << " to be an outlier, "
        << "but is " << clustering.clusters_by_vertex[v] << " instead";
      return error_msg.str();
    }
  }

  std::vector<std::set<uintE>> actual_clusters_vector(num_clusters);
  for (size_t v = 0; v < clustering.clusters_by_vertex.size(); v++) {
    const i::VertexType& vertex_type{clustering.clusters_by_vertex[v]};
    const i::ClusterMember* const cluster_member{
      std::get_if<i::ClusterMember>(&vertex_type)};
    if (cluster_member != nullptr) {
      for (const uintE cluster_id : cluster_member->clusters) {
        if (cluster_id >= num_clusters) {
          error_msg << "Invalid cluster ID " << cluster_id << " out of "
            << num_clusters << " clusters";
          return error_msg.str();
        }
        actual_clusters_vector[cluster_id].insert(v);
      }
    }
  }
  ClusterList actual_clusters{
    std::make_move_iterator(actual_clusters_vector.begin()),
    std::make_move_iterator(actual_clusters_vector.end())};
  if (actual_clusters != expected_clusters) {
    error_msg << "Clusters don't match";
    return error_msg.str();
  }

  return {};
}

// Wrapper for `CheckClusteringForErrors` that prints out the error message and
// the clustering if there's an error.
bool CheckClustering(
    const i::Clustering& clustering,
    const ClusterList& expected_clusters,
    const VertexList& expected_hubs,
    const VertexList& expected_outliers) {
  const std::optional<std::string> error_msg{
    CheckClusteringForErrors(
        clustering, expected_clusters, expected_hubs, expected_outliers)};
  if (error_msg.has_value()) {
    std::cerr << "Clustering error: " << *error_msg << '\n';
    std::cerr << "Actual clustering: " << clustering << '\n';
    return false;
  } else {
    return true;
  }
}

}  // namespace

TEST(ScanSubroutines, NullGraph) {
  const size_t kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const ii::StructuralSimilarities similarity_table{
    ii::ComputeStructuralSimilarities(&graph)};
  EXPECT_THAT(similarity_table.entries(), IsEmpty());

  const ii::NeighborOrder neighbor_order{
    ii::ComputeNeighborOrder(&graph, similarity_table)};
  EXPECT_THAT(neighbor_order, IsEmpty());

  const auto core_order{ii::ComputeCoreOrder(neighbor_order)};
  EXPECT_THAT(core_order, IsEmpty());
}

TEST(ScanSubroutines, EmptyGraph) {
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const ii::StructuralSimilarities similarity_table{
    ii::ComputeStructuralSimilarities(&graph)};
  EXPECT_THAT(similarity_table.entries(), IsEmpty());

  const ii::NeighborOrder neighbor_order{
    ii::ComputeNeighborOrder(&graph, similarity_table)};
  EXPECT_EQ(neighbor_order.size(), kNumVertices);
  for (const auto& vertex_order : neighbor_order) {
    EXPECT_THAT(vertex_order, IsEmpty());
  }

  const auto core_order{ii::ComputeCoreOrder(neighbor_order)};
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
    {0, 1},
    {1, 2},
    {1, 3},
    {2, 3},
    {2, 4},
    {2, 5},
    {3, 4},
  };
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const ii::StructuralSimilarities similarity_table{
    ii::ComputeStructuralSimilarities(&graph)};
  const auto similarities{similarity_table.entries()};
  EXPECT_THAT(similarities, UnorderedElementsAre(
        std::make_tuple(UndirectedEdge{0, 1}, 2.0 / sqrt(8)   /* .71 */),
        std::make_tuple(UndirectedEdge{1, 2}, 3.0 / sqrt(20)  /* .67 */),
        std::make_tuple(UndirectedEdge{1, 3}, 3.0 / sqrt(16)  /* .75 */),
        std::make_tuple(UndirectedEdge{2, 3}, 4.0 / sqrt(20)  /* .89 */),
        std::make_tuple(UndirectedEdge{2, 4}, 3.0 / sqrt(15)  /* .77 */),
        std::make_tuple(UndirectedEdge{2, 5}, 2.0 / sqrt(10)  /* .63 */),
        std::make_tuple(UndirectedEdge{3, 4}, 3.0 / sqrt(12)  /* .87 */)));

  const ii::NeighborOrder neighbor_order{
    ii::ComputeNeighborOrder(&graph, similarity_table)};
  ASSERT_EQ(neighbor_order.size(), kNumVertices);
  EXPECT_THAT(
      neighbor_order[0],
      ElementsAre(MakeNeighborSimilarity(1, 2.0 / sqrt(8))));
  EXPECT_THAT(
      neighbor_order[1],
      ElementsAre(
        MakeNeighborSimilarity(3, 3.0 / sqrt(16)),
        MakeNeighborSimilarity(0, 2.0 / sqrt(8)),
        MakeNeighborSimilarity(2, 3.0 / sqrt(20))));
  EXPECT_THAT(
      neighbor_order[2],
      ElementsAre(
        MakeNeighborSimilarity(3, 4.0 / sqrt(20)),
        MakeNeighborSimilarity(4, 3.0 / sqrt(15)),
        MakeNeighborSimilarity(1, 3.0 / sqrt(20)),
        MakeNeighborSimilarity(5, 2.0 / sqrt(10))));
  EXPECT_THAT(
      neighbor_order[3],
      ElementsAre(
        MakeNeighborSimilarity(2, 4.0 / sqrt(20)),
        MakeNeighborSimilarity(4, 3.0 / sqrt(12)),
        MakeNeighborSimilarity(1, 3.0 / sqrt(16))));
  EXPECT_THAT(
      neighbor_order[4],
      ElementsAre(
        MakeNeighborSimilarity(3, 3.0 / sqrt(12)),
        MakeNeighborSimilarity(2, 3.0 / sqrt(15))));
  EXPECT_THAT(
      neighbor_order[5],
      ElementsAre(
        MakeNeighborSimilarity(2, 2.0 / sqrt(10))));

  const auto core_order{ii::ComputeCoreOrder(neighbor_order)};
  EXPECT_EQ(core_order.size(), 6);
  EXPECT_THAT(core_order[0], IsEmpty());
  EXPECT_THAT(core_order[1], IsEmpty());
  ASSERT_EQ(core_order[2].size(), 6);
  EXPECT_THAT(
      core_order[2].slice(0, 2),
      UnorderedElementsAre(
        MakeCoreThreshold(2, 4.0 / sqrt(20)),
        MakeCoreThreshold(3, 4.0 / sqrt(20))));
  EXPECT_THAT(
      core_order[2].slice(2, core_order[2].size()),
      ElementsAre(
        MakeCoreThreshold(4, 3.0 / sqrt(12)),
        MakeCoreThreshold(1, 3.0 / sqrt(16)),
        MakeCoreThreshold(0, 2.0 / sqrt(8)),
        MakeCoreThreshold(5, 2.0 / sqrt(10))));
  ASSERT_EQ(core_order[3].size(), 4);
  EXPECT_EQ(core_order[3][0], MakeCoreThreshold(3, 3.0 / sqrt(12)));
  EXPECT_THAT(
      core_order[3].slice(1, 3),
      UnorderedElementsAre(
        MakeCoreThreshold(2, 3.0 / sqrt(15)),
        MakeCoreThreshold(4, 3.0 / sqrt(15))));
  EXPECT_THAT(
      core_order[3].slice(3, core_order[3].size()),
      ElementsAre(MakeCoreThreshold(1, 2.0 / sqrt(8))));
  ASSERT_EQ(core_order[4].size(), 3);
  EXPECT_EQ(core_order[4][0], MakeCoreThreshold(3, 3.0 / sqrt(16)));
  EXPECT_THAT(
      core_order[4].slice(1, core_order[4].size()),
      UnorderedElementsAre(
        MakeCoreThreshold(1, 3.0 / sqrt(20)),
        MakeCoreThreshold(2, 3.0 / sqrt(20))));
  EXPECT_THAT(core_order[5], ElementsAre(MakeCoreThreshold(2, 2.0 / sqrt(10))));
}

TEST(ScanSubroutines, DisconnectedGraph) {
  // Graph diagram:
  //     0 -- 1    2    3 -- 4 -- 5

  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{
    {0, 1},
    {3, 4},
    {4, 5},
  };
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const ii::StructuralSimilarities similarity_table{
    ii::ComputeStructuralSimilarities(&graph)};
  const auto similarities{similarity_table.entries()};
  EXPECT_THAT(similarities, UnorderedElementsAre(
        std::make_tuple(UndirectedEdge{0, 1}, 1.0),
        std::make_tuple(UndirectedEdge{3, 4}, 2.0 / sqrt(6)  /* .82 */),
        std::make_tuple(UndirectedEdge{4, 5}, 2.0 / sqrt(6)  /* .82 */)));

  const ii::NeighborOrder neighbor_order{
    ii::ComputeNeighborOrder(&graph, similarity_table)};
  ASSERT_EQ(neighbor_order.size(), kNumVertices);
  EXPECT_THAT(neighbor_order[0], ElementsAre(MakeNeighborSimilarity(1, 1.0)));
  EXPECT_THAT(neighbor_order[1], ElementsAre(MakeNeighborSimilarity(0, 1.0)));
  EXPECT_THAT(neighbor_order[2], IsEmpty());
  EXPECT_THAT(
      neighbor_order[3],
      ElementsAre(MakeNeighborSimilarity(4, 2.0 / sqrt(6))));
  EXPECT_THAT(
      neighbor_order[4],
      UnorderedElementsAre(
        MakeNeighborSimilarity(3, 2.0 / sqrt(6)),
        MakeNeighborSimilarity(5, 2.0 / sqrt(6))));
  EXPECT_THAT(
      neighbor_order[5],
      ElementsAre(MakeNeighborSimilarity(4, 2.0 / sqrt(6))));

  const auto core_order{ii::ComputeCoreOrder(neighbor_order)};
  EXPECT_EQ(core_order.size(), 4);
  EXPECT_THAT(core_order[0], IsEmpty());
  EXPECT_THAT(core_order[1], IsEmpty());
  ASSERT_EQ(core_order[2].size(), 5);
  EXPECT_THAT(
      core_order[2].slice(0, 2),
      UnorderedElementsAre(
        MakeCoreThreshold(0, 1.0),
        MakeCoreThreshold(1, 1.0)));
  EXPECT_THAT(
      core_order[2].slice(2, core_order[2].size()),
      UnorderedElementsAre(
        MakeCoreThreshold(3, 2.0 / sqrt(6)),
        MakeCoreThreshold(4, 2.0 / sqrt(6)),
        MakeCoreThreshold(5, 2.0 / sqrt(6))));
  EXPECT_THAT(core_order[3], ElementsAre(MakeCoreThreshold(4, 2.0 / sqrt(6))));
}

TEST(Cluster, NullGraph) {
  const size_t kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const i::Index index{&graph};
  constexpr float kMu{2};
  constexpr float kEpsilon{0.5};
  i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

  EXPECT_EQ(clustering.num_clusters, 0);
  EXPECT_THAT(clustering.clusters_by_vertex, IsEmpty());
}

TEST(Cluster, EmptyGraph) {
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{gt::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const i::Index index{&graph};
  constexpr float kMu{2};
  constexpr float kEpsilon{0.5};
  i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

  EXPECT_EQ(clustering.num_clusters, 0);
  ASSERT_THAT(clustering.clusters_by_vertex.size(), kNumVertices);
  for (const auto& vertex_type : clustering.clusters_by_vertex) {
    EXPECT_TRUE(std::holds_alternative<i::Outlier>(vertex_type));
  }
}

TEST(Cluster, BasicUsage) {
  // Graph diagram with similarity scores labeled:
  //       .71    .67  .63
  //     0 --- 1 ---- 2 -- 5
  //           |    / |
  //       .75 |  .89 | .77
  //           | /    |
  //           3 ---- 4
  //             .87
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{
    {0, 1},
    {1, 2},
    {1, 3},
    {2, 3},
    {2, 4},
    {2, 5},
    {3, 4},
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
    EXPECT_TRUE(CheckClustering(
        clustering, kExpectedClusters, kExpectedHubs, kExpectedOutliers));
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.7};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{0, 1, 2, 3, 4}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{5};
    EXPECT_TRUE(CheckClustering(
        clustering, kExpectedClusters, kExpectedHubs, kExpectedOutliers));
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.73};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{1, 2, 3, 4}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 5};
    EXPECT_TRUE(CheckClustering(
        clustering, kExpectedClusters, kExpectedHubs, kExpectedOutliers));
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.88};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{2, 3}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 1, 4, 5};
    EXPECT_TRUE(CheckClustering(
        clustering, kExpectedClusters, kExpectedHubs, kExpectedOutliers));
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.95};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 1, 2, 3, 4, 5};
    EXPECT_TRUE(CheckClustering(
        clustering, kExpectedClusters, kExpectedHubs, kExpectedOutliers));
  }
  {
    constexpr uint64_t kMu{4};
    constexpr float kEpsilon{0.7};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{1, 2, 3, 4}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 5};
    EXPECT_TRUE(CheckClustering(
        clustering, kExpectedClusters, kExpectedHubs, kExpectedOutliers));
  }
}

TEST(Cluster, DisconnectedGraph) {
  // Graph diagram with similarity scores labeled:
  //       1.0            .82  .82
  //     0 -- 1    2    3 -- 4 -- 5

  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{
    {0, 1},
    {3, 4},
    {4, 5},
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
    EXPECT_TRUE(CheckClustering(
        clustering, kExpectedClusters, kExpectedHubs, kExpectedOutliers));
  }
  {
    constexpr uint64_t kMu{2};
    constexpr float kEpsilon{0.8};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{0, 1}, {3, 4, 5}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{2};
    EXPECT_TRUE(CheckClustering(
        clustering, kExpectedClusters, kExpectedHubs, kExpectedOutliers));
  }
  {
    constexpr uint64_t kMu{3};
    constexpr float kEpsilon{0.8};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{{3, 4, 5}};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 1, 2};
    EXPECT_TRUE(CheckClustering(
        clustering, kExpectedClusters, kExpectedHubs, kExpectedOutliers));
  }
  {
    constexpr uint64_t kMu{4};
    constexpr float kEpsilon{0.8};
    const i::Clustering clustering{index.Cluster(kMu, kEpsilon)};

    const ClusterList kExpectedClusters{};
    const VertexList kExpectedHubs{};
    const VertexList kExpectedOutliers{0, 1, 2, 3, 4, 5};
    EXPECT_TRUE(CheckClustering(
        clustering, kExpectedClusters, kExpectedHubs, kExpectedOutliers));
  }
}

TEST(Cluster, TwoClusterGraph) {
  // Graph diagram with similarity scores labeled:
  //           2                   6
  //     .87 /   \ .87       .87 /   \ .87
  //        /     \             /     \
  // 0 --- 1 ----- 3 --- 4 --- 5 ----- 7 --- 8
  //   .71    .75    .58   .58    .75   .71
  const size_t kNumVertices{9};
  const std::unordered_set<UndirectedEdge> kEdges{
    {0, 1},
    {1, 2},
    {1, 3},
    {2, 3},
    {3, 4},
    {4, 5},
    {5, 6},
    {5, 7},
    {6, 7},
    {7, 8},
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
    EXPECT_TRUE(CheckClustering(
        clustering, kExpectedClusters, kExpectedHubs, kExpectedOutliers));
  }
}
