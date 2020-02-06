#include "benchmarks/SCAN/scan.h"

#include <cmath>
#include <unordered_set>
#include <utility>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ligra/graph.h"
#include "ligra/vertex.h"
#include "pbbslib/seq.h"

using ::testing::ElementsAre;
using ::testing::IsEmpty;
using ::testing::UnorderedElementsAre;

namespace si = scan::internal;

namespace {

// Make an undirected graph from a list of edges.
symmetric_graph<symmetric_vertex, pbbslib::empty> MakeGraph(
    const uintE num_vertices,
    const std::unordered_set<UndirectedEdge>& edges) {
  pbbs::sequence<std::tuple<uintE, uintE, pbbslib::empty>> edge_sequence(
      edges.size() * 2);
  auto edges_it{edges.cbegin()};
  for (size_t i = 0; i < edges.size(); i++) {
    edge_sequence[2 * i] =
      std::make_tuple(edges_it->from(), edges_it->to(), pbbs::empty{});
    edge_sequence[2 * i + 1] =
      std::make_tuple(edges_it->to(), edges_it->from(), pbbs::empty{});
    ++edges_it;
  }
  return sym_graph_from_edges(edge_sequence, num_vertices);
}

si::NeighborSimilarity
MakeNeighborSimilarity(const uintE neighbor, const double similarity) {
  return {.neighbor = neighbor, .similarity = static_cast<float>(similarity)};
}

si::CoreThreshold
MakeCoreThreshold(const uintE vertex_id, const double threshold) {
  return {.vertex_id = vertex_id, .threshold = static_cast<float>(threshold)};
}

}  // namespace

TEST(ScanSubroutines, NullGraph) {
  const size_t kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{MakeGraph(kNumVertices, kEdges)};

  const si::StructuralSimilarities similarity_table{
    si::ComputeStructuralSimilarities(&graph)};
  EXPECT_THAT(similarity_table.entries(), IsEmpty());

  const si::NeighborOrder neighbor_order{
    si::ComputeNeighborOrder(&graph, similarity_table)};
  EXPECT_THAT(neighbor_order, IsEmpty());

  const si::CoreOrder core_order{si::ComputeCoreOrder(neighbor_order)};
  EXPECT_THAT(core_order, IsEmpty());
}

TEST(ScanSubroutines, EmptyGraph) {
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{MakeGraph(kNumVertices, kEdges)};

  const si::StructuralSimilarities similarity_table{
    si::ComputeStructuralSimilarities(&graph)};
  EXPECT_THAT(similarity_table.entries(), IsEmpty());

  const si::NeighborOrder neighbor_order{
    si::ComputeNeighborOrder(&graph, similarity_table)};
  EXPECT_EQ(neighbor_order.size(), kNumVertices);
  for (const auto& vertex_order : neighbor_order) {
    EXPECT_THAT(vertex_order, IsEmpty());
  }

  const si::CoreOrder core_order{si::ComputeCoreOrder(neighbor_order)};
  ASSERT_EQ(core_order.size(), 2);
  EXPECT_THAT(core_order[0], IsEmpty());
  EXPECT_THAT(core_order[1], IsEmpty());
}

TEST(ScanSubroutines, BasicUsage) {
  // Graph diagram:
  //
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
  auto graph{MakeGraph(kNumVertices, kEdges)};

  const si::StructuralSimilarities similarity_table{
    si::ComputeStructuralSimilarities(&graph)};
  const auto similarities{similarity_table.entries()};
  EXPECT_THAT(similarities.slice(), UnorderedElementsAre(
        std::make_tuple(UndirectedEdge{0, 1}, 2.0 / sqrt(8)   /* .71 */),
        std::make_tuple(UndirectedEdge{1, 2}, 3.0 / sqrt(20)  /* .67 */),
        std::make_tuple(UndirectedEdge{1, 3}, 3.0 / sqrt(16)  /* .75 */),
        std::make_tuple(UndirectedEdge{2, 3}, 4.0 / sqrt(20)  /* .89 */),
        std::make_tuple(UndirectedEdge{2, 4}, 3.0 / sqrt(15)  /* .77 */),
        std::make_tuple(UndirectedEdge{2, 5}, 2.0 / sqrt(10)  /* .63 */),
        std::make_tuple(UndirectedEdge{3, 4}, 3.0 / sqrt(12)  /* .87 */)));

  const si::NeighborOrder neighbor_order{
    si::ComputeNeighborOrder(&graph, similarity_table)};
  ASSERT_EQ(neighbor_order.size(), kNumVertices);
  EXPECT_THAT(
      neighbor_order[0].slice(),
      ElementsAre(MakeNeighborSimilarity(1, 2.0 / sqrt(8))));
  EXPECT_THAT(
      neighbor_order[1].slice(),
      ElementsAre(
        MakeNeighborSimilarity(3, 3.0 / sqrt(16)),
        MakeNeighborSimilarity(0, 2.0 / sqrt(8)),
        MakeNeighborSimilarity(2, 3.0 / sqrt(20))));
  EXPECT_THAT(
      neighbor_order[2].slice(),
      ElementsAre(
        MakeNeighborSimilarity(3, 4.0 / sqrt(20)),
        MakeNeighborSimilarity(4, 3.0 / sqrt(15)),
        MakeNeighborSimilarity(1, 3.0 / sqrt(20)),
        MakeNeighborSimilarity(5, 2.0 / sqrt(10))));
  EXPECT_THAT(
      neighbor_order[3].slice(),
      ElementsAre(
        MakeNeighborSimilarity(2, 4.0 / sqrt(20)),
        MakeNeighborSimilarity(4, 3.0 / sqrt(12)),
        MakeNeighborSimilarity(1, 3.0 / sqrt(16))));
  EXPECT_THAT(
      neighbor_order[4].slice(),
      ElementsAre(
        MakeNeighborSimilarity(3, 3.0 / sqrt(12)),
        MakeNeighborSimilarity(2, 3.0 / sqrt(15))));
  EXPECT_THAT(
      neighbor_order[5].slice(),
      ElementsAre(
        MakeNeighborSimilarity(2, 2.0 / sqrt(10))));

  const si::CoreOrder core_order{si::ComputeCoreOrder(neighbor_order)};
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
  EXPECT_THAT(
      core_order[5].slice(),
      ElementsAre(MakeCoreThreshold(2, 2.0 / sqrt(10))));
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
  auto graph{MakeGraph(kNumVertices, kEdges)};

  const si::StructuralSimilarities similarity_table{
    si::ComputeStructuralSimilarities(&graph)};
  const auto similarities{similarity_table.entries()};
  EXPECT_THAT(similarities.slice(), UnorderedElementsAre(
        std::make_tuple(UndirectedEdge{0, 1}, 1.0),
        std::make_tuple(UndirectedEdge{3, 4}, 2.0 / sqrt(6)  /* .82 */),
        std::make_tuple(UndirectedEdge{4, 5}, 2.0 / sqrt(6)  /* .82 */)));

  const si::NeighborOrder neighbor_order{
    si::ComputeNeighborOrder(&graph, similarity_table)};
  ASSERT_EQ(neighbor_order.size(), kNumVertices);
  EXPECT_THAT(
      neighbor_order[0].slice(), ElementsAre(MakeNeighborSimilarity(1, 1.0)));
  EXPECT_THAT(
      neighbor_order[1].slice(), ElementsAre(MakeNeighborSimilarity(0, 1.0)));
  EXPECT_THAT(neighbor_order[2], IsEmpty());
  EXPECT_THAT(
      neighbor_order[3].slice(),
      ElementsAre(MakeNeighborSimilarity(4, 2.0 / sqrt(6))));
  EXPECT_THAT(
      neighbor_order[4].slice(),
      UnorderedElementsAre(
        MakeNeighborSimilarity(3, 2.0 / sqrt(6)),
        MakeNeighborSimilarity(5, 2.0 / sqrt(6))));
  EXPECT_THAT(
      neighbor_order[5].slice(),
      ElementsAre(MakeNeighborSimilarity(4, 2.0 / sqrt(6))));

  const si::CoreOrder core_order{si::ComputeCoreOrder(neighbor_order)};
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
  EXPECT_THAT(
      core_order[3].slice(),
      ElementsAre(MakeCoreThreshold(4, 2.0 / sqrt(6))));
}

TEST(Cluster, NullGraph) {
  const size_t kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{MakeGraph(kNumVertices, kEdges)};

  const scan::ScanIndex index{&graph};
  constexpr float kEpsilon{0.5};
  constexpr float kMu{2};
  scan::Clustering clustering{index.Cluster(kEpsilon, kMu)};

  EXPECT_EQ(clustering.num_clusters, 0);
  EXPECT_THAT(clustering.clusters_by_vertex, IsEmpty());
}

TEST(Cluster, EmptyGraph) {
  const size_t kNumVertices{6};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{MakeGraph(kNumVertices, kEdges)};

  const scan::ScanIndex index{&graph};
  constexpr float kEpsilon{0.5};
  constexpr float kMu{2};
  scan::Clustering clustering{index.Cluster(kEpsilon, kMu)};

  EXPECT_EQ(clustering.num_clusters, 0);
  ASSERT_THAT(clustering.clusters_by_vertex.size(), kNumVertices);
  for (const auto& vertex_type : clustering.clusters_by_vertex) {
    EXPECT_TRUE(std::holds_alternative<scan::Outlier>(vertex_type));
  }
}
