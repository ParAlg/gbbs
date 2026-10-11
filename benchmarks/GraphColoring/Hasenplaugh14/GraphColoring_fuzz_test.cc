#include <algorithm>
#include <cstdint>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "fuzztest/fuzztest.h"
#include "benchmarks/GraphColoring/Hasenplaugh14/GraphColoring.h"
#include "gbbs/bridge.h"
#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "gbbs/vertex.h"

namespace gbbs {
namespace {

using Graph = symmetric_graph<symmetric_vertex, empty>;

// Constructs a symmetric graph from fuzz-generated edges. The fuzz engine
// provides up to 500 arbitrary (u, v) pairs with endpoints in [0, 99]. The
// graph is constructed as follows:
//   - The number of vertices is max(endpoints) + 1.
//   - Self-loops (u == v) are dropped.
//   - Each remaining edge is added in both directions to ensure symmetry.
// Returns std::nullopt if the resulting graph would have no edges.
std::optional<Graph> BuildGraphFromEdges(
    const std::vector<std::pair<uint32_t, uint32_t>>& edges) {
  if (edges.empty()) {
    return std::nullopt;
  }

  uint32_t max_vertex = 0;
  for (const auto& [u, v] : edges) {
    max_vertex = std::max(max_vertex, std::max(u, v));
  }
  uintE n = static_cast<uintE>(max_vertex) + 1;

  std::vector<std::tuple<uintE, uintE, empty>> edge_list;
  for (const auto& [u_raw, v_raw] : edges) {
    uintE u = static_cast<uintE>(u_raw);
    uintE v = static_cast<uintE>(v_raw);
    if (u != v) {
      edge_list.push_back(std::make_tuple(u, v, empty{}));
      edge_list.push_back(std::make_tuple(v, u, empty{}));
    }
  }

  if (edge_list.empty()) {
    return std::nullopt;
  }

  auto symmetric_edges = sequence<std::tuple<uintE, uintE, empty>>(
      edge_list.begin(), edge_list.end());
  return Graph::from_edges(symmetric_edges, n);
}

class GraphColoringFuzzTest {
 public:
  void ColoringIsValidForArbitraryEdges(
      const std::vector<std::pair<uint32_t, uint32_t>>& edges, bool use_lf) {
    auto graph = BuildGraphFromEdges(edges);
    if (!graph.has_value()) return;
    auto colors = Coloring(*graph, use_lf);
    EXPECT_EQ(coloring::CountColoringViolations(*graph, colors), 0);
  }

  void ColoringWithSeedIsValidForArbitraryEdges(
      const std::vector<std::pair<uint32_t, uint32_t>>& edges, bool use_lf,
      int32_t seed) {
    auto graph = BuildGraphFromEdges(edges);
    if (!graph.has_value()) return;
    auto colors = coloring::ColoringWithSeed(*graph, use_lf, seed);
    EXPECT_EQ(coloring::CountColoringViolations(*graph, colors), 0);
  }
};

FUZZ_TEST_F(GraphColoringFuzzTest, ColoringIsValidForArbitraryEdges)
    .WithDomains(
        fuzztest::VectorOf(fuzztest::PairOf(fuzztest::InRange<uint32_t>(0, 99),
                                            fuzztest::InRange<uint32_t>(0, 99)))
            .WithMaxSize(500),
        fuzztest::Arbitrary<bool>());

FUZZ_TEST_F(GraphColoringFuzzTest, ColoringWithSeedIsValidForArbitraryEdges)
    .WithDomains(
        fuzztest::VectorOf(fuzztest::PairOf(fuzztest::InRange<uint32_t>(0, 99),
                                            fuzztest::InRange<uint32_t>(0, 99)))
            .WithMaxSize(500),
        fuzztest::Arbitrary<bool>(), fuzztest::Arbitrary<int32_t>());

}  // namespace
}  // namespace gbbs
