#include "benchmarks/CliqueCounting/Clique.h"

#include <unordered_set>

#include "gbbs/graph.h"
#include "gbbs/macros.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using ::testing::AnyOf;
using ::testing::ElementsAre;

namespace gbbs {

TEST(CliqueCounting, EdgelessGraph) {
  constexpr uintE kNumVertices{10};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  size_t clique_count = Clique(graph,
                               /* k = */ 4,
                               /* order_type = */ 0,
                               /* epsilon = */ 0.01,
                               /* space_type = */ 5,
                               /* label = */ false,
                               /* filter = */ false,
                               /* use_base = */ false,
                               /* recursive_level = */ 0,
                               /* approx_peel = */ false,
                               /* approx_eps = */ 0.01);
  EXPECT_EQ(clique_count, 0);
}

// TEST(CliqueCounting, BasicUsage) {
//  // Graph diagram:
//  //     K5 on (0,1,2,3,4)    5 - 6 - 9
//  //                           \  | / |
//  //                             7 -- 8 -- 10
//  constexpr uintE kNumVertices{11};
//  const std::unordered_set<UndirectedEdge> kEdges{
//    {0, 1}, {0, 2}, {0, 3}, {0, 4},
//    {1, 2}, {1, 3}, {1, 4},
//    {2, 3}, {2, 4},
//    {3, 4},
//    {5, 6}, {5, 7},
//    {6, 7}, {6, 9},
//    {7, 8}, {7, 9},
//    {8, 9}, {8, 10}
//  };
//  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
//  std::cout << "Before running clique" << std::endl;
//
//  auto print = [&] (const uintE& u, const uintE& v, const gbbs::empty& wgh) {
//    std::cout << u << " " << v << std::endl;
//  };
//  // graph.get_vertex(0).out_neighbors().map(print, false);
//  // graph.get_vertex(1).out_neighbors().map(print, false);
//  // graph.get_vertex(2).out_neighbors().map(print, false);
//  // graph.get_vertex(3).out_neighbors().map(print, false);
//  // graph.get_vertex(4).out_neighbors().map(print, false);
//  // graph.get_vertex(5).out_neighbors().map(print, false);
//  // graph.get_vertex(6).out_neighbors().map(print, false);
//  // graph.get_vertex(7).out_neighbors().map(print, false);
//  // graph.get_vertex(8).out_neighbors().map(print, false);
//  // graph.get_vertex(9).out_neighbors().map(print, false);
//  // graph.get_vertex(10).out_neighbors().map(print, false);
//
//  size_t triangle_count = Clique(graph,
//      /* k = */ 3,
//      /* order_type = */ 0,
//      /* epsilon = */ 0.01,
//      /* space_type = */ 5,
//      /* label = */ false,
//      /* filter = */ false,
//      /* use_base = */ false,
//      /* recursive_level = */ 0,
//      /* approx_peel = */ false,
//      /* approx_eps = */ 0.01);
//  EXPECT_EQ(triangle_count, 13);
//
//  size_t four_clique_count = Clique(graph,
//      /* k = */ 4,
//      /* order_type = */ 0,
//      /* epsilon = */ 0.01,
//      /* space_type = */ 5,
//      /* label = */ false,
//      /* filter = */ false,
//      /* use_base = */ false,
//      /* recursive_level = */ 0,
//      /* approx_peel = */ false,
//      /* approx_eps = */ 0.01);
//  EXPECT_EQ(four_clique_count, 5);
//}

}  // namespace gbbs
