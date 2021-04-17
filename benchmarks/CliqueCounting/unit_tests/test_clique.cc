#include "benchmarks/CliqueCounting/Clique.h"

#include <unordered_set>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gbbs/graph.h"
#include "gbbs/graph_test_utils.h"
#include "gbbs/macros.h"

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

TEST(CliqueCounting, BasicUsage) {
  // Graph diagram:
  //     K5 on (0,1,2,3,4)    5 - 6 - 9
  //                           \  | / |
  //                             7 -- 8 -- 10
  constexpr uintE kNumVertices{5};
  const std::unordered_set<UndirectedEdge> kEdges{
    UndirectedEdge(0, 1), UndirectedEdge(0, 2), UndirectedEdge(0, 3), UndirectedEdge(0, 4),
    UndirectedEdge(1, 2), UndirectedEdge(1, 3), UndirectedEdge(1, 4),
    UndirectedEdge(2, 3), UndirectedEdge(2, 4),
    UndirectedEdge(3, 4),
    //UndirectedEdge{5, 6}, UndirectedEdge{5, 7},
    //UndirectedEdge{6, 7}, UndirectedEdge{6, 9},
    //UndirectedEdge{7, 8}, UndirectedEdge{7, 9},
    //UndirectedEdge{8, 9}, UndirectedEdge{8, 10}
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};
//  std::cout << "Before running clique" << std::endl;
//
  auto print = [&] (const uintE& u, const uintE& v, const gbbs::empty& wgh) {
    std::cout << u << " " << v << std::endl;
  };
   graph.get_vertex(0).out_neighbors().map(print, false);
   graph.get_vertex(1).out_neighbors().map(print, false);
   graph.get_vertex(2).out_neighbors().map(print, false);
   graph.get_vertex(3).out_neighbors().map(print, false);
   graph.get_vertex(4).out_neighbors().map(print, false);
   std::cout << "Deg 0: " << graph.get_vertex(0).out_degree() << std::endl;
   std::cout << "Deg 1: " << graph.get_vertex(1).out_degree() << std::endl;
   std::cout << "Deg 2: " << graph.get_vertex(2).out_degree() << std::endl;
   std::cout << "Deg 3: " << graph.get_vertex(3).out_degree() << std::endl;
   std::cout << "Deg 4: " << graph.get_vertex(4).out_degree() << std::endl;
   std::cout << "Max deg: " << gbbs::induced_hybrid::get_max_deg(graph) << std::endl;
//  // graph.get_vertex(5).out_neighbors().map(print, false);
//  // graph.get_vertex(6).out_neighbors().map(print, false);
//  // graph.get_vertex(7).out_neighbors().map(print, false);
//  // graph.get_vertex(8).out_neighbors().map(print, false);
//  // graph.get_vertex(9).out_neighbors().map(print, false);
//  // graph.get_vertex(10).out_neighbors().map(print, false);

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
//  EXPECT_EQ(triangle_count, 10);

  size_t four_clique_count = Clique(graph,
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
  EXPECT_EQ(four_clique_count, 5);
}

}  // namespace gbbs
