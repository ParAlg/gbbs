
#include "benchmarks/SSWidestPath/JulienneDBS17/SSWidestPath.h"
#include "gtest/gtest.h"

TEST(SSWidestPath, FirstReachBeyondBucketWindow) {
  using Graph = gbbs::symmetric_graph<gbbs::symmetric_vertex,gbbs::intE>;
  // A disconnected heavy edge sets the global weight range. Both the bogus
  // previous bucket and the first reachable bucket then map to overflow.
  gbbs::sequence<Graph::edge> edges{{0,1,1},{1,2,1},{3,4,1000}};
  auto graph = Graph::from_edges(edges,5);
  auto widths = gbbs::SSWidestPath(graph,0,4);
  EXPECT_EQ(widths[1],1u);
  EXPECT_EQ(widths[2],1u);
  EXPECT_EQ(widths[3],0u);
  EXPECT_EQ(widths[4],0u);
}
