
#include "benchmarks/CliqueCounting/Clique.h"
#include "gtest/gtest.h"

TEST(ApproxPeel, LaterRoundsKeepCliqueDensityUnits) {
  using namespace gbbs;
  using Graph = symmetric_graph<symmetric_vertex,empty>;
  sequence<Graph::edge> edges{{0,1,{}},{0,2,{}},{0,3,{}},
                            {1,2,{}},{1,3,{}},{2,3,{}}};
  auto graph=Graph::from_edges(edges,5);
  sequence<uintE> rank{0,1,2,3,4};
  auto orient=[&](uintE u,uintE v,empty) {return rank[u]<rank[v];};
  auto directed=filterGraph(graph,orient);
  sequence<size_t> counts{3,3,3,3,0};
  // K4 has four triangles; no induced subgraph has triangle density above 1.
  double density=ApproxPeel(graph,directed,2,counts.begin(),4,true,rank,0.1);
  EXPECT_DOUBLE_EQ(density,1.0);
}
