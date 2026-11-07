
#include "benchmarks/CliqueCounting/Clique.h"
#include "gtest/gtest.h"

TEST(Peel, ZeroContactDoesNotAppendBeforePositiveContribution) {
  using namespace gbbs;
  using Graph = symmetric_graph<symmetric_vertex,empty>;
  sequence<Graph::edge> edges{{0,1,{}},{0,2,{}},{1,2,{}}};
  auto graph=Graph::from_edges(edges,3);
  sequence<uintE> rank{0,2,1};
  auto orient=[&](uintE u,uintE v,empty) { return rank[u]<rank[v]; };
  auto directed=filterGraph(graph,orient);
  sequence<size_t> counts(3,0);
  sequence<char> active(3,0);
  // The spare entry detects duplicate appends without invoking an OOB write.
  sequence<uintE> updates(3,UINT_E_MAX);
  auto count=triUpdate_serial(graph,directed,[](size_t) {return uintE{0};},1,
      active.begin(),rank,counts,false,[](auto&,size_t,uintE) {},updates);
  EXPECT_EQ(count,2u);
  EXPECT_EQ(counts[1],1u);
  EXPECT_EQ(counts[2],1u);
  EXPECT_EQ(updates[2],UINT_E_MAX);
}
