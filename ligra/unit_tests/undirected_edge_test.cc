#include "ligra/undirected_edge.h"

#include <gtest/gtest.h>

TEST(UndirectedEdge, Equality) {
  EXPECT_EQ((UndirectedEdge{1, 2}), (UndirectedEdge{1, 2}));
  EXPECT_EQ((UndirectedEdge{1, 2}), (UndirectedEdge{2, 1}));
  EXPECT_NE((UndirectedEdge{1, 2}), (UndirectedEdge{3, 1}));
}
