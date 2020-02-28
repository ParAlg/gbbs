// Utility functions for working with graphs that are useful for writing unit
// tests.
#pragma once

#include <unordered_set>

#include "ligra/graph.h"
#include "ligra/macros.h"
#include "ligra/undirected_edge.h"
#include "ligra/vertex.h"
#include "pbbslib/utilities.h"

namespace graph_test {

// Make an undirected, unweighted graph from a list of edges.
symmetric_graph<symmetric_vertex, pbbslib::empty> MakeUnweightedSymmetricGraph(
    const uintE num_vertices,
    const std::unordered_set<UndirectedEdge>& edges);

}  // namespace graph_test
