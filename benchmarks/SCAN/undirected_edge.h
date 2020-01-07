#pragma once

#include <utility>

#include "ligra/macros.h"

// Represents an undirected edge in a graph.
class UndirectedEdge {
 public:
  UndirectedEdge(uintE u, uintE v);

  bool operator==(const UndirectedEdge& other);

 private:
  friend uint64_t HashUndirectedEdge(const UndirectedEdge& edge);

  std::pair<uintE, uintE> edge;
};

// Hashes an `UndirectedEdge`.
uint64_t HashUndirectedEdge(const UndirectedEdge& edge);
