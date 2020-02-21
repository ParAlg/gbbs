#pragma once

#include <utility>

#include "ligra/macros.h"

class UndirectedEdge;

namespace std {

template <>
struct hash<UndirectedEdge> {
  size_t operator()(const UndirectedEdge& edge) const;
};

}  // namespace std

// Represents an unweighted, undirected edge in a graph.
// This has the property that `UndirectedEdge({u, v}) == Undirected({v, u})`, so
// it can be used for, e.g., storing and retrieving edges from a hash table
// without worrying about the direction of the edge.
class UndirectedEdge {
 public:
  UndirectedEdge(uintE u, uintE v);
  explicit UndirectedEdge(const std::pair<uintE, uintE>& edge);

  bool operator==(const UndirectedEdge& other) const;
  bool operator!=(const UndirectedEdge& other) const;
  bool operator<(const UndirectedEdge& other) const;

  // Returns the two vertices that the edge connects.
  const std::pair<uintE, uintE>& endpoints() const;

 private:
  friend size_t
  std::hash<UndirectedEdge>::operator()(const UndirectedEdge&) const;

  std::pair<uintE, uintE> edge_;
};
