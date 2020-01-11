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

// Represents an undirected edge in a graph.
class UndirectedEdge {
 public:
  UndirectedEdge(uintE u, uintE v);
  explicit UndirectedEdge(std::pair<uintE, uintE> edge);

  bool operator==(const UndirectedEdge& other) const;
  bool operator!=(const UndirectedEdge& other) const;

  uintE from() const;
  uintE to() const;

 private:
  friend size_t std::hash<UndirectedEdge>::operator()(
      const UndirectedEdge&) const;

  std::pair<uintE, uintE> edge_;
};
