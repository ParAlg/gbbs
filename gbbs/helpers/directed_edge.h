#pragma once

#include <ostream>
#include <utility>

#include "gbbs/macros.h"

namespace gbbs {

class DirectedEdge;

}  // namespace gbbs

namespace std {

template <>
struct hash<gbbs::DirectedEdge> {
  size_t operator()(const gbbs::DirectedEdge& edge) const;
};

}  // namespace std

namespace gbbs {

class DirectedEdge;

// Represents an unweighted, directed edge in a graph. This class is primarily
// used for testing.
class DirectedEdge {
 public:
  DirectedEdge(uintE u, uintE v);
  explicit DirectedEdge(const std::pair<uintE, uintE>& edge);

  bool operator==(const DirectedEdge& other) const;
  bool operator!=(const DirectedEdge& other) const;
  bool operator<(const DirectedEdge& other) const;

  // Returns the two vertices that the edge connects.
  const std::pair<uintE, uintE>& endpoints() const;

 private:
  friend size_t std::hash<DirectedEdge>::operator()(
      const DirectedEdge&) const;

  std::pair<uintE, uintE> edge_;
};

std::ostream& operator<<(std::ostream& os, const DirectedEdge&);

}  // namespace gbbs
