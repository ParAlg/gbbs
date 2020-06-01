#include "gbbs/undirected_edge.h"

#include <algorithm>

#include "pbbslib/utilities.h"

namespace std {

size_t hash<UndirectedEdge>::operator()(const UndirectedEdge& edge) const {
  return pbbs::hash_combine(
      pbbs::hash64_2(edge.edge_.first), pbbs::hash64_2(edge.edge_.second));
}

}  // namespace std

UndirectedEdge::UndirectedEdge(const uintE u, const uintE v)
  : edge_{std::minmax(u, v)} {}

UndirectedEdge::UndirectedEdge(const std::pair<uintE, uintE>& edge)
  : edge_{std::minmax(edge.first, edge.second)} {}

bool UndirectedEdge::operator==(const UndirectedEdge& other) const {
  return edge_ == other.edge_;
}

bool UndirectedEdge::operator!=(const UndirectedEdge& other) const {
  return !(*this == other);
}

const std::pair<uintE, uintE>& UndirectedEdge::endpoints() const {
  return edge_;
}

std::ostream& operator<<(std::ostream& os, const UndirectedEdge& edge) {
  const std::pair<uintE, uintE>& endpoints{edge.endpoints()};
  os<< "{" << endpoints.first << ", " << endpoints.second << "}";
  return os;
}
