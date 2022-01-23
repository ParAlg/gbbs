#include "gbbs/helpers/undirected_edge.h"

#include <algorithm>

namespace std {

size_t hash<gbbs::UndirectedEdge>::operator()(
    const gbbs::UndirectedEdge& edge) const {
  return gbbs::hash_combine(parlay::hash64_2(edge.edge_.first),
                            parlay::hash64_2(edge.edge_.second));
}

}  // namespace std

namespace gbbs {
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
  os << "{" << endpoints.first << ", " << endpoints.second << "}";
  return os;
}
}  // namespace gbbs
