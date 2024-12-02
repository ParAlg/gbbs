#include "gbbs/helpers/directed_edge.h"

#include <algorithm>

namespace std {

size_t hash<gbbs::DirectedEdge>::operator()(
    const gbbs::DirectedEdge& edge) const {
  return gbbs::hash_combine(parlay::hash64_2(edge.edge_.first),
                            parlay::hash64_2(edge.edge_.second));
}

}  // namespace std

namespace gbbs {
DirectedEdge::DirectedEdge(const uintE u, const uintE v)
    : edge_{u, v} {}

DirectedEdge::DirectedEdge(const std::pair<uintE, uintE>& edge)
    : edge_{edge.first, edge.second} {}

bool DirectedEdge::operator==(const DirectedEdge& other) const {
  return edge_ == other.edge_;
}

bool DirectedEdge::operator!=(const DirectedEdge& other) const {
  return !(*this == other);
}

const std::pair<uintE, uintE>& DirectedEdge::endpoints() const {
  return edge_;
}

std::ostream& operator<<(std::ostream& os, const DirectedEdge& edge) {
  const std::pair<uintE, uintE>& endpoints{edge.endpoints()};
  os << "{" << endpoints.first << ", " << endpoints.second << "}";
  return os;
}
}  // namespace gbbs
