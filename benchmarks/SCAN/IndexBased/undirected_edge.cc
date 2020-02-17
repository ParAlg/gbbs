#define NOTMAIN

#include "benchmarks/SCAN/undirected_edge.h"

#include <algorithm>

#include "pbbslib/utilities.h"

namespace std {

size_t hash<UndirectedEdge>::operator()(const UndirectedEdge& edge) const {
  const uint64_t hash{pbbslib::hash64_2(edge.edge_.first)};
  return hash ^ (pbbslib::hash64_2(edge.edge_.second)
      + 0x9e3779b97f4a7c15 + (hash << 6) + (hash >> 2));
}

}  // namespace std

UndirectedEdge::UndirectedEdge(const uintE u, const uintE v)
  : edge_{std::minmax(u, v)} {}

UndirectedEdge::UndirectedEdge(const std::pair<uintE, uintE> edge)
  : edge_{std::minmax(edge.first, edge.second)} {}

bool UndirectedEdge::operator==(const UndirectedEdge& other) const {
  return edge_ == other.edge_;
}

bool UndirectedEdge::operator!=(const UndirectedEdge& other) const {
  return !(*this == other);
}

uintE UndirectedEdge::from() const {
  return edge_.first;
}

uintE UndirectedEdge::to() const {
  return edge_.second;
}
