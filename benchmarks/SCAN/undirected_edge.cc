#define NOTMAIN

#include "benchmarks/SCAN/undirected_edge.h"

#include <algorithm>

#include "pbbslib/utilities.h"

UndirectedEdge::UndirectedEdge(uintE u, uintE v)
  : edge{std::minmax(u, v)} {}

bool UndirectedEdge::operator==(const UndirectedEdge& other) {
  return edge == other.edge;
}

uint64_t HashUndirectedEdge(const UndirectedEdge& edge) {
  const uint64_t hash{pbbslib::hash64_2(edge.edge.first)};
  return hash ^ (pbbslib::hash64_2(edge.edge.second)
      + 0x9e3779b97f4a7c15 + (hash << 6) + (hash >> 2));
}
