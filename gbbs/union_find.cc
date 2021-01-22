#include "union_find.h"

namespace gbbs {
UnionFind::UnionFind(size_t _n) : n(_n) {
  parents = parlay::sequence<uintE>(n);
  parallel_for(0, n, [&] (size_t i) { parents[i] = -1; });
}

uintE UnionFind::find(int32_t i) {
  if (parents[i] < 0) return i;
  uintE j = parents[i];
  if (parents[j] < 0) return j;
  do
    j = parents[j];
  while (parents[j] >= 0);
  uintE tmp;
  while ((tmp = parents[i]) != j) {
    parents[i] = j;
    i = tmp;
  }
  return j;
}

void UnionFind::link(uintE u, uintE v) { parents[u] = v; }

}  // namespace gbbs
