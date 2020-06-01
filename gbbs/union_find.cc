#include "union_find.h"

UnionFind::UnionFind(size_t _n) : n(_n) {
  parents = pbbslib::new_array_no_init<intT>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { parents[i] = -1; });
}

intT UnionFind::find(int32_t i) {
  if (parents[i] < 0) return i;
  intT j = parents[i];
  if (parents[j] < 0) return j;
  do
    j = parents[j];
  while (parents[j] >= 0);
  intT tmp;
  while ((tmp = parents[i]) != j) {
    parents[i] = j;
    i = tmp;
  }
  return j;
}

void UnionFind::link(intT u, intT v) { parents[u] = v; }

void UnionFind::clear() { pbbslib::free_array(parents); }

