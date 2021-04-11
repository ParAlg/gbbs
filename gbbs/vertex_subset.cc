#include "vertex_subset.h"

namespace gbbs {
void add_to_vsubset(vertexSubset& vs, uintE* new_verts, uintE num_new_verts) {
  if (vs.isDense) {
    parallel_for(0, num_new_verts, [&] (size_t i)
                    { vs.d[new_verts[i]] = true; });
    vs.m += num_new_verts;
  } else {
    const size_t vs_size = vs.numNonzeros();
    const size_t new_size = num_new_verts + vs_size;
    auto all_verts = sequence<uintE>(new_size);
    par_for(0, new_size, kDefaultGranularity, [&] (size_t i)
                    {
                      if (i < vs_size) {
                        all_verts[i] = vs.s[i];
                      } else {
                        all_verts[i] = new_verts[i - vs_size];
                      }
                    });
    auto old_s = std::move(vs.s);
    vs.s = all_verts;
    vs.m = new_size;
  }
}
}  // namespace gbbs
