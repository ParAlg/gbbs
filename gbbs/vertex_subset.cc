#include "vertex_subset.h"

namespace gbbs {

// void add_to_vsubset(vertexSubset& vs, uintE* new_verts, uintE num_new_verts) {
//   // TODO: update
//   assert(false);
//   exit(-1);
// //  if (vs.isDense) {
// //    parallel_for(0, num_new_verts, [&] (size_t i)
// //                    { vs.d[new_verts[i]] = true; });
// //    vs.m += num_new_verts;
// //  } else {
// //    const size_t vs_size = vs.numNonzeros();
// //    const size_t new_size = num_new_verts + vs_size;
// //    uintE* all_verts = pbbslib::new_array_no_init<uintE>(new_size);
// //    par_for(0, new_size, gbbs::kSequentialForThreshold, [&] (size_t i)
// //                    {
// //                      if (i < vs_size) {
// //                        all_verts[i] = vs.s[i];
// //                      } else {
// //                        all_verts[i] = new_verts[i - vs_size];
// //                      }
// //                    });
// //    uintE* old_s = vs.s;
// //    vs.s = all_verts;
// //    vs.m = new_size;
// //    if (old_s) {
// //      pbbslib::free_array(old_s);
// //    }
// //  }
// }

}  // namespace gbbs
