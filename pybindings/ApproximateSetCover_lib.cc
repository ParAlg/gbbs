#include "ApproximateSetCover_lib.h"

#include "benchmarks/ApproximateSetCover/MANISBPT11/ApproximateSetCover.h"

#include "gbbs/gbbs.h"

namespace gbbs {
namespace compiled {

sequence<uintE> ApproximateSetCover(symmetric_unweighted_graph& G,
                                    size_t num_buckets) {
  return gbbs::SetCover<symmetric_unweighted_graph>(G, num_buckets);
}

sequence<uintE> ApproximateSetCover(symmetric_uint32_graph& G,
                                    size_t num_buckets) {
  return gbbs::SetCover<symmetric_uint32_graph>(G, num_buckets);
}

sequence<uintE> ApproximateSetCover(symmetric_float_graph& G,
                                    size_t num_buckets) {
  return gbbs::SetCover<symmetric_float_graph>(G, num_buckets);
}

sequence<uintE> ApproximateSetCover(symmetric_double_graph& G,
                                    size_t num_buckets) {
  return gbbs::SetCover<symmetric_double_graph>(G, num_buckets);
}

sequence<uintE> ApproximateSetCover(symmetric_unweighted_cgraph& G,
                                    size_t num_buckets) {
  return gbbs::SetCover<symmetric_unweighted_cgraph>(G, num_buckets);
}

sequence<uintE> ApproximateSetCover(symmetric_uint32_cgraph& G,
                                    size_t num_buckets) {
  return gbbs::SetCover<symmetric_uint32_cgraph>(G, num_buckets);
}

sequence<uintE> ApproximateSetCover(symmetric_float_cgraph& G,
                                    size_t num_buckets) {
  return gbbs::SetCover<symmetric_float_cgraph>(G, num_buckets);
}

sequence<uintE> ApproximateSetCover(symmetric_double_cgraph& G,
                                    size_t num_buckets) {
  return gbbs::SetCover<symmetric_double_cgraph>(G, num_buckets);
}

}  // namespace compiled
}  // namespace gbbs
