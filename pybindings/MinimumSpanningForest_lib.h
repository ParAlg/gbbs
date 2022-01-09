#pragma once

#include "gbbs/gbbs.h"

namespace gbbs {
namespace compiled {

using symmetric_unweighted_graph =
    symmetric_graph<symmetric_vertex, gbbs::empty>;
using symmetric_uint32_graph = symmetric_graph<symmetric_vertex, uint32_t>;
using symmetric_float_graph = symmetric_graph<symmetric_vertex, float>;
using symmetric_double_graph = symmetric_graph<symmetric_vertex, double>;

using symmetric_unweighted_cgraph =
    symmetric_graph<csv_bytepd_amortized, gbbs::empty>;
using symmetric_uint32_cgraph = symmetric_graph<csv_bytepd_amortized, uint32_t>;
using symmetric_float_cgraph = symmetric_graph<csv_bytepd_amortized, float>;
using symmetric_double_cgraph = symmetric_graph<csv_bytepd_amortized, double>;

using edge = std::pair<uintE, uintE>;
using uint_edge = std::tuple<uintE, uintE, uintE>;
using float_edge = std::tuple<uintE, uintE, float>;
using double_edge = std::tuple<uintE, uintE, double>;

sequence<edge> MinimumSpanningForest(symmetric_unweighted_graph& G);
sequence<uint_edge> MinimumSpanningForest(symmetric_uint32_graph& G);
sequence<float_edge> MinimumSpanningForest(symmetric_float_graph& G);
sequence<double_edge> MinimumSpanningForest(symmetric_double_graph& G);

sequence<edge> MinimumSpanningForest(symmetric_unweighted_cgraph& G);
sequence<uint_edge> MinimumSpanningForest(symmetric_uint32_cgraph& G);
sequence<float_edge> MinimumSpanningForest(symmetric_float_cgraph& G);
sequence<double_edge> MinimumSpanningForest(symmetric_double_cgraph& G);

}  // namespace compiled
}  // namespace gbbs
