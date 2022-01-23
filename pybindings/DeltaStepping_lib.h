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

using asymmetric_unweighted_graph =
    asymmetric_graph<asymmetric_vertex, gbbs::empty>;
using asymmetric_uint32_graph = asymmetric_graph<asymmetric_vertex, uint32_t>;
using asymmetric_float_graph = asymmetric_graph<asymmetric_vertex, float>;
using asymmetric_double_graph = asymmetric_graph<asymmetric_vertex, double>;

using asymmetric_unweighted_cgraph =
    asymmetric_graph<cav_bytepd_amortized, gbbs::empty>;
using asymmetric_uint32_cgraph =
    asymmetric_graph<cav_bytepd_amortized, uint32_t>;
using asymmetric_float_cgraph = asymmetric_graph<cav_bytepd_amortized, float>;
using asymmetric_double_cgraph = asymmetric_graph<cav_bytepd_amortized, double>;

sequence<uintE> DeltaStepping(symmetric_unweighted_graph& G, uintE src,
                              uintE delta);
sequence<uintE> DeltaStepping(symmetric_uint32_graph& G, uintE src,
                              uintE delta);
sequence<float> DeltaStepping(symmetric_float_graph& G, uintE src, float delta);
sequence<double> DeltaStepping(symmetric_double_graph& G, uintE src,
                               double delta);

sequence<uintE> DeltaStepping(symmetric_unweighted_cgraph& G, uintE src,
                              uintE delta);
sequence<uintE> DeltaStepping(symmetric_uint32_cgraph& G, uintE src,
                              uintE delta);
sequence<float> DeltaStepping(symmetric_float_cgraph& G, uintE src,
                              float delta);
sequence<double> DeltaStepping(symmetric_double_cgraph& G, uintE src,
                               double delta);

sequence<uintE> DeltaStepping(asymmetric_unweighted_graph& G, uintE src,
                              uintE delta);
sequence<uintE> DeltaStepping(asymmetric_uint32_graph& G, uintE src,
                              uintE delta);
sequence<float> DeltaStepping(asymmetric_float_graph& G, uintE src,
                              float delta);
sequence<double> DeltaStepping(asymmetric_double_graph& G, uintE src,
                               double delta);

sequence<uintE> DeltaStepping(asymmetric_unweighted_cgraph& G, uintE src,
                              uintE delta);
sequence<uintE> DeltaStepping(asymmetric_uint32_cgraph& G, uintE src,
                              uintE delta);
sequence<float> DeltaStepping(asymmetric_float_cgraph& G, uintE src,
                              float delta);
sequence<double> DeltaStepping(asymmetric_double_cgraph& G, uintE src,
                               double delta);

}  // namespace compiled
}  // namespace gbbs
