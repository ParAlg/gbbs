// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "DeltaStepping_lib.h"

#include "benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaStepping.h"

#include "gbbs/gbbs.h"

namespace gbbs {
namespace compiled {

sequence<uintE> DeltaStepping(symmetric_unweighted_graph& G, uintE src,
                              uintE delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<uintE> DeltaStepping(symmetric_uint32_graph& G, uintE src,
                              uintE delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<float> DeltaStepping(symmetric_float_graph& G, uintE src,
                              float delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<double> DeltaStepping(symmetric_double_graph& G, uintE src,
                               double delta) {
  return gbbs::DeltaStepping(G, src, delta);
}

sequence<uintE> DeltaStepping(symmetric_unweighted_cgraph& G, uintE src,
                              uintE delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<uintE> DeltaStepping(symmetric_uint32_cgraph& G, uintE src,
                              uintE delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<float> DeltaStepping(symmetric_float_cgraph& G, uintE src,
                              float delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<double> DeltaStepping(symmetric_double_cgraph& G, uintE src,
                               double delta) {
  return gbbs::DeltaStepping(G, src, delta);
}

sequence<uintE> DeltaStepping(asymmetric_unweighted_graph& G, uintE src,
                              uintE delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<uintE> DeltaStepping(asymmetric_uint32_graph& G, uintE src,
                              uintE delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<float> DeltaStepping(asymmetric_float_graph& G, uintE src,
                              float delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<double> DeltaStepping(asymmetric_double_graph& G, uintE src,
                               double delta) {
  return gbbs::DeltaStepping(G, src, delta);
}

sequence<uintE> DeltaStepping(asymmetric_unweighted_cgraph& G, uintE src,
                              uintE delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<uintE> DeltaStepping(asymmetric_uint32_cgraph& G, uintE src,
                              uintE delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<float> DeltaStepping(asymmetric_float_cgraph& G, uintE src,
                              float delta) {
  return gbbs::DeltaStepping(G, src, delta);
}
sequence<double> DeltaStepping(asymmetric_double_cgraph& G, uintE src,
                               double delta) {
  return gbbs::DeltaStepping(G, src, delta);
}

}  // namespace compiled
}  // namespace gbbs
