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

#include "HAC_lib.h"

#include "benchmarks/Clustering/SeqHAC/HAC_api.h"

#include "gbbs/gbbs.h"

namespace gbbs {
namespace compiled {

sequence<uint_parent> HAC(symmetric_uint32_graph& G, std::string linkage,
                          bool similarity) {
  return gbbs::HAC(G, linkage, similarity);
}
sequence<float_parent> HAC(symmetric_float_graph& G, std::string linkage,
                           bool similarity) {
  return gbbs::HAC(G, linkage, similarity);
}
sequence<double_parent> HAC(symmetric_double_graph& G, std::string linkage,
                            bool similarity) {
  return gbbs::HAC(G, linkage, similarity);
}

sequence<uint_parent> HAC(symmetric_uint32_cgraph& G, std::string linkage,
                          bool similarity) {
  return gbbs::HAC(G, linkage, similarity);
}
sequence<float_parent> HAC(symmetric_float_cgraph& G, std::string linkage,
                           bool similarity) {
  return gbbs::HAC(G, linkage, similarity);
}
sequence<double_parent> HAC(symmetric_double_cgraph& G, std::string linkage,
                            bool similarity) {
  return gbbs::HAC(G, linkage, similarity);
}

}  // namespace compiled
}  // namespace gbbs
