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

// Usage:
// numactl -i all ./CC -rounds 1 -pack -s -m twitter_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -pack : delete edges from G
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc

#include "CC.h"
#include "ligra.h"

template <class vertex>
void CC_runner(graph<vertex>& GA, commandLine P) {
  auto beta = P.getOptionDoubleValue("-beta", 0.2);
  auto pack = P.getOption("-pack");
  auto sym = P.getOption("-s");
  assert(sym);
  timer t;
  t.start();
  auto components = cc::CC(GA, beta, pack);
  if (P.getOption("-stats")) {
    auto cc_f = [&](size_t i) { return components[i]; };
    auto cc_im =
        pbbslib::make_sequence<uintE>(GA.n, cc_f);
    cc::num_cc(cc_im);
    cc::largest_cc(cc_im);
  }
  t.stop(); t.reportTotal("CC time");
  components.clear();
  if (pack) {
    // packing mutates the graph, packing out all intra-cluster edges, and can
    // only be run once unless the input graph is copied.
    exit(0);
  }
}

generate_main(CC_runner, true);
