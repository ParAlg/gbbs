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
// numactl -i all ./BellmanFord -src 10012 -s -m -rounds 3 twitter_wgh_SJ
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//     -w: indicate that the graph is weighted
//   optional:
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -specfor : run the speculative_for (union-find) based algorithm from pbbs
//     -largemem : set the sampling thresholds to utilize less memory
//
// Note: in our experiments we set -largemem when running MST on the weighted
// hyperlink2012 graph.

#include "MST.h"
#include "ligra.h"

template <class vertex>
void MST_runner(graph<vertex>& GA, commandLine P) {
  bool spec_for = P.getOption("-specfor");
  bool largemem = P.getOption("-largemem");
  timer mst_t;
  mst_t.start();
  if (spec_for) {
    MST_spec_for::MST(GA);
  } else {
    MST_boruvka::MST(GA, largemem);
  }
  mst_t.stop();
  mst_t.reportTotal("MST time");
  // MST mutates the underlying graph (unless it is copied, which we don't do to
  // prevent memory issues), so we make sure the algorithm is run exactly once.
  exit(0);
}

generate_weighted_main(MST_runner, true);
