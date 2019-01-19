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
// numactl -i all ./SetCover -nb 512 -rounds 1 -s -c clueweb_sym.bytepda
// flags:
//   optional:
//     -s : indicate that the graph is symmetric
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -nb : the number of buckets to use in the bucketing implementation

#include "SetCover.h"
#include "ligra.h"

template <class vertex>
void SetCover_runner(graph<vertex>& GA, commandLine P) {
  timer t;
  t.start();
  size_t num_buckets = P.getOptionLongValue("-nb", 128);
  std::cout << "here, nb = " << num_buckets << "\n";
  auto cover = SetCover(GA, num_buckets);
  cover.del();
  t.stop();
  t.reportTotal("SetCover Time");

  // Set-cover mutates the underlying graph (unless it is copied, which
  // we don't do to prevent memory issues), so we make sure the algorithm is run
  // exactly once.
  exit(0);
}

generate_main(SetCover_runner, true)
