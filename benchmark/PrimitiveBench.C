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
// numactl -i all ./PrimitiveBench -src 10012 -s -m -rounds 3 twitter_SJ
// flags:
//   required:
//     -src: the source to compute the PrimitiveBench from
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#include "ligra.h"

template <class G>
double PrimitiveBench_runner(G& GA, commandLine P) {
  using W = typename G::weight_type;
  timer t; t.start();

  auto ins = pbbs::sequence<uintE>(GA.n, (uintE)1);
  auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
    return ins[v];
  };
  auto outs = pbbs::sequence<uintE>(GA.n);
  auto red_mon = pbbs::addm<size_t>();
  parallel_for(0, GA.n, [&] (size_t i) {
    outs[i] = GA.V[i].template reduceOutNgh<size_t>(i, map_f, red_mon);
  });

  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

generate_main(PrimitiveBench_runner, false);
