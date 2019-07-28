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
#include "packed_graph.h"

template <class G>
void bench_mapreduce(G& GA) {
  using W = typename G::weight_type;

  timer t; t.start();
  auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
    return v; // no fetch/edge
  };
  auto outs = pbbs::sequence<uintE>(GA.n);
  auto red_m = pbbs::addm<size_t>();
  parallel_for(0, GA.n, [&] (size_t i) {
    outs[i] = GA.get_vertex(i).reduceOutNgh(i, map_f, red_m);
  });
  t.stop(); t.reportTotal("bench mapreduce time");
}

template <class G>
void bench_mapreduce_fetch(G& GA) {
  using W = typename G::weight_type;

  auto ins = pbbs::sequence<uintE>(GA.n, (uintE)1);
  timer t; t.start();
  auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
    return ins[v]; // no fetch/edge
  };
  auto outs = pbbs::sequence<uintE>(GA.n);
  auto red_m = pbbs::addm<size_t>();
  parallel_for(0, GA.n, [&] (size_t i) {
    outs[i] = GA.get_vertex(i).reduceOutNgh(i, map_f, red_m);
  });
  t.stop(); t.reportTotal("bench mapreduce_fetch time");
}

template <class G>
void bench_packedgraph(G& GA) {
  using W = typename G::weight_type;

  timer t; t.start();
  auto PG = build_packed_graph(GA);
  PG.del();
  t.stop(); t.reportTotal("bench packedgraph time");
}

template <class G>
double PrimitiveBench_runner(G& GA, commandLine P) {
  size_t test = P.getOptionLongValue("-t", 0);
  switch(test) {
    case 0:
      bench_mapreduce(GA);
      break;
    case 1:
      bench_mapreduce_fetch(GA);
      break;
    case 2:
      bench_packedgraph(GA);
      break;
    default:
      cout << "-t : test_num" << endl;
      exit(0);
  }
}

generate_main(PrimitiveBench_runner, false);
