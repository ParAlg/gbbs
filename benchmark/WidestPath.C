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
// numactl -i all ./WidestPath -src 10012 -s -m -rounds 3 twitter_wgh_SJ
// flags:
//   required:
//     -src: the source to compute shortest path distances from
//     -w: indicate that the graph is weighted
//   optional:
//     -rounds : the number of times to run the algorithm
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -s : indicate that the graph is symmetric

#define WEIGHTED 1

#include "WidestPath.h"
#include "weight_utils.h"

template <class G>
double WidestPath_runner(G& GA, commandLine P) {
  uintE src = P.getOptionLongValue("-src", 0);
  size_t num_buckets = P.getOptionLongValue("-nb", 32);
  bool no_blocked = P.getOptionValue("-noblocked");
  bool largemem = P.getOptionValue("-largemem");

  std::cout << "### Application: WidestPath (Single Source Widest-Path)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: -src = " << src << " -nb (num_buckets) = " << num_buckets << std::endl;
  std::cout << "### ------------------------------------" << endl;

  // Define weight type mapping from edge -> custom_weight
  using edge_weight_type = uintE;
  auto degree_f = [&](size_t i) {
    auto vtx = GA.get_vertex(i);
    return std::max(vtx.getInDegree(), vtx.getOutDegree());
  };
  auto degree_im = pbbslib::make_sequence<size_t>(GA.n, degree_f);
  size_t max_degree = pbbslib::reduce_max(degree_im);
  size_t normalize = 2*max_degree+1;
  auto unweighted_weights = [&] (const uintE& u, const uintE& v, const W& wgh) -> uintE {
    uintE deg_u = degree_f(u);
    uintE deg_v = degree_f(v);
    return pbbs::log2_up((size_t)((1/static_cast<double>(deg_u + deg_v + 1))*normalize));
  };
  auto get_weight = gw<G>(unweighted_weights);

  if (num_buckets != (((uintE)1) << pbbslib::log2_up(num_buckets))) {
    std::cout << "Please specify a number of buckets that is a power of two"
              << "\n";
    exit(-1);
  }
  timer t; t.start();
  if (P.getOptionValue("-bf")) {
    auto widths = WidestPathBF(GA, get_weight, src);
  } else {
    auto widths = WidestPath(GA, get_weight, src, num_buckets, largemem, no_blocked);
  }
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

generate_weighted_main(WidestPath_runner, false);
