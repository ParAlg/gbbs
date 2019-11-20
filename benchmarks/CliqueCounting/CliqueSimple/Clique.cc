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

#include "Clique.h"

template <class Graph>
double Clique_runner(Graph& GA, commandLine P) {
  long k = P.getOptionLongValue("-k", 3);
  std::string order = P.getOptionValue("-order", "goodrichpszona");
  std::cout << "### Application: Clique" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: -k = " << k << std::endl;
  std::cout << "### ------------------------------------" << endl;
  assert(P.getOption("-s"));

  timer t; t.start();
  size_t count;
  if (order == "goodrichpszona") {
    count = KClique<Graph, GoodrichPszonaOrder>(GA, P, k);
  } else if (order == "kcore") {
    count = KClique<Graph, KCoreOrder>(GA, P, k);
  } else if (order == "barenboimelkin") {
    count = KClique<Graph, BarenboimElkinOrder>(GA, P, k);
  }
  double tt = t.stop();
  std::cout << "count: " << count << std::endl;
  std::cout << "### Running Time: " << tt << std::endl;

  return tt;
}

generate_symmetric_main(Clique_runner, false);
