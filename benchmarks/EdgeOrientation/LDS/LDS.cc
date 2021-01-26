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

#include "LDS.h"
using namespace std;

namespace gbbs {
template<class Graph>
inline parlay::sequence<tuple<uintE, uintE>> shuffle_edges(Graph &G){
  
    size_t m = G.num_edges();
    auto perm = parlay::random_permutation<uintE>(m/2);
    parlay::sequence<tuple<uintE, uintE, typename Graph::weight_type>> edge_list = G.edges();
    auto edge_list_dedup = parlay::filter(edge_list, [&](const tuple<uintE, uintE, typename Graph::weight_type> & e){
      return get<0>(e) < get<1>(e);
    });
    edge_list.clear();
    // assert(edge_list_dedup.size() == m/2);
    parlay::sequence<tuple<uintE, uintE>> updates_shuffled = parlay::sequence<tuple<uintE, uintE>>(edge_list_dedup.size());
    parallel_for(0, edge_list_dedup.size(), [&](size_t i) {
      updates_shuffled[i] = make_tuple(get<0>(edge_list_dedup[perm[i]]), get<1>(edge_list_dedup[perm[i]]));
    });
    edge_list_dedup.clear();
    cout << "shuffled and deduped" << endl;
    cout << updates_shuffled.size() << " deduped edges" << endl;
    return updates_shuffled;    
}



template <class Graph>
double LDS_runner(Graph& G, commandLine P) {
  cout << "### Application: LDS" << endl;
  cout << "### Graph: " << P.getArgument(0) << endl;
  cout << "### Threads: " << num_workers() << endl;
  cout << "### n: " << G.n << endl;
  cout << "### m: " << G.m << endl;
  cout << "### Params: " <<  endl;
  cout << "### ------------------------------------" << endl;
  assert(P.getOption("-s"));

  const string kInputFlag{"-i"};
  const char* const input_file{P.getOptionValue(kInputFlag)};
  bool is_ordered{P.getOption("-o")};

  using W = typename Graph::weight_type;

  // BatchDynamicEdges<W> batch_edge_list = (input_file && input_file[0]) ? 
  //   read_batch_dynamic_edge_list<W>(input_file, is_ordered) : 
  //   BatchDynamicEdges<W>{};

  parlay::sequence<tuple<uintE, uintE>> batch_edge_list = shuffle_edges(G);

  timer t; t.start();

  RunLDS(G, batch_edge_list);
//  auto cores = (fa) ? LDS_FA(G, num_buckets) : LDS(G, num_buckets);
  double tt = t.stop();

  cout << "### Running Time: " << tt << endl;

  return tt;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::LDS_runner, false);
