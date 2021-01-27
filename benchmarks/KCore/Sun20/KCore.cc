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
// ./KCore -s -rounds 1 -o ./out.txt ../../../inputs/star.txt 

#define OUTPUT_RESULT

#include "gbbs/gbbs.h"
#include <unordered_set>
using namespace std;

namespace gbbs {
struct pair_hash
{
    template <class T1, class T2>
    size_t operator () (pair<T1, T2> const &pair) const
    {
        return gbbs::hash_combine(parlay::hash64_2(pair.first), parlay::hash64_2(pair.second));
	}
};

template<class Graph>
inline void output_edges_sun(Graph &G, std::ofstream& out){
	using Edge = pair<uintE, uintE>;
	using Update = tuple<uintE, uintE, bool>;
	using Unordered_map = unordered_map<Edge, size_t, pair_hash>;
	using W = typename Graph::weight_type;
	Unordered_map set;
  
    size_t m = G.num_edges();
	size_t n = G.n;
    auto perm = parlay::random_permutation<uintE>(m);
	auto updates_shuffled = parlay::sequence<Update>(m);

	size_t idx = 0;
	for (size_t i = 0; i < n; i++) {
		auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
			size_t permed_idx = perm[idx];
			if (u < v) { 
				// processed first, insert to table with permuted idx
				set.insert(pair(Edge(u,v), permed_idx));
			}else if(u > v){ 
				// cout << u << " " << v << endl;
				// check table for permuted idx, write the first edge as insert, larger one as delete
				Unordered_map::const_iterator permed_idx_other_itr = set.find(Edge(v,u));
				assert(permed_idx_other_itr != set.end());//must exist in map
				size_t permed_idx_other = permed_idx_other_itr->second;
				size_t small_permed_idx = min(permed_idx_other, permed_idx);
				size_t large_permed_idx = max(permed_idx_other, permed_idx);

				updates_shuffled[small_permed_idx] = Update(v,u,true);// order of edges have to be the same
				updates_shuffled[large_permed_idx] = Update(v,u,false);
			}else{
				// u==v, will be filtered out
				updates_shuffled[permed_idx] = Update(v,u,true);
			}
			idx++;
		};
		G.get_vertex(i).out_neighbors().map(map_f, /* parallel = */ false);
	}
	// 

	// filter out self-edge
	auto updates_cleaned = parlay::filter(updates_shuffled, [&](const Update& e){
      return get<0>(e) != get<1>(e) ;
    });
	// 


	cout << "writing edges to outputfile" << endl;

	for (size_t i = 0; i< updates_cleaned.size(); ++i) {
		gbbs::uintE u = get<0>(updates_cleaned[i]);
		gbbs::uintE v = get<1>(updates_cleaned[i]);
		if(get<2>(updates_cleaned[i])){
			out << "+ " << u << " " << v << " " << i << endl;
		}else{
			out << "- " << u << " " << v << endl;
		}
	}
	cout << "finish writing" << endl;
    cout << updates_cleaned.size() << " non-self edges" << endl;

	perm.clear();
	updates_shuffled.clear();
	updates_cleaned.clear();  
}


template <class Graph>
double KCore_runner(Graph& G, commandLine P) {
  cout << "### Application: KCore Sun20" << endl;
  cout << "### Graph: " << P.getArgument(0) << endl;
  cout << "### Threads: " << num_workers() << endl;
  cout << "### n: " << G.n << endl;
  cout << "### m: " << G.m << endl;
  cout << "### Params: " <<  endl;
  cout << "### ------------------------------------" << endl;
  assert(P.getOption("-s"));
  assert(P.getOption("-rounds") == 1);

//   const string kOutputFlag{"-of"};
  auto outfile = P.getOptionValue("-o", "");
  std::cout << "Outfile: " << outfile << std::endl;
  if (outfile == "") {
    std::cout << "Please specify an output file" << std::endl;
    exit(0);
  }
  std::ofstream out(outfile.c_str());

  timer t; t.start();

  output_edges_sun(G, out);

  double tt = t.stop();

  return tt;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::KCore_runner, false);