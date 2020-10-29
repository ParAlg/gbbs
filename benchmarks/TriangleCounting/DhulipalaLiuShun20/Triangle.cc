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
// numactl -i all ./Triangle -rounds 2 -s -c -m clueweb_sym.bytepda
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm

#include "shared.h"
#include "Triangle.h"
#include "makkar.h"


namespace gbbs {
template <class Graph, class UT>
double Dynamic_Triangle_runner(Graph& G, UT& updates, size_t batch_size, commandLine P) {
  // auto ordering = P.getOptionValue("-ordering", "degree");
  std::cout << "### Application: Dynamic Triangle Counting" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  // std::cout << "### n: " << G.n << std::endl;
  // std::cout << "### m: " << G.m << std::endl;
  // std::cout << "### b: " << updates.size() << std::endl;
  // // std::cout << "### Params: ordering=" << ordering << std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  assert(P.getOption("-s"));
  size_t count = 0;
  auto f = [&] (uintE u, uintE v, uintE w) { };
  timer t; t.start();
  count = Dynamic_Triangle(G, updates, f, batch_size, P);
  std::cout << "### Num triangles = " << count << "\n";
  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

#define run_dynamic_app(G, updates, APP, rounds, batch_size)                                            \
  auto before_state = gbbs::get_pcm_state();                               \
  pbbs::timer st;                                                                \
  double total_time = 0.0;                                                 \
  for (size_t r = 0; r < rounds; r++) {                                    \
    total_time += APP(G, updates, batch_size, P);                                               \
  }                                                                        \
  auto time_per_iter = total_time / rounds;                                \
  std::cout << "# time per iter: " << time_per_iter << "\n";               \
  auto after_state = gbbs::get_pcm_state();                                \
  gbbs::print_pcm_stats(before_state, after_state, rounds, time_per_iter); \
  G.del();


template<class Graph>
inline vector<gbbs::gbbs_io::Edge<int>> shuffle_edges(Graph G, int weight){
  vector<gbbs::gbbs_io::Edge<int>> updates_shuffled;
    size_t m = G.num_edges();
    auto perm = pbbs::random_permutation<gbbs::uintE>(m/2);
    pbbs::sequence<std::tuple<gbbs::uintE, gbbs::uintE, typename Graph::weight_type>> edge_list = G.edges();
    auto edge_list_dedup = pbbs::filter(edge_list, [&](const std::tuple<gbbs::uintE, gbbs::uintE, typename Graph::weight_type> & e){
      return std::get<0>(e) < std::get<1>(e);
    });
    edge_list.clear();
    // assert(edge_list_dedup.size() == m/2);
    for (size_t i = 0; i< edge_list_dedup.size(); ++i) {
      updates_shuffled.emplace_back(gbbs::gbbs_io::Edge<int>(std::get<0>(edge_list_dedup[perm[i]]), std::get<1>(edge_list_dedup[perm[i]]), weight));
    }
    edge_list_dedup.clear();
    std::cout << "shuffled and deduped" << std::endl;
    std::cout << updates_shuffled.size() << " deduped edges" << std::endl;
    return updates_shuffled;    
}


template<class Graph>
inline tuple<vector<gbbs::gbbs_io::Edge<int>>, vector<gbbs::gbbs_io::Edge<pbbs::empty>>> shuffle_edges_mix(Graph G, size_t batch_size){
  vector<gbbs::gbbs_io::Edge<int>> updates_shuffled;
  vector<gbbs::gbbs_io::Edge<pbbs::empty>> base_graph;
    size_t m = G.num_edges();
    auto perm = pbbs::random_permutation<gbbs::uintE>(m/2);
    pbbs::sequence<std::tuple<gbbs::uintE, gbbs::uintE, typename Graph::weight_type>> edge_list = G.edges();
    auto edge_list_dedup = pbbs::filter(edge_list, [&](const std::tuple<gbbs::uintE, gbbs::uintE, typename Graph::weight_type> & e){
      return std::get<0>(e) < std::get<1>(e);
    });
    edge_list.clear();
    // assert(edge_list_dedup.size() == m/2);
    size_t batch_num = edge_list_dedup.size() / batch_size;
    size_t end = batch_num * batch_size;
    for (size_t j = 0; j< batch_num; ++j) {
      for (size_t i = j * batch_size/2; i< (j+1)*batch_size/2; ++i) {
        updates_shuffled.emplace_back(gbbs::gbbs_io::Edge<int>(std::get<0>(edge_list_dedup[perm[i]]), std::get<1>(edge_list_dedup[perm[i]]), 1));
      }
      for (size_t i = end/2 + j * batch_size/2; i< end/2 +  (j+1)*batch_size/2; ++i) {
        updates_shuffled.emplace_back(gbbs::gbbs_io::Edge<int>(std::get<0>(edge_list_dedup[perm[i]]), std::get<1>(edge_list_dedup[perm[i]]), 0));
      }
    }
    for (size_t i = end/2; i< end; ++i) {
      base_graph.emplace_back(gbbs::gbbs_io::Edge<pbbs::empty>(std::get<0>(edge_list_dedup[perm[i]]), std::get<1>(edge_list_dedup[perm[i]])));
    }

    edge_list_dedup.clear();
    std::cout << "shuffled and deduped" << std::endl;
    std::cout << updates_shuffled.size() << " deduped edges" << std::endl;
    std::cout << end/2 << " : " <<  end <<std::endl;
    return make_tuple(updates_shuffled, base_graph);    
}

// Read first [partial_end] edges from a file that has the following format:
//     # There can be comments at the top of the file as long as each line of
//     # the comment starts with '#'.
//     <edge 1 first endpoint> <edge 1 second endpoint>
//     <edge 2 first endpoint> <edge 2 second endpoint>
//     <edge 3 first endpoint> <edge 3 second endpoint>
//     ...
//     <edge m first endpoint> <edge m second endpoint>
// give a weighted edge list with weight being "weight" for all edges
template <class weight_type>
std::vector<gbbs::gbbs_io::Edge<weight_type>> DBT_read_edge_list(const char* filename, bool file_weighted, weight_type my_weight, size_t partial_end) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }
  gbbs::gbbs_io::internal::skip_ifstream_comments(&file);

  std::vector<gbbs::gbbs_io::Edge<weight_type>> edge_list;
  gbbs::uintE from;
  gbbs::uintE to;
  size_t line_ct = 0;
  if(file_weighted){
    weight_type weight;
    while (file >> from >> to >> weight) {
      edge_list.emplace_back(from, to, weight);
      line_ct ++;
      if(partial_end  != 0 && line_ct  >= partial_end) break;
    }
  }else{
    while (file >> from >> to) {
      edge_list.emplace_back(from, to, my_weight);
      line_ct ++;
      if(partial_end  != 0 && line_ct  >= partial_end) break;
    }
  }
  return edge_list;
}

/* Macro to generate binary for unweighted graph applications that can ingest
 * only
 * symmetric graph inputs and weighted edge updates input (weight is int type) */
// #define generate_symmetric_dynamic_main(APP, mutates)
  int main(int argc, char* argv[]) {
    gbbs::commandLine P(argc, argv, " [-s] <inFile> <updateFile1>");
    char* iFile = P.getArgument(1);
    char* uFile1 = P.getArgument(0);
    int weight = P.getOptionIntValue("-w", 1);
    size_t batch_size =  P.getOptionLongValue("-batchsize", 10000);
    bool symmetric = P.getOptionValue("-s");
    bool compressed = P.getOptionValue("-c");
    bool mmap = P.getOptionValue("-m");
    bool mmapcopy = false;//mutates;
    bool shuffle = P.getOptionValue("-shuffle");
    bool start_graph = P.getOptionValue("-sg");
    size_t partial_end = P.getOptionLongValue("-partial", 0);
    if(shuffle && start_graph){
      std::cout << "can't shuffle start graph, can use only one of -shuffle and -eg" << std::endl;
      abort();
    }

    if (!symmetric) {
      std::cout << "# The application expects the input graph to be symmetric (-s "
              "flag)."
           << std::endl;
      std::cout << "# Please run on a symmetric input." << std::endl;
    }

    size_t rounds = P.getOptionLongValue("-rounds", 1);
    gbbs::pcm_init();
    vector<gbbs::gbbs_io::Edge<int>> updates;
    if(!shuffle){
      if(weight  == 0){
        updates = DBT_read_edge_list<int>(uFile1, true, 0, partial_end);
        cout << "# update size: " << updates.size() << endl;
      }else if (weight == 1){
        updates = DBT_read_edge_list<int>(uFile1, false, 1, partial_end);
        cout << "# update size: " << updates.size() << endl;
      }else if (weight == 2){\
        updates = DBT_read_edge_list<int>(uFile1, false, 0, partial_end);
        cout << "# update size: " << updates.size() << endl;
      }else{
        std::cout << "# wrong  weighted flag. use 0 for weighted, 1 for inserts , 2 for deletes"  << std::endl;
      }
    }

    bool makkar = P.getOptionValue("-makkar");
    if(makkar){
      gbbs::symmetric_graph<gbbs::symmetric_vertex, pbbslib::empty> G;
      if(shuffle) {
        G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);
        updates = shuffle_edges(G, weight);
      }
      gbbs::alloc_init(G);
      gbbs::Makkar_Dynamic_Triangle(G, updates, batch_size, weight, P);
    }else{


    if (compressed) {
       cout << "not supporting compressed format" << endl;
      //gbbs::symmetric_graph<gbbs::csv_bytepd_amortized, pbbslib::empty> G;
      //if(shuffle || start_graph) G = gbbs::gbbs_io::read_compressed_symmetric_graph<pbbslib::empty>(iFile, mmap, mmapcopy);
      //if(shuffle) {
      //  updates = shuffle_edges(G, weight);
     // }
     // gbbs::alloc_init(G);
      //run_dynamic_app(G, updates, gbbs::Dynamic_Triangle_runner, rounds, batch_size)
    } else {
      gbbs::symmetric_graph<gbbs::symmetric_vertex, pbbslib::empty> G;
      if(shuffle || start_graph) G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap);
      if(shuffle) {
        if(weight == 0){
          vector<gbbs::gbbs_io::Edge<pbbs::empty>> base_graph;
          tie(updates, base_graph) = shuffle_edges_mix(G, batch_size);
          G.del();
          G = gbbs::gbbs_io::edge_list_to_symmetric_graph(base_graph);
        }else{
          updates = shuffle_edges(G, weight);
        }
      }
      gbbs::alloc_init(G);
      run_dynamic_app(G, updates, gbbs::Dynamic_Triangle_runner, rounds, batch_size)
    }

    }// end of else markkar
    gbbs::alloc_finish();
  }

  //
  //   "Usage: ./Triangle [-s] <inFile> <updateFile> \n"
  //   "Methods:\n"
  //   "  -sg: inFile is the original graph"
  //   "  -trict: triangle counts in <inFile>"

  //   "  -shuffle: shuffle edges in <inFile> as updates, do not use <updateFile>"

  //   "  -n: number of vertices"

  //   "  -static: statically count"
  //   "  -partial_end: when reading edge list files, only read first [partial_end] edges"

  //  OPTIONAL:
  //   "  -c: compressed <inFile>"
  //   "  -w: 0 if the edge list is weighted with 32-bit integers., 1 if unweighted inserts. 2 if unweighted deletes\n"
  //   "  -nb: number of batches"
  //   "  -bo: updates start eith [bo]th batch", if [-eg], first [bo] batches are statically counted
  //   "  -blocksize: blocksize to use"

  // if neither -sg or -shuffle is used, <inFile> is not used
  // [-sg]  [-trict] : inFile is the original graph, has trict triangles. start with [bo]th batch
  // [-n] : first [bo] batches of <updateFile> are counted statically, start with [bo]th batch
  // [-shuffle]: first [bo] batches of shuffled <inFile> edges are counted statically, start with [bo]th batch
// generate_symmetric_dynamic_main(gbbs::Dynamic_Triangle_runner, false);
